//! A small columnar table type and its Parquet writer.
//!
//! Each dataset converts its records into a [`Table`]. The same column
//! definitions drive the Parquet schema (including per-field descriptions as
//! Arrow field metadata) and the schema section of the dataset card, so the
//! documentation cannot drift from the files.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use anyhow::{bail, Context, Result};
use arrow::array::{
    ArrayRef, BooleanArray, Float32Builder, Float64Array, Int32Array, Int64Array, Int8Builder, ListBuilder,
    StringArray,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use serde::Serialize;

use crate::cache::hash_file;

#[derive(Debug, Clone)]
pub enum ColumnData {
    Utf8(Vec<String>),
    OptUtf8(Vec<Option<String>>),
    Int32(Vec<i32>),
    OptInt64(Vec<Option<i64>>),
    Float64(Vec<f64>),
    OptFloat64(Vec<Option<f64>>),
    Bool(Vec<bool>),
    ListInt8(Vec<Vec<i8>>),
    /// Per-residue floats. f32 rather than f64: these are measured
    /// quantities with a few significant figures, and halving the file costs
    /// nothing real.
    ListFloat32(Vec<Vec<f32>>),
}

impl ColumnData {
    fn len(&self) -> usize {
        match self {
            ColumnData::Utf8(v) => v.len(),
            ColumnData::OptUtf8(v) => v.len(),
            ColumnData::Int32(v) => v.len(),
            ColumnData::OptInt64(v) => v.len(),
            ColumnData::Float64(v) => v.len(),
            ColumnData::OptFloat64(v) => v.len(),
            ColumnData::Bool(v) => v.len(),
            ColumnData::ListInt8(v) => v.len(),
            ColumnData::ListFloat32(v) => v.len(),
        }
    }

    fn data_type(&self) -> DataType {
        match self {
            ColumnData::Utf8(_) | ColumnData::OptUtf8(_) => DataType::Utf8,
            ColumnData::Int32(_) => DataType::Int32,
            ColumnData::OptInt64(_) => DataType::Int64,
            ColumnData::Float64(_) | ColumnData::OptFloat64(_) => DataType::Float64,
            ColumnData::Bool(_) => DataType::Boolean,
            ColumnData::ListInt8(_) => DataType::List(Arc::new(Field::new("item", DataType::Int8, true))),
            ColumnData::ListFloat32(_) => {
                DataType::List(Arc::new(Field::new("item", DataType::Float32, true)))
            }
        }
    }

    fn nullable(&self) -> bool {
        matches!(self, ColumnData::OptUtf8(_) | ColumnData::OptInt64(_) | ColumnData::OptFloat64(_))
    }

    /// Type name used in dataset cards.
    pub fn type_name(&self) -> &'static str {
        match self {
            ColumnData::Utf8(_) => "string",
            ColumnData::OptUtf8(_) => "string (nullable)",
            ColumnData::Int32(_) => "int32",
            ColumnData::OptInt64(_) => "int64 (nullable)",
            ColumnData::Float64(_) => "float64",
            ColumnData::OptFloat64(_) => "float64 (nullable)",
            ColumnData::Bool(_) => "bool",
            ColumnData::ListInt8(_) => "list<int8>",
            ColumnData::ListFloat32(_) => "list<float32>",
        }
    }

    fn to_array(&self) -> ArrayRef {
        match self {
            ColumnData::Utf8(v) => Arc::new(StringArray::from_iter_values(v)),
            ColumnData::OptUtf8(v) => Arc::new(StringArray::from(v.clone())),
            ColumnData::Int32(v) => Arc::new(Int32Array::from(v.clone())),
            ColumnData::OptInt64(v) => Arc::new(Int64Array::from(v.clone())),
            ColumnData::Float64(v) => Arc::new(Float64Array::from(v.clone())),
            ColumnData::OptFloat64(v) => Arc::new(Float64Array::from(v.clone())),
            ColumnData::Bool(v) => Arc::new(BooleanArray::from(v.clone())),
            ColumnData::ListInt8(v) => {
                let mut b = ListBuilder::new(Int8Builder::new());
                for row in v {
                    b.values().append_slice(row);
                    b.append(true);
                }
                Arc::new(b.finish())
            }
            ColumnData::ListFloat32(v) => {
                let mut b = ListBuilder::new(Float32Builder::new());
                for row in v {
                    b.values().append_slice(row);
                    b.append(true);
                }
                Arc::new(b.finish())
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct Column {
    pub name: &'static str,
    pub description: &'static str,
    pub data: ColumnData,
}

#[derive(Debug, Clone, Default)]
pub struct Table {
    pub columns: Vec<Column>,
}

impl Table {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, name: &'static str, description: &'static str, data: ColumnData) {
        self.columns.push(Column { name, description, data });
    }

    pub fn n_rows(&self) -> usize {
        self.columns.first().map_or(0, |c| c.data.len())
    }

    fn to_batch(&self) -> Result<RecordBatch> {
        let n = self.n_rows();
        for c in &self.columns {
            if c.data.len() != n {
                bail!("column {} has {} rows, expected {n}", c.name, c.data.len());
            }
        }
        let fields: Vec<Field> = self
            .columns
            .iter()
            .map(|c| {
                Field::new(c.name, c.data.data_type(), c.data.nullable())
                    .with_metadata(HashMap::from([("description".to_owned(), c.description.to_owned())]))
            })
            .collect();
        let arrays: Vec<ArrayRef> = self.columns.iter().map(|c| c.data.to_array()).collect();
        Ok(RecordBatch::try_new(Arc::new(Schema::new(fields)), arrays)?)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct WrittenFile {
    pub path: String,
    pub rows: u64,
    pub bytes: u64,
    pub sha256: String,
}

/// Writes `table` to `path` as zstd-compressed Parquet and returns its checksum.
pub fn write_parquet(table: &Table, path: &Path, display_path: &str) -> Result<WrittenFile> {
    let batch = table.to_batch()?;
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::try_new(9)?))
        .set_max_row_group_row_count(Some(50_000))
        .build();
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    w.write(&batch)?;
    w.close()?;
    let (sha256, bytes) = hash_file(path)?;
    Ok(WrittenFile { path: display_path.to_owned(), rows: table.n_rows() as u64, bytes, sha256 })
}

#[cfg(test)]
mod tests {
    use super::*;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    fn sample() -> Table {
        let mut t = Table::new();
        t.push("id", "row id", ColumnData::Utf8(vec!["a".into(), "b".into()]));
        t.push("x", "a value", ColumnData::Float64(vec![1.5, -2.0]));
        t.push("n", "maybe", ColumnData::OptInt64(vec![Some(3), None]));
        t.push("labels", "per residue", ColumnData::ListInt8(vec![vec![0, 1, 1], vec![]]));
        t
    }

    #[test]
    fn writes_readable_parquet_with_descriptions() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("t.parquet");
        let info = write_parquet(&sample(), &p, "t.parquet").unwrap();
        assert_eq!(info.rows, 2);
        let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(&p).unwrap()).unwrap();
        let schema = reader.schema().clone();
        assert_eq!(schema.field_with_name("x").unwrap().metadata()["description"], "a value");
        let batches: Vec<_> = reader.build().unwrap().collect::<Result<_, _>>().unwrap();
        assert_eq!(batches.iter().map(|b| b.num_rows()).sum::<usize>(), 2);
    }

    #[test]
    fn output_is_byte_identical_across_runs() {
        let dir = tempfile::tempdir().unwrap();
        let a = write_parquet(&sample(), &dir.path().join("a.parquet"), "a").unwrap();
        let b = write_parquet(&sample(), &dir.path().join("b.parquet"), "b").unwrap();
        assert_eq!(a.sha256, b.sha256);
    }

    #[test]
    fn ragged_table_is_rejected() {
        let mut t = sample();
        t.push("bad", "", ColumnData::Int32(vec![1]));
        let dir = tempfile::tempdir().unwrap();
        assert!(write_parquet(&t, &dir.path().join("x.parquet"), "x").is_err());
    }
}
