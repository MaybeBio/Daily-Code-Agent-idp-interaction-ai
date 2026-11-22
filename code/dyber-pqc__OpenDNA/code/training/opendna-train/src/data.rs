//! Reading the Parquet datasets produced by the `datasets/` workspace.

use anyhow::{anyhow, bail, Context, Result};
use arrow::array::{
    Array, BooleanArray, Float32Array, Float64Array, Int16Array, Int32Array, Int64Array, Int8Array,
    ListArray, StringArray, UInt8Array,
};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::path::Path;

/// Which portion of a dataset to read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Split {
    /// Training rows.
    Train,
    /// Model selection and early stopping.
    Validation,
    /// Held out; read once, for the final numbers.
    Test,
}

impl Split {
    /// File stem used by the dataset builder.
    pub fn name(self) -> &'static str {
        match self {
            Self::Train => "train",
            Self::Validation => "validation",
            Self::Test => "test",
        }
    }

    /// All three, in order.
    pub fn all() -> [Split; 3] {
        [Split::Train, Split::Validation, Split::Test]
    }
}

fn read_batches(root: &Path, dataset: &str, split: Split) -> Result<Vec<RecordBatch>> {
    let path = root.join(dataset).join("data").join(format!("{}.parquet", split.name()));
    let file = File::open(&path).with_context(|| {
        format!("opening {}; build the datasets first (see datasets/README.md)", path.display())
    })?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .with_context(|| format!("reading {}", path.display()))?
        .build()?;
    reader.collect::<Result<Vec<_>, _>>().with_context(|| format!("decoding {}", path.display()))
}

fn strings(batch: &RecordBatch, column: &str) -> Result<Vec<String>> {
    let array = batch.column_by_name(column).ok_or_else(|| anyhow!("missing column {column}"))?;
    let array = array
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| anyhow!("column {column} is not a string column"))?;
    (0..array.len())
        .map(|i| {
            if array.is_null(i) {
                bail!("column {column} has a null at row {i}");
            }
            Ok(array.value(i).to_string())
        })
        .collect()
}

fn floats(batch: &RecordBatch, column: &str) -> Result<Vec<f64>> {
    let array = batch.column_by_name(column).ok_or_else(|| anyhow!("missing column {column}"))?;
    let array = array
        .as_any()
        .downcast_ref::<Float64Array>()
        .ok_or_else(|| anyhow!("column {column} is not a float column"))?;
    (0..array.len())
        .map(|i| {
            if array.is_null(i) {
                bail!("column {column} has a null at row {i}");
            }
            Ok(array.value(i))
        })
        .collect()
}

fn integers(batch: &RecordBatch, column: &str) -> Result<Vec<i64>> {
    let array = batch.column_by_name(column).ok_or_else(|| anyhow!("missing column {column}"))?;
    if let Some(a) = array.as_any().downcast_ref::<Int64Array>() {
        return Ok((0..a.len()).map(|i| a.value(i)).collect());
    }
    if let Some(a) = array.as_any().downcast_ref::<Int32Array>() {
        return Ok((0..a.len()).map(|i| i64::from(a.value(i))).collect());
    }
    bail!("column {column} is not an integer column")
}

/// One point mutation with its measured ΔΔG.
#[derive(Debug, Clone)]
pub struct StabilityRow {
    /// Wild-type sequence the mutation is applied to.
    pub wt_sequence: String,
    /// 0-based position of the substitution.
    pub position: usize,
    /// Wild-type residue.
    pub wild_type: char,
    /// Substituted residue.
    pub mutant: char,
    /// Measured change in folding free energy.
    pub ddg: f64,
}

/// Read the stability dataset.
///
/// # Errors
/// If the file is missing, malformed, or a mutation does not match its
/// wild-type sequence.
pub fn stability(root: &Path, split: Split) -> Result<Vec<StabilityRow>> {
    let mut out = Vec::new();
    for batch in read_batches(root, "stability", split)? {
        let sequences = strings(&batch, "wt_sequence")?;
        // The dataset stores 1-based positions.
        let positions = integers(&batch, "position")?;
        let wild = strings(&batch, "wt_aa")?;
        let mutant = strings(&batch, "mut_aa")?;
        let ddg = floats(&batch, "ddg")?;
        for i in 0..batch.num_rows() {
            let position = usize::try_from(positions[i] - 1)
                .with_context(|| format!("row {i}: position {} is out of range", positions[i]))?;
            let wt = wild[i].chars().next().ok_or_else(|| anyhow!("row {i}: empty wt_aa"))?;
            let mt = mutant[i].chars().next().ok_or_else(|| anyhow!("row {i}: empty mut_aa"))?;
            // Guard against an off-by-one in the position convention, which
            // would silently train the head on the wrong residue.
            let actual = sequences[i].as_bytes().get(position).map(|b| *b as char);
            if actual != Some(wt) {
                bail!("row {i}: sequence has {actual:?} at position {position} but the mutation says {wt:?}");
            }
            out.push(StabilityRow {
                wt_sequence: sequences[i].clone(),
                position,
                wild_type: wt,
                mutant: mt,
                ddg: ddg[i],
            });
        }
    }
    Ok(out)
}

/// One protein with its measured solubility.
#[derive(Debug, Clone)]
pub struct SolubilityRow {
    /// Amino-acid sequence.
    pub sequence: String,
    /// Solubility percentage as reported by the source assay.
    pub solubility: f64,
}

/// Read the solubility dataset.
///
/// # Errors
/// If the file is missing or malformed.
pub fn solubility(root: &Path, split: Split) -> Result<Vec<SolubilityRow>> {
    let mut out = Vec::new();
    for batch in read_batches(root, "solubility", split)? {
        let sequences = strings(&batch, "sequence")?;
        let values = floats(&batch, "solubility")?;
        for (sequence, solubility) in sequences.into_iter().zip(values) {
            out.push(SolubilityRow { sequence, solubility });
        }
    }
    Ok(out)
}

/// One protein with a per-residue disorder annotation.
#[derive(Debug, Clone)]
pub struct DisorderRow {
    /// Amino-acid sequence.
    pub sequence: String,
    /// One label per residue: `true` where the residue is disordered.
    pub disorder: Vec<bool>,
}

/// One domain's per-residue solvent accessibility.
pub struct BurialRow {
    /// Amino-acid sequence.
    pub sequence: String,
    /// Relative solvent accessibility per residue, 0 buried to about 1
    /// fully exposed. Can exceed 1: the Tien et al. maxima are theoretical.
    pub rsa: Vec<f32>,
}

/// Read the burial dataset.
///
/// # Errors
/// If the file is missing or a row's label count does not match its sequence.
pub fn burial(root: &Path, split: Split) -> Result<Vec<BurialRow>> {
    let mut out = Vec::new();
    for batch in read_batches(root, "burial", split)? {
        let sequences = strings(&batch, "sequence")?;
        let lists = batch
            .column_by_name("rsa")
            .ok_or_else(|| anyhow!("missing column rsa"))?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| anyhow!("column rsa is not a list column"))?
            .clone();
        for (i, sequence) in sequences.into_iter().enumerate() {
            let row = lists.value(i);
            let values = accessibility(row.as_ref(), i)?;
            if values.len() != sequence.len() {
                bail!("row {i}: {} values for {} residues", values.len(), sequence.len());
            }
            out.push(BurialRow { sequence, rsa: values });
        }
    }
    Ok(out)
}

/// Read one domain's accessibility values.
///
/// The builder writes `list<float32>`, but float width is an encoding choice
/// a writer is free to make, so a `float64` column means the same thing.
/// Accepting both keeps the reader from breaking on a detail that carries no
/// meaning — the same reasoning as [`labels`] below, which was widened after
/// exactly that happened.
fn accessibility(values: &dyn Array, row: usize) -> Result<Vec<f32>> {
    if let Some(a) = values.as_any().downcast_ref::<Float32Array>() {
        return Ok((0..a.len()).map(|j| a.value(j)).collect());
    }
    if let Some(a) = values.as_any().downcast_ref::<Float64Array>() {
        return Ok((0..a.len()).map(|j| a.value(j) as f32).collect());
    }
    Err(anyhow!("row {row}: accessibility values are {:?}, which is not a float list", values.data_type()))
}

/// Read one protein's per-residue labels as booleans.
///
/// The dataset builder writes `list<int8>`, but the width of an integer
/// column is an encoding choice a writer is free to change, and a boolean
/// column would represent the same thing just as faithfully. Accepting any of
/// those spellings keeps the reader from breaking on a detail that carries no
/// meaning. Anything else is a genuine mismatch and is reported with the type
/// that was actually found, rather than a bare "not an integer list".
fn labels(values: &dyn Array, row: usize) -> Result<Vec<bool>> {
    macro_rules! nonzero {
        ($ty:ty) => {
            values
                .as_any()
                .downcast_ref::<$ty>()
                .map(|a| (0..a.len()).map(|j| a.value(j) != 0).collect::<Vec<bool>>())
        };
    }

    if let Some(a) = values.as_any().downcast_ref::<BooleanArray>() {
        return Ok((0..a.len()).map(|j| a.value(j)).collect());
    }
    nonzero!(Int8Array)
        .or_else(|| nonzero!(Int16Array))
        .or_else(|| nonzero!(Int32Array))
        .or_else(|| nonzero!(Int64Array))
        .or_else(|| nonzero!(UInt8Array))
        .ok_or_else(|| {
            anyhow!(
                "row {row}: disorder labels are {:?}, which is neither an integer nor a boolean \
                 list",
                values.data_type()
            )
        })
}

/// Read the disorder dataset.
///
/// # Errors
/// If the file is missing, malformed, or a label array does not match the
/// length of its sequence.
pub fn disorder(root: &Path, split: Split) -> Result<Vec<DisorderRow>> {
    let mut out = Vec::new();
    for batch in read_batches(root, "disorder", split)? {
        let sequences = strings(&batch, "sequence")?;
        let lists = batch
            .column_by_name("disorder")
            .ok_or_else(|| anyhow!("missing column disorder"))?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| anyhow!("column disorder is not a list column"))?
            .clone();
        for (i, sequence) in sequences.into_iter().enumerate() {
            let row = lists.value(i);
            let values = labels(row.as_ref(), i)?;
            if values.len() != sequence.len() {
                bail!("row {i}: {} labels for {} residues", values.len(), sequence.len());
            }
            out.push(DisorderRow { sequence, disorder: values });
        }
    }
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn labels_are_read_at_any_integer_width() {
        let expected = vec![false, true, true, false];
        let raw = [0_i64, 1, 1, 0];
        assert_eq!(labels(&Int8Array::from(vec![0_i8, 1, 1, 0]), 0).unwrap(), expected);
        assert_eq!(labels(&Int16Array::from(vec![0_i16, 1, 1, 0]), 0).unwrap(), expected);
        assert_eq!(labels(&Int32Array::from(vec![0_i32, 1, 1, 0]), 0).unwrap(), expected);
        assert_eq!(labels(&Int64Array::from(raw.to_vec()), 0).unwrap(), expected);
        assert_eq!(labels(&UInt8Array::from(vec![0_u8, 1, 1, 0]), 0).unwrap(), expected);
        assert_eq!(labels(&BooleanArray::from(expected.clone()), 0).unwrap(), expected);
    }

    #[test]
    fn a_non_integer_label_column_is_rejected_and_names_the_type() {
        let err = labels(&Float64Array::from(vec![0.0, 1.0]), 7).unwrap_err().to_string();
        assert!(err.contains("row 7"), "{err}");
        assert!(err.contains("Float64"), "{err}");
    }

    #[test]
    fn any_nonzero_label_counts_as_disordered() {
        // DisProt writes 0/1, but a 2 would be a label, not a no.
        assert_eq!(labels(&Int8Array::from(vec![0_i8, 2, -1]), 0).unwrap(), vec![false, true, true]);
    }
}
