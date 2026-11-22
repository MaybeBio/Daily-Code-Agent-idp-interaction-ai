//! An on-disk store of per-sequence matrices.
//!
//! Embedding a dataset takes minutes to hours, so the results are cached and
//! reused across training runs. One cache holds matrices of a fixed width:
//! per-residue embeddings (`length × hidden_size`) or masked log-probability
//! rows (`length × vocab_size`).
//!
//! Values are stored as `f16`. The heads consume them as `f32`; half the
//! precision is ample for an input feature and halves both file size and the
//! time spent reading them back.
//!
//! Entries are keyed by a hash of the sequence, so duplicate sequences are
//! embedded once, and a run interrupted part-way can be resumed.

use anyhow::{bail, Context, Result};
use half::f16;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

/// Content hash of a sequence, used as the cache key.
pub fn key_for(sequence: &str) -> String {
    use std::fmt::Write;
    let digest = Sha256::digest(sequence.as_bytes());
    digest[..12].iter().fold(String::with_capacity(24), |mut s, b| {
        let _ = write!(s, "{b:02x}");
        s
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Entry {
    key: String,
    /// Offset into the data file, in values (not bytes).
    offset: u64,
    rows: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct Index {
    /// Model that produced these matrices; a mismatch invalidates the cache.
    model: String,
    /// Width of every matrix.
    cols: usize,
    entries: Vec<Entry>,
}

/// A readable and appendable matrix cache.
pub struct MatrixCache {
    index_path: PathBuf,
    data_path: PathBuf,
    index: Index,
    lookup: HashMap<String, usize>,
    next_offset: u64,
}

impl MatrixCache {
    /// Open an existing cache or create an empty one.
    ///
    /// # Errors
    /// If the files cannot be read, or an existing cache was built with a
    /// different model or width.
    pub fn open(dir: &Path, name: &str, model: &str, cols: usize) -> Result<Self> {
        std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
        let index_path = dir.join(format!("{name}.json"));
        let data_path = dir.join(format!("{name}.bin"));

        let index = if index_path.exists() {
            let text = std::fs::read_to_string(&index_path)
                .with_context(|| format!("reading {}", index_path.display()))?;
            let index: Index =
                serde_json::from_str(&text).with_context(|| format!("parsing {}", index_path.display()))?;
            if index.model != model {
                bail!(
                    "{} was built with {} but {model} was requested; delete it to rebuild",
                    index_path.display(),
                    index.model
                );
            }
            if index.cols != cols {
                bail!(
                    "{} holds {}-wide matrices but {cols} was requested; delete it to rebuild",
                    index_path.display(),
                    index.cols
                );
            }
            index
        } else {
            Index { model: model.to_string(), cols, entries: Vec::new() }
        };

        let lookup = index.entries.iter().enumerate().map(|(i, e)| (e.key.clone(), i)).collect();
        let next_offset = index.entries.last().map_or(0, |e| e.offset + (e.rows * cols) as u64);
        Ok(Self { index_path, data_path, index, lookup, next_offset })
    }

    /// Number of cached matrices.
    pub fn len(&self) -> usize {
        self.index.entries.len()
    }

    /// Whether the cache holds no matrices.
    #[allow(dead_code, reason = "clippy requires is_empty wherever len exists")]
    pub fn is_empty(&self) -> bool {
        self.index.entries.is_empty()
    }

    /// Whether `sequence` has already been embedded.
    pub fn contains(&self, sequence: &str) -> bool {
        self.lookup.contains_key(&key_for(sequence))
    }

    /// Append a matrix. Appending a sequence that is already present is a
    /// no-op, so an interrupted run can simply be restarted.
    ///
    /// # Errors
    /// If `values` does not match `rows × cols`, or the write fails.
    pub fn append(&mut self, sequence: &str, rows: usize, values: &[f32]) -> Result<()> {
        if values.len() != rows * self.index.cols {
            bail!("expected {} values for {rows} rows, got {}", rows * self.index.cols, values.len());
        }
        let key = key_for(sequence);
        if self.lookup.contains_key(&key) {
            return Ok(());
        }
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.data_path)
            .with_context(|| format!("opening {}", self.data_path.display()))?;
        let mut writer = BufWriter::new(&mut file);
        for v in values {
            writer.write_all(&f16::from_f32(*v).to_le_bytes())?;
        }
        writer.flush()?;

        self.lookup.insert(key.clone(), self.index.entries.len());
        self.index.entries.push(Entry { key, offset: self.next_offset, rows });
        self.next_offset += values.len() as u64;
        Ok(())
    }

    /// Write the index to disk. Call periodically during a long run so an
    /// interruption loses only the work since the last save.
    ///
    /// # Errors
    /// If the index cannot be written.
    pub fn save(&self) -> Result<()> {
        let text = serde_json::to_string(&self.index)?;
        std::fs::write(&self.index_path, text)
            .with_context(|| format!("writing {}", self.index_path.display()))
    }

    /// Read one matrix back as `f32`, with its row count.
    ///
    /// # Errors
    /// If the sequence is absent or the data file cannot be read.
    pub fn get(&self, sequence: &str) -> Result<(Vec<f32>, usize)> {
        let key = key_for(sequence);
        let index = *self
            .lookup
            .get(&key)
            .with_context(|| format!("sequence {key} is not in {}", self.data_path.display()))?;
        let entry = &self.index.entries[index];
        let count = entry.rows * self.index.cols;
        let mut file =
            File::open(&self.data_path).with_context(|| format!("opening {}", self.data_path.display()))?;
        file.seek(SeekFrom::Start(entry.offset * 2))?;
        let mut bytes = vec![0u8; count * 2];
        file.read_exact(&mut bytes)
            .with_context(|| format!("reading {count} values for {key} from {}", self.data_path.display()))?;
        let values =
            bytes.as_chunks::<2>().0.iter().map(|c| f16::from_le_bytes([c[0], c[1]]).to_f32()).collect();
        Ok((values, entry.rows))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("opendna-cache-test-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    #[test]
    fn round_trips_through_f16() {
        let dir = temp_dir("roundtrip");
        let mut cache = MatrixCache::open(&dir, "test", "model-a", 3).unwrap();
        assert!(cache.is_empty());
        cache.append("MKT", 2, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        cache.save().unwrap();

        let (values, rows) = cache.get("MKT").unwrap();
        assert_eq!(rows, 2);
        assert_eq!(values, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert!(cache.contains("MKT"));
        assert!(!cache.contains("ACDE"));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn reopening_sees_previous_entries() {
        let dir = temp_dir("reopen");
        {
            let mut cache = MatrixCache::open(&dir, "test", "model-a", 2).unwrap();
            cache.append("AA", 1, &[1.0, 2.0]).unwrap();
            cache.append("BB", 2, &[3.0, 4.0, 5.0, 6.0]).unwrap();
            cache.save().unwrap();
        }
        let cache = MatrixCache::open(&dir, "test", "model-a", 2).unwrap();
        assert_eq!(cache.len(), 2);
        assert_eq!(cache.get("BB").unwrap().0, vec![3.0, 4.0, 5.0, 6.0]);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn rejects_a_cache_from_a_different_model_or_width() {
        let dir = temp_dir("mismatch");
        MatrixCache::open(&dir, "test", "model-a", 2).unwrap().save().unwrap();
        assert!(MatrixCache::open(&dir, "test", "model-b", 2).is_err());
        assert!(MatrixCache::open(&dir, "test", "model-a", 4).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn appending_twice_keeps_the_first_copy() {
        let dir = temp_dir("dedup");
        let mut cache = MatrixCache::open(&dir, "test", "m", 1).unwrap();
        cache.append("AA", 1, &[1.0]).unwrap();
        cache.append("AA", 1, &[9.0]).unwrap();
        assert_eq!(cache.len(), 1);
        assert_eq!(cache.get("AA").unwrap().0, vec![1.0]);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn wrong_value_count_is_rejected() {
        let dir = temp_dir("badlen");
        let mut cache = MatrixCache::open(&dir, "test", "m", 3).unwrap();
        assert!(cache.append("AA", 2, &[1.0, 2.0]).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
