//! Picking a long run up where it stopped.
//!
//! A run over fifty thousand sequences takes hours, and the things that
//! interrupt it are not exotic: a scheduler's wall-clock limit, a power cut, a
//! laptop lid, somebody needing the GPU. Without a checkpoint the only option
//! is to start again, which on a queue with a four-hour limit and a five-hour
//! job is not an option at all.
//!
//! So a run can write its results as it goes and be told to continue. The hard
//! part is not the file; it is making sure that continuing produces the same
//! answer as never having stopped.
//!
//! Two things make that true.
//!
//! The first is that work is done in fixed chunks, recorded in the checkpoint.
//! Embedding batches are composed from the sequences in a chunk, and batch
//! composition changes the order of floating-point sums, so a run cut in two
//! at an arbitrary point would not agree with an uninterrupted one in the last
//! digit of a few values. Cutting at the same boundaries every time removes
//! the question: a run resumed from any chunk is byte-identical to one that
//! was never interrupted, provided the chunk size matches. That is why
//! `--chunk` is recorded and why resuming with a different one is refused
//! rather than accommodated.
//!
//! The second is that everything else which could change a number — the input
//! bytes, the head's build, the embedding model, the device, the batch budget,
//! the mutation list, the confidence cut-off — is recorded in the header and
//! compared on resume. A mismatch names the field and stops. Half a table
//! produced one way and half another is the mixture that [`crate::merge`]
//! exists to refuse, and it would be odd to refuse it between files while
//! creating it inside one.
//!
//! The file is line-delimited JSON, appended to and flushed to disk after
//! every chunk. That shape is chosen for the failure that actually happens: a
//! process killed partway through a write leaves a truncated last line, which
//! is detectable precisely because the lines are independent. A reader keeps
//! the records up to the last one that parses and the writer truncates the
//! file to that point before appending, so a torn tail costs one chunk rather
//! than the run.

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::io::{BufRead as _, Seek as _, Write as _};
use std::path::Path;

/// The checkpoint format's version.
///
/// Bumped if the shape changes. A checkpoint from another version is refused
/// rather than read optimistically: the whole purpose of the file is to be
/// trusted about what has already been computed.
pub const FORMAT: u32 = 1;

/// Everything about a run that has to match for its checkpoint to be resumed.
///
/// The rule for what belongs here is the same as [`crate::merge`]'s: anything
/// that changes a number, plus anything that changes the shape of the table.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Header {
    /// The checkpoint format's version.
    pub format: u32,
    /// Which property is being predicted.
    pub task: String,
    /// SHA-256 of the input's bytes.
    ///
    /// The bytes, not the path: the same filename holding a different file is
    /// exactly the mistake worth catching, and an input that arrived on
    /// standard input has no path to compare.
    pub input_sha256: String,
    /// How many units of work the run has. Sequences, or for stability the
    /// sequences that have mutations requested against them.
    pub units: usize,
    /// Units per chunk.
    pub chunk: usize,
    /// The head's build id.
    pub head_build: Option<String>,
    /// The head's recorded weights digest.
    pub head_weights_sha256: Option<String>,
    /// The embedding model's repository.
    pub embedding_model: String,
    /// The device, as it is described to a reader.
    pub device: String,
    /// Padded tokens per embedding batch.
    pub batch_tokens: usize,
    /// SHA-256 of the mutation file's bytes, for a stability run.
    pub mutations_sha256: Option<String>,
    /// The confidence cut-off, which adds a column.
    pub min_confidence: Option<f64>,
}

impl Header {
    /// Which fields differ from `other`, in words a user can act on.
    ///
    /// Returns the differences rather than a boolean because "this checkpoint
    /// does not match this run" is not a useful thing to be told. Which field
    /// differs is the difference between a typo in a filename and a head that
    /// was re-downloaded halfway through a week of work.
    fn differences(&self, other: &Self) -> Vec<String> {
        let mut out = Vec::new();
        let mut text = |what: &str, a: &str, b: &str| {
            if a != b {
                out.push(format!("{what}: the checkpoint says {a}, this run has {b}"));
            }
        };
        text("format version", &self.format.to_string(), &other.format.to_string());
        text("task", &self.task, &other.task);
        text("input digest", &self.input_sha256, &other.input_sha256);
        text("units of work", &self.units.to_string(), &other.units.to_string());
        text("chunk size", &self.chunk.to_string(), &other.chunk.to_string());
        text("head build", &shown(self.head_build.as_deref()), &shown(other.head_build.as_deref()));
        text(
            "head weights digest",
            &shown(self.head_weights_sha256.as_deref()),
            &shown(other.head_weights_sha256.as_deref()),
        );
        text("embedding model", &self.embedding_model, &other.embedding_model);
        text("device", &self.device, &other.device);
        text("batch tokens", &self.batch_tokens.to_string(), &other.batch_tokens.to_string());
        text(
            "mutation file digest",
            &shown(self.mutations_sha256.as_deref()),
            &shown(other.mutations_sha256.as_deref()),
        );
        // Compared by bits rather than by value, so that two runs are either
        // the same run or reported as different. No cut-off and a cut-off of
        // zero are different runs and produce different columns.
        let cut = |v: Option<f64>| v.map_or_else(|| "none".to_string(), |f| format!("{f}"));
        if self.min_confidence.map(f64::to_bits) != other.min_confidence.map(f64::to_bits) {
            out.push(format!(
                "confidence cut-off: the checkpoint says {}, this run has {}",
                cut(self.min_confidence),
                cut(other.min_confidence)
            ));
        }
        out
    }
}

/// An optional field as the difference report states it.
fn shown(value: Option<&str>) -> String {
    value.map_or_else(|| "nothing recorded".to_string(), str::to_owned)
}

/// One chunk's results.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Record {
    /// Which chunk, counting from zero.
    pub chunk: usize,
    /// The table's columns.
    ///
    /// Repeated on every record rather than stated once in the header,
    /// because the header is written before any work is done and the columns
    /// are a property of the result. Carrying them per chunk also makes the
    /// file self-checking: chunks whose columns disagree cannot be assembled
    /// into one table, and this is where that is noticed.
    pub columns: Vec<String>,
    /// The rows, as the table holds them.
    pub rows: Vec<Vec<String>>,
}

/// How many chunks a run has.
///
/// # Errors
/// If the chunk size is zero, which would be an infinite number of chunks.
pub fn chunks(units: usize, chunk: usize) -> Result<usize> {
    if chunk == 0 {
        bail!("--chunk must be at least 1");
    }
    Ok(units.div_ceil(chunk))
}

/// What a checkpoint file holds, and how much of it is usable.
struct Existing {
    /// The header as recorded.
    header: Header,
    /// The records that are contiguous from chunk zero.
    records: Vec<Record>,
    /// How many bytes of the file those records occupy, including the header
    /// line. Anything after this is a torn or out-of-order tail.
    good_bytes: u64,
    /// Whether anything was discarded.
    discarded: bool,
}

/// Read a checkpoint, keeping what is usable and saying what was not.
///
/// A line that does not parse ends the usable part. So does a chunk index that
/// is not the next one expected: records are appended in order by a single
/// process, so a gap means the file was concatenated or edited, and guessing
/// which rows belong where is not something a results file should do.
fn load(path: &Path) -> Result<Existing> {
    let file = std::fs::File::open(path).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = std::io::BufReader::new(file).lines();
    let first =
        lines.next().transpose().with_context(|| format!("reading {}", path.display()))?.ok_or_else(
            || anyhow::anyhow!("{} is empty. Delete it and start the run again.", path.display()),
        )?;
    let header: Header = serde_json::from_str(&first).with_context(|| {
        format!(
            "parsing the first line of {} as a checkpoint header. If this is not an OpenDNA \
             checkpoint, point --checkpoint somewhere else.",
            path.display()
        )
    })?;

    // +1 for the newline the reader strips. Byte offsets, because the file is
    // truncated to this point before anything is appended.
    let mut good_bytes = first.len() as u64 + 1;
    let mut records: Vec<Record> = Vec::new();
    let mut discarded = false;
    for line in lines {
        let Ok(text) = line else {
            discarded = true;
            break;
        };
        let Ok(record) = serde_json::from_str::<Record>(&text) else {
            discarded = true;
            break;
        };
        if record.chunk != records.len() {
            discarded = true;
            break;
        }
        good_bytes += text.len() as u64 + 1;
        records.push(record);
    }
    Ok(Existing { header, records, good_bytes, discarded })
}

/// A checkpoint open for appending.
#[derive(Debug)]
pub struct Writer {
    file: std::fs::File,
    /// The chunks already recorded, in order.
    records: Vec<Record>,
}

impl Writer {
    /// The chunks already done, in order from chunk zero.
    pub fn records(&self) -> &[Record] {
        &self.records
    }

    /// Record a chunk and get it onto the disk.
    ///
    /// Flushed and synced before returning, because a checkpoint that is still
    /// in a buffer when the power goes is not a checkpoint. The cost is one
    /// sync per chunk, which against the minutes a chunk of embedding takes is
    /// not measurable.
    ///
    /// # Errors
    /// If the record cannot be written or the file cannot be synced, or if the
    /// chunk's columns disagree with the chunks already recorded.
    pub fn append(&mut self, record: Record) -> Result<()> {
        if let Some(first) = self.records.first() {
            if first.columns != record.columns {
                bail!(
                    "chunk {} produced columns {:?} where the checkpoint's earlier chunks have \
                     {:?}. These cannot be one table.",
                    record.chunk,
                    record.columns,
                    first.columns
                );
            }
        }
        let line = serde_json::to_string(&record).context("encoding a checkpoint record")?;
        writeln!(self.file, "{line}").context("appending to the checkpoint")?;
        self.file.flush().context("flushing the checkpoint")?;
        self.file.sync_all().context("syncing the checkpoint to disk")?;
        self.records.push(record);
        Ok(())
    }
}

/// Open a checkpoint for this run, resuming an existing one where it matches.
///
/// A file that is not there is created with `header` as its first line. A file
/// that is there is read, compared against `header`, and truncated to its last
/// usable record so that appending cannot leave a torn line in the middle.
///
/// # Errors
/// If the file cannot be read or written, if it is not a checkpoint, or if it
/// describes a different run.
pub fn open(path: &Path, header: &Header) -> Result<Writer> {
    if !path.exists() {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
            }
        }
        let mut file = std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let line = serde_json::to_string(header).context("encoding the checkpoint header")?;
        writeln!(file, "{line}").with_context(|| format!("writing {}", path.display()))?;
        file.sync_all().with_context(|| format!("syncing {}", path.display()))?;
        return Ok(Writer { file, records: Vec::new() });
    }

    let Existing { header: found, records, good_bytes, discarded } = load(path)?;
    let differences = found.differences(header);
    if !differences.is_empty() {
        bail!(
            "{} is a checkpoint for a different run, so continuing it would produce a table \
             that is partly one run and partly another:\n  {}\n\
             Either put the earlier run's arguments back, or point --checkpoint at a new file.",
            path.display(),
            differences.join("\n  ")
        );
    }
    if discarded {
        eprintln!(
            "The end of {} was incomplete, which is what an interrupted write looks like. \
             Dropping it and continuing from chunk {}.",
            path.display(),
            records.len()
        );
    }

    // Opened for writing and positioned by hand rather than opened in append
    // mode. A torn tail has to be cut off before anything is added after it,
    // and on Windows a handle opened for appending is not allowed to shorten
    // the file it is appending to.
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .open(path)
        .with_context(|| format!("opening {} to append", path.display()))?;
    // Only after the header matched: truncating a file that turns out to
    // describe another run would destroy somebody else's checkpoint.
    file.set_len(good_bytes).with_context(|| format!("truncating {}", path.display()))?;
    file.seek(std::io::SeekFrom::Start(good_bytes))
        .with_context(|| format!("seeking to the end of {}", path.display()))?;
    Ok(Writer { file, records })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::fmt::Write as _;

    fn header() -> Header {
        Header {
            format: FORMAT,
            task: "disorder".into(),
            input_sha256: "ab".repeat(32),
            units: 10,
            chunk: 4,
            head_build: Some("f7d99273e4b9b961".into()),
            head_weights_sha256: Some("cd".repeat(32)),
            embedding_model: "facebook/esm2_t30_150M_UR50D".into(),
            device: "CPU".into(),
            batch_tokens: 8192,
            mutations_sha256: None,
            min_confidence: None,
        }
    }

    fn record(chunk: usize) -> Record {
        Record {
            chunk,
            columns: vec!["name".into(), "position".into()],
            rows: vec![vec![format!("seq{chunk}"), "1".into()]],
        }
    }

    fn dir(name: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("opendna-checkpoint-{name}"));
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn the_number_of_chunks_covers_every_unit() {
        assert_eq!(chunks(10, 4).unwrap(), 3);
        assert_eq!(chunks(8, 4).unwrap(), 2);
        assert_eq!(chunks(1, 4).unwrap(), 1);
        assert_eq!(chunks(0, 4).unwrap(), 0);
        assert!(chunks(10, 0).is_err());
    }

    #[test]
    fn a_fresh_checkpoint_starts_empty_and_a_reopened_one_remembers() {
        let path = dir("fresh").join("run.jsonl");
        let mut w = open(&path, &header()).unwrap();
        assert!(w.records().is_empty());
        w.append(record(0)).unwrap();
        w.append(record(1)).unwrap();
        drop(w);

        let again = open(&path, &header()).unwrap();
        assert_eq!(again.records().len(), 2);
        assert_eq!(again.records()[1].rows[0][0], "seq1");
    }

    /// The failure this file format exists for: a process killed partway
    /// through writing a record. The completed chunks must survive.
    #[test]
    fn a_torn_last_line_costs_one_chunk_and_not_the_run() {
        let path = dir("torn").join("run.jsonl");
        let mut w = open(&path, &header()).unwrap();
        w.append(record(0)).unwrap();
        w.append(record(1)).unwrap();
        drop(w);

        // Half a line, as an interrupted write leaves behind.
        let mut text = std::fs::read_to_string(&path).unwrap();
        text.push_str("{\"chunk\":2,\"columns\":[\"na");
        std::fs::write(&path, &text).unwrap();

        let mut w = open(&path, &header()).unwrap();
        assert_eq!(w.records().len(), 2);
        // Appending after the truncation must leave a file that reads back
        // cleanly, rather than one with rubbish in the middle.
        w.append(record(2)).unwrap();
        drop(w);
        let reread = open(&path, &header()).unwrap();
        assert_eq!(reread.records().len(), 3);
        assert_eq!(reread.records()[2].rows[0][0], "seq2");
    }

    #[test]
    fn a_checkpoint_for_another_input_is_refused_and_says_which_field() {
        let path = dir("other-input").join("run.jsonl");
        let mut w = open(&path, &header()).unwrap();
        w.append(record(0)).unwrap();
        drop(w);

        let mut mine = header();
        mine.input_sha256 = "ef".repeat(32);
        let error = format!("{:#}", open(&path, &mine).unwrap_err());
        assert!(error.contains("input digest"), "{error}");
        assert!(error.contains("partly one run and partly another"), "{error}");
    }

    /// Chunk size changes batch composition, which changes the last digit of
    /// some values, which is exactly the mixture this refuses.
    #[test]
    fn a_checkpoint_at_another_chunk_size_is_refused() {
        let path = dir("other-chunk").join("run.jsonl");
        drop(open(&path, &header()).unwrap());
        let mut mine = header();
        mine.chunk = 8;
        let error = format!("{:#}", open(&path, &mine).unwrap_err());
        assert!(error.contains("chunk size"), "{error}");
    }

    #[test]
    fn a_checkpoint_from_another_device_or_batch_size_is_refused() {
        let path = dir("other-device").join("run.jsonl");
        drop(open(&path, &header()).unwrap());
        let mut mine = header();
        mine.device = "CUDA GPU".into();
        mine.batch_tokens = 0;
        let error = format!("{:#}", open(&path, &mine).unwrap_err());
        assert!(error.contains("device"), "{error}");
        assert!(error.contains("batch tokens"), "{error}");
    }

    /// A confidence cut-off adds a column, so a checkpoint written without one
    /// cannot be continued by a run that has one.
    #[test]
    fn a_confidence_cut_off_is_part_of_the_run_s_identity() {
        let path = dir("cut").join("run.jsonl");
        drop(open(&path, &header()).unwrap());
        let mut mine = header();
        mine.min_confidence = Some(0.1);
        let error = format!("{:#}", open(&path, &mine).unwrap_err());
        assert!(error.contains("confidence cut-off"), "{error}");
    }

    #[test]
    fn a_file_that_is_not_a_checkpoint_is_refused_rather_than_overwritten() {
        let path = dir("not-a-checkpoint").join("notes.txt");
        std::fs::write(&path, "these are my notes\n").unwrap();
        let error = format!("{:#}", open(&path, &header()).unwrap_err());
        assert!(error.contains("checkpoint header"), "{error}");
        // The file is still there, which matters: a mistyped --checkpoint
        // should not destroy what it was pointed at.
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "these are my notes\n");
    }

    #[test]
    fn chunks_that_disagree_about_their_columns_cannot_be_one_table() {
        let path = dir("columns").join("run.jsonl");
        let mut w = open(&path, &header()).unwrap();
        w.append(record(0)).unwrap();
        let mut odd = record(1);
        odd.columns.push("low_confidence".into());
        let error = format!("{:#}", w.append(odd).unwrap_err());
        assert!(error.contains("cannot be one table"), "{error}");
    }

    /// A gap in the chunk indices means the file was edited or concatenated.
    /// Keeping the prefix is safe; guessing where the rest belongs is not.
    #[test]
    fn a_gap_in_the_chunks_ends_the_usable_part() {
        let path = dir("gap").join("run.jsonl");
        let mut w = open(&path, &header()).unwrap();
        w.append(record(0)).unwrap();
        drop(w);
        let mut text = std::fs::read_to_string(&path).unwrap();
        let _ = writeln!(text, "{}", serde_json::to_string(&record(7)).unwrap());
        std::fs::write(&path, text).unwrap();
        assert_eq!(open(&path, &header()).unwrap().records().len(), 1);
    }
}
