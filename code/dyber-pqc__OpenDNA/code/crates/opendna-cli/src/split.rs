//! Cutting a FASTA into pieces a scheduler can run in parallel.
//!
//! A pipeline needs chunks, and a lab that staged one binary to run offline
//! should not also need `seqkit` to produce them. Splitting here means the
//! chunks are parsed by exactly the code that will predict on them: a sequence
//! `predict` accepts is a sequence `split` preserves, including the
//! selenocysteine and ambiguity codes that a stricter tool would reject or
//! silently rewrite.
//!
//! Two things matter beyond correctness.
//!
//! **Determinism.** The same input and the same part count give the same
//! chunks, byte for byte, on any machine. A pipeline that resumes after a
//! failure has to produce the same chunk boundaries or it is resuming into a
//! different job.
//!
//! **Balance.** Runtime scales with residues, not with sequences, so a part
//! holding ten 1000-residue proteins takes far longer than one holding ten
//! 70-residue proteins. Cutting by sequence count leaves a scheduler waiting
//! on one straggler while the rest of the cluster idles, so parts are balanced
//! by residue count instead, keeping the input order within each.

use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};

/// How to cut the file.
#[derive(Debug, Clone, Copy)]
pub enum Cut {
    /// Into this many parts, balanced by residue count.
    Parts(usize),
    /// Into parts of at most this many sequences each.
    PerChunk(usize),
}

/// Decide which sequences go in which part.
///
/// Returns one index range per part. Parts are contiguous in input order, so a
/// chunk is a slice of the original file rather than a scatter of it, which
/// keeps a chunk readable next to the input it came from.
///
/// Balancing by residues is greedy and deliberately simple: walk the
/// sequences, close the current part once it has reached its share of the
/// total, and never leave fewer parts than asked for by overfilling early.
fn plan(lengths: &[usize], cut: Cut) -> Vec<(usize, usize)> {
    if lengths.is_empty() {
        return Vec::new();
    }
    match cut {
        Cut::PerChunk(size) => {
            let size = size.max(1);
            (0..lengths.len()).step_by(size).map(|s| (s, (s + size).min(lengths.len()))).collect()
        }
        Cut::Parts(parts) => {
            let parts = parts.max(1).min(lengths.len());
            let total: usize = lengths.iter().sum();
            let mut out = Vec::with_capacity(parts);
            let mut start = 0;
            let mut carried = 0usize;
            for part in 0..parts {
                if start >= lengths.len() {
                    break;
                }
                // The share this part should carry, measured against what is
                // left rather than the original total, so an early part that
                // overshot does not push the error all the way to the end.
                let remaining_parts = parts - part;
                if remaining_parts == 1 {
                    out.push((start, lengths.len()));
                    break;
                }
                let remaining: usize = total - carried;
                let target = remaining / remaining_parts;
                let mut end = start;
                let mut taken = 0usize;
                // Always take at least one, and leave at least one sequence
                // for each part still to come.
                let latest = lengths.len() - (remaining_parts - 1);
                while end < latest && (taken < target || end == start) {
                    taken += lengths[end];
                    end += 1;
                }
                out.push((start, end));
                carried += taken;
                start = end;
            }
            out
        }
    }
}

/// Write the parts of `input` into `out_dir`.
///
/// # Errors
/// If the input cannot be read or parsed, or a part cannot be written.
pub fn run(input: &Path, out_dir: &Path, cut: Cut) -> Result<Vec<PathBuf>> {
    let text = std::fs::read_to_string(input).with_context(|| format!("reading {}", input.display()))?;
    let proteins = opendna_io::fasta::read_raw(&text)
        .with_context(|| format!("parsing {} as FASTA", input.display()))?;
    if proteins.is_empty() {
        bail!("{} contains no sequences", input.display());
    }
    std::fs::create_dir_all(out_dir).with_context(|| format!("creating {}", out_dir.display()))?;

    let lengths: Vec<usize> = proteins.iter().map(|(_, s)| s.chars().count()).collect();
    let parts = plan(&lengths, cut);
    let mut written = Vec::with_capacity(parts.len());

    for (index, (start, end)) in parts.iter().enumerate() {
        // Zero-padded and fixed width, so the parts sort the same way in a
        // shell glob as they do numerically.
        let path = out_dir.join(format!("part-{index:04}.fasta"));
        let mut body = String::new();
        for (name, sequence) in &proteins[*start..*end] {
            body.push('>');
            body.push_str(name);
            body.push('\n');
            body.push_str(sequence);
            body.push('\n');
        }
        std::fs::write(&path, body).with_context(|| format!("writing {}", path.display()))?;
        let residues: usize = lengths[*start..*end].iter().sum();
        println!("{}\t{} sequences\t{residues} residues", path.display(), end - start);
        written.push(path);
    }
    eprintln!("Wrote {} parts into {}", written.len(), out_dir.display());
    Ok(written)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn spans(lengths: &[usize], cut: Cut) -> Vec<usize> {
        plan(lengths, cut).iter().map(|(s, e)| lengths[*s..*e].iter().sum()).collect()
    }

    /// Every sequence lands in exactly one part, which is the property a
    /// pipeline silently depends on: a dropped chunk is a result table
    /// missing rows nobody counted.
    #[test]
    fn the_parts_cover_every_sequence_exactly_once() {
        let lengths: Vec<usize> = (0..37).map(|i| 50 + (i * 73) % 500).collect();
        for cut in [Cut::Parts(1), Cut::Parts(4), Cut::Parts(37), Cut::Parts(100), Cut::PerChunk(5)] {
            let parts = plan(&lengths, cut);
            let mut covered: Vec<usize> = Vec::new();
            for (start, end) in &parts {
                assert!(start < end, "{cut:?} produced an empty part");
                covered.extend(*start..*end);
            }
            assert_eq!(covered, (0..lengths.len()).collect::<Vec<_>>(), "{cut:?}");
        }
    }

    /// The reason for balancing at all: a scheduler waits for the slowest
    /// part, so cutting by count would leave one part doing most of the work.
    #[test]
    fn parts_are_balanced_by_residues_not_by_count() {
        // Ten short sequences then two very long ones.
        let mut lengths = vec![50; 10];
        lengths.extend([2000, 2000]);
        let residues = spans(&lengths, Cut::Parts(2));
        assert_eq!(residues.len(), 2);
        let (low, high) = (residues[0].min(residues[1]), residues[0].max(residues[1]));
        assert!(high <= low * 2, "parts should be within a factor of two by residues, got {residues:?}");

        // Cutting by count would put all the work in one part, which is the
        // behaviour being avoided.
        let by_count = spans(&lengths, Cut::PerChunk(6));
        assert!(by_count[1] > by_count[0] * 5, "{by_count:?}");
    }

    #[test]
    fn asking_for_more_parts_than_sequences_gives_one_each() {
        let parts = plan(&[10, 20, 30], Cut::Parts(99));
        assert_eq!(parts, vec![(0, 1), (1, 2), (2, 3)]);
    }

    #[test]
    fn a_single_part_is_the_whole_file() {
        assert_eq!(plan(&[10, 20, 30], Cut::Parts(1)), vec![(0, 3)]);
        assert_eq!(plan(&[10, 20, 30], Cut::PerChunk(99)), vec![(0, 3)]);
    }

    #[test]
    fn an_empty_input_plans_nothing() {
        assert_eq!(plan(&[], Cut::Parts(4)), [] as [(usize, usize); 0]);
        assert_eq!(plan(&[], Cut::PerChunk(4)), [] as [(usize, usize); 0]);
    }

    /// A pipeline that resumes has to produce the same boundaries, or it is
    /// resuming into a different job.
    #[test]
    fn the_plan_is_deterministic() {
        let lengths: Vec<usize> = (0..61).map(|i| 40 + (i * 157) % 900).collect();
        for cut in [Cut::Parts(7), Cut::PerChunk(9)] {
            assert_eq!(plan(&lengths, cut), plan(&lengths, cut), "{cut:?}");
        }
    }
}
