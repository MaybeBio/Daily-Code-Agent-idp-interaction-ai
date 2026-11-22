//! Minimal FASTA reading and writing.

use std::io::{BufRead, Write};

use anyhow::{bail, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Header line without the leading `>`.
    pub header: String,
    /// Sequence lines concatenated, whitespace removed.
    pub sequence: String,
}

impl Record {
    /// First whitespace-delimited token of the header.
    pub fn id(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }
}

/// Parses a multi-line FASTA stream.
pub fn read(r: impl BufRead) -> Result<Vec<Record>> {
    let mut out = Vec::new();
    let mut cur: Option<Record> = None;
    for (lineno, line) in r.lines().enumerate() {
        let line = line?;
        let line = line.trim_end();
        if let Some(h) = line.strip_prefix('>') {
            if let Some(rec) = cur.take() {
                out.push(rec);
            }
            cur = Some(Record { header: h.trim().to_owned(), sequence: String::new() });
        } else if !line.is_empty() {
            match cur.as_mut() {
                Some(rec) => rec.sequence.extend(line.chars().filter(|c| !c.is_whitespace())),
                None => bail!("FASTA line {} has sequence data before any header", lineno + 1),
            }
        }
    }
    if let Some(rec) = cur {
        out.push(rec);
    }
    Ok(out)
}

/// Writes records with 80-column sequence lines.
pub fn write<'a>(mut w: impl Write, records: impl IntoIterator<Item = (&'a str, &'a str)>) -> Result<()> {
    for (id, seq) in records {
        writeln!(w, ">{id}")?;
        for chunk in seq.as_bytes().chunks(80) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multiline_records() {
        let text = ">a desc\nMKV\nLLA\n\n>b\nGG\n";
        let recs = read(text.as_bytes()).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id(), "a");
        assert_eq!(recs[0].sequence, "MKVLLA");
        assert_eq!(recs[1].sequence, "GG");
    }

    #[test]
    fn rejects_headerless_data() {
        assert!(read("MKV\n>a\nG\n".as_bytes()).is_err());
    }

    #[test]
    fn write_then_read_roundtrips() {
        let long = "A".repeat(170);
        let mut buf = Vec::new();
        write(&mut buf, [("x", long.as_str()), ("y", "MK")]).unwrap();
        let recs = read(buf.as_slice()).unwrap();
        assert_eq!(recs[0].sequence, long);
        assert_eq!(recs[1].id(), "y");
    }
}
