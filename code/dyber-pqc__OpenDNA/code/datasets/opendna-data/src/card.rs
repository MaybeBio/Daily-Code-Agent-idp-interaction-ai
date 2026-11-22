//! Hugging Face dataset card rendering.
//!
//! Cards are generated from the same data the pipeline writes: split sizes
//! and the schema come from the actual tables, and source checksums come from
//! the lockfile. Prose sections are supplied by each dataset module.

use std::fmt::Write as _;

use crate::cache::Lockfile;
use crate::parquet_out::Table;
use crate::sources::{ExcludedSource, Source};
use crate::split::{Split, SplitSummary};

/// License declared in the card's YAML header.
#[derive(Debug, Clone)]
pub enum CardLicense {
    /// A license id from the Hugging Face list, e.g. `cc-by-4.0`.
    Known(&'static str),
    /// A license not on the Hugging Face list (`license: other`).
    Other { name: &'static str, link: &'static str },
}

pub struct CardInput<'a> {
    pub pretty_name: &'static str,
    pub license: CardLicense,
    /// One paragraph explaining what the license of the derived dataset is and why.
    pub license_note: &'static str,
    pub tags: &'static [&'static str],
    pub summary: String,
    pub sources: &'a [Source],
    pub excluded: Vec<&'static ExcludedSource>,
    pub construction: String,
    pub schema: &'a Table,
    /// What a row is, plural ("mutations", "proteins", "domains").
    pub row_unit: &'static str,
    /// What a split group is, plural ("clusters", "components").
    pub group_unit: &'static str,
    pub splits: &'a SplitSummary,
    /// Named tables, with row counts, for a dataset that has no train/test split.
    ///
    /// Some data cannot honestly be split. The antibody cohorts are the case
    /// this exists for: inbred mice share exact rearrangements across animals
    /// and studies, so every mouse donor is one indivisible group and there is
    /// no grouped split that puts mouse rows on both sides. A card for such a
    /// dataset has to describe cohorts and say why there are no splits, rather
    /// than print a split table with everything in `train`.
    ///
    /// When this is non-empty, `splits` is expected to be empty and these
    /// names become the Hugging Face split names.
    pub cohorts: &'a [(&'static str, u64)],
    pub statistics: Vec<(String, String)>,
    pub citations: &'static str,
    pub limitations: &'static [&'static str],
    pub lock: &'a Lockfile,
}

pub fn size_category(n: u64) -> &'static str {
    match n {
        0..1_000 => "n<1K",
        1_000..10_000 => "1K<n<10K",
        10_000..100_000 => "10K<n<100K",
        100_000..1_000_000 => "100K<n<1M",
        1_000_000..10_000_000 => "1M<n<10M",
        _ => "10M<n<100M",
    }
}

pub fn render(c: &CardInput<'_>) -> String {
    let mut s = String::new();
    // YAML front matter.
    s.push_str("---\n");
    match &c.license {
        CardLicense::Known(id) => {
            let _ = writeln!(s, "license: {id}");
        }
        CardLicense::Other { name, link } => {
            let _ = writeln!(s, "license: other\nlicense_name: {name}\nlicense_link: {link}");
        }
    }
    let _ = writeln!(s, "pretty_name: {}", c.pretty_name);
    s.push_str("tags:\n");
    for t in c.tags {
        let _ = writeln!(s, "- {t}");
    }
    let cohort_rows: u64 = c.cohorts.iter().map(|(_, n)| n).sum();
    let _ = writeln!(s, "size_categories:\n- {}", size_category(c.splits.total_rows() + cohort_rows));
    s.push_str("configs:\n- config_name: default\n  data_files:\n");
    for sp in Split::ALL {
        if c.splits.get(sp).rows > 0 {
            let _ = writeln!(s, "  - split: {sp}\n    path: data/{sp}.parquet");
        }
    }
    for (name, rows) in c.cohorts {
        if *rows > 0 {
            let _ = writeln!(s, "  - split: {name}\n    path: data/{name}.parquet");
        }
    }
    s.push_str("---\n\n");

    let _ = writeln!(s, "# {}\n\n{}\n", c.pretty_name, c.summary.trim());

    s.push_str("## Sources and licenses\n\n");
    s.push_str("| Source | License | File | SHA-256 |\n|---|---|---|---|\n");
    for src in c.sources {
        let sha = c.lock.sources.get(src.id).map_or("(not locked)".to_owned(), |e| format!("`{}`", e.sha256));
        let _ = writeln!(
            s,
            "| {} | [{}]({}) | [{}]({}) | {} |",
            src.description,
            src.license.name(),
            src.license.url(),
            src.file_name,
            src.url,
            sha
        );
    }
    let _ = writeln!(s, "\n{}\n", c.license_note.trim());
    if !c.excluded.is_empty() {
        s.push_str("### Sources considered and excluded\n\n");
        for e in &c.excluded {
            let _ = writeln!(s, "- **{}**: {}", e.name, e.reason);
        }
        s.push('\n');
    }

    let _ = writeln!(s, "## Construction\n\n{}\n", c.construction.trim());

    s.push_str("## Schema\n\n| Column | Type | Description |\n|---|---|---|\n");
    for col in &c.schema.columns {
        let _ = writeln!(s, "| `{}` | {} | {} |", col.name, col.data.type_name(), col.description);
    }

    if !c.cohorts.is_empty() {
        let _ = writeln!(s, "\n## Cohorts\n\n| Cohort | {} |\n|---|---:|", c.row_unit);
        for (name, rows) in c.cohorts {
            let _ = writeln!(s, "| {name} | {rows} |");
        }
        let _ = writeln!(s, "| total | {cohort_rows} |\n");
        let _ = writeln!(
            s,
            "There is no train/validation/test split. Why not is under Construction, and it is \
             a fact about the data rather than a choice: read it before splitting these cohorts \
             yourself.\n"
        );
        return finish(s, c);
    }

    // When every row is its own split group the second column repeats the
    // first, which reads as a mistake and hides the fact that makes them
    // equal. One column, and a sentence saying why.
    let one_per_group = Split::ALL.iter().all(|sp| {
        let k = c.splits.get(*sp);
        k.rows == k.groups
    });
    if one_per_group {
        let _ = writeln!(s, "\n## Splits\n\n| Split | {} |\n|---|---:|", c.row_unit);
        for sp in Split::ALL {
            let _ = writeln!(s, "| {sp} | {} |", c.splits.get(sp).rows);
        }
        let _ = writeln!(s, "| total | {} |\n", c.splits.total_rows());
        let _ = writeln!(
            s,
            "Each of these {} is its own split unit, so the no-shared-fold guarantee comes from \
             the grouping described under Construction rather than from grouping rows together \
             here.\n",
            c.row_unit
        );
    } else {
        let _ = writeln!(s, "\n## Splits\n\n| Split | {} | {} |\n|---|---:|---:|", c.row_unit, c.group_unit);
        for sp in Split::ALL {
            let k = c.splits.get(sp);
            let _ = writeln!(s, "| {sp} | {} | {} |", k.rows, k.groups);
        }
        let _ = writeln!(s, "| total | {} | |\n", c.splits.total_rows());
    }

    finish(s, c)
}

/// The sections every card ends with, whether it describes splits or cohorts.
fn finish(mut s: String, c: &CardInput<'_>) -> String {
    use std::fmt::Write as _;
    if !c.statistics.is_empty() {
        s.push_str("### Cleaning statistics\n\n| Step | Count |\n|---|---:|\n");
        for (k, v) in &c.statistics {
            let _ = writeln!(s, "| {k} | {v} |");
        }
        s.push('\n');
    }

    s.push_str("## Known limitations\n\n");
    for l in c.limitations {
        let _ = writeln!(s, "- {l}");
    }
    let _ = writeln!(
        s,
        "\n## Citation\n\nIf you use this dataset, cite the upstream sources:\n\n```bibtex\n{}\n```\n",
        c.citations.trim()
    );
    s.push_str(
        "## Reproduction\n\nBuilt with `opendna-data` from the OpenDNA repository \
         (`datasets/`). The exact upstream files are pinned by SHA-256 in \
         `datasets/sources.lock`; `manifest.json` next to this card records the \
         pipeline parameters and the checksum of every Parquet file.\n",
    );
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parquet_out::ColumnData;
    use crate::split::SplitCounts;

    #[test]
    fn renders_front_matter_and_sections() {
        let mut t = Table::new();
        t.push("id", "identifier", ColumnData::Utf8(vec![]));
        let splits = SplitSummary {
            train: SplitCounts { rows: 8, groups: 4 },
            validation: SplitCounts { rows: 1, groups: 1 },
            test: SplitCounts::default(),
        };
        let lock = Lockfile::default();
        let card = render(&CardInput {
            pretty_name: "Example",
            license: CardLicense::Other { name: "cc-by-sa-2.1-jp", link: "https://example.invalid" },
            license_note: "note",
            tags: &["biology"],
            summary: "Summary.".into(),
            sources: &[],
            excluded: vec![],
            construction: "Steps.".into(),
            schema: &t,
            row_unit: "rows",
            group_unit: "clusters",
            splits: &splits,
            cohorts: &[],
            statistics: vec![("dropped".into(), "3".into())],
            citations: "@misc{x}",
            limitations: &["small"],
            lock: &lock,
        });
        assert!(card.starts_with("---\nlicense: other\nlicense_name: cc-by-sa-2.1-jp\n"));
        assert!(card.contains("- split: train\n    path: data/train.parquet"));
        // Empty splits are not declared, so the HF loader does not look for missing files.
        assert!(!card.contains("path: data/test.parquet"));
        assert!(card.contains("| `id` | string | identifier |"));
        assert!(card.contains("| train | 8 | 4 |"));
        assert!(card.contains("- n<1K"));
    }

    /// A dataset whose split unit is its row prints one column, not two of
    /// the same.
    ///
    /// Burial is that case: every sampled domain is its own split unit, so the
    /// card read "| Split | domains | domains |" with identical figures, which
    /// looks like a mistake and hides the fact that makes them identical. The
    /// no-shared-fold guarantee there comes from the grouping its Construction
    /// section describes, which is why the card has to say so instead.
    #[test]
    fn a_dataset_grouped_one_row_at_a_time_prints_one_column() {
        let mut t = Table::new();
        t.push("id", "identifier", ColumnData::Utf8(vec![]));
        let lock = Lockfile::default();
        let input = |splits: &SplitSummary| {
            render(&CardInput {
                pretty_name: "Example",
                license: CardLicense::Known("cc-by-4.0"),
                license_note: "note",
                tags: &["biology"],
                summary: "Summary.".into(),
                sources: &[],
                excluded: vec![],
                construction: "Steps.".into(),
                schema: &t,
                row_unit: "domains",
                group_unit: "domains",
                splits,
                cohorts: &[],
                statistics: vec![],
                citations: "@misc{x}",
                limitations: &["small"],
                lock: &lock,
            })
        };

        let one_each = SplitSummary {
            train: SplitCounts { rows: 751, groups: 751 },
            validation: SplitCounts { rows: 95, groups: 95 },
            test: SplitCounts { rows: 93, groups: 93 },
        };
        let card = input(&one_each);
        assert!(card.contains("| Split | domains |\n|---|---:|"), "{card}");
        assert!(card.contains("| train | 751 |"), "{card}");
        assert!(!card.contains("| train | 751 | 751 |"), "the column is printed twice:\n{card}");
        assert!(card.contains("is its own split unit"), "{card}");

        // And a dataset that really does group rows keeps both columns.
        let grouped = SplitSummary {
            train: SplitCounts { rows: 314_716, groups: 123 },
            validation: SplitCounts { rows: 37_983, groups: 18 },
            test: SplitCounts { rows: 36_369, groups: 24 },
        };
        let card = input(&grouped);
        assert!(card.contains("| train | 314716 | 123 |"), "{card}");
        assert!(!card.contains("is its own split unit"), "{card}");
    }

    #[test]
    fn size_buckets() {
        assert_eq!(size_category(999), "n<1K");
        assert_eq!(size_category(1_000), "1K<n<10K");
        assert_eq!(size_category(400_000), "100K<n<1M");
    }
}
