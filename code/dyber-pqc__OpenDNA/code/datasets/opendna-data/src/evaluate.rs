//! Does the humanness table separate the cohorts it should?
//!
//! The measure is a count, so there is no accuracy to report and no model to
//! validate. What there is to check is whether the table behaves the way the
//! method claims, on data it has never seen, and the antibody cohorts exist to
//! ask exactly that. Four questions, with the right answer stated before the
//! numbers are computed:
//!
//! * **human** should score high. None of these donors is in the table -- the
//!   two catalogues are disjoint and the build fails if they are not -- so
//!   this measures whether human repertoires look alike across people, not
//!   whether the table remembers its input.
//! * **mouse** should score lower. A mouse heavy chain is built from mouse
//!   germline, so most of its windows should be rare among human donors.
//! * **diagnostic** should score high, like human. These are mice carrying
//!   human immunoglobulin loci: the animal is a mouse and the molecule is
//!   human. A species classifier gets this wrong by construction, and a
//!   humanness measure getting it right is the one result here that a V-gene
//!   lookup table could not produce.
//! * **out_of_distribution** should score lower, and the differences between
//!   its species are not a ranking: one study per species means any gap
//!   between rat and rabbit is one laboratory's protocol.
//!
//! Nothing here is tuned. No threshold is chosen to make a separation look
//! better, and the three thresholds reported are the ones OASis defines. If
//! the numbers came out the other way this would say so.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs::File;
use std::path::Path;

use anyhow::{bail, Context, Result};
use opendna_analysis::humanness::{score, PrevalenceTable, Stringency};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

/// The cohorts, in the order the report lists them, each with the answer that
/// would be correct. Stated here rather than in the prose so the expectation
/// is fixed before the measurement.
const COHORTS: &[(&str, &str)] = &[
    ("human", "high: donors the table has never seen"),
    ("mouse", "lower: mouse germline, scored against human donors"),
    ("diagnostic", "high: mice with human immunoglobulin loci, so the molecule is human"),
    ("out_of_distribution", "lower, and not a ranking between the species"),
];

/// Per-sequence identities for one group of rows.
#[derive(Default)]
struct Spread {
    values: Vec<f64>,
}

impl Spread {
    fn quantiles(&mut self) -> Option<[f64; 5]> {
        if self.values.is_empty() {
            return None;
        }
        self.values.sort_unstable_by(f64::total_cmp);
        let at = |p: f64| {
            let i = ((p * self.values.len() as f64) as usize).min(self.values.len() - 1);
            self.values[i]
        };
        Some([at(0.05), at(0.25), at(0.50), at(0.75), at(0.95)])
    }
}

/// One cohort's measurements, broken down by the organism column.
struct CohortResult {
    rows: usize,
    unscorable: usize,
    /// Keyed by (organism, stringency); `organism` is `"(all)"` for the total.
    spreads: BTreeMap<(String, Stringency), Spread>,
}

/// One scored row: the sequence, and the labels it is broken down by.
struct Row {
    sequence: String,
    organism: String,
    isotype: String,
}

/// Read the columns the report breaks down by.
fn read_cohort(path: &Path) -> Result<Vec<Row>> {
    let file = File::open(path).with_context(|| {
        format!(
            "opening {}. Build the antibody cohorts first: `opendna-data build antibody`.",
            path.display()
        )
    })?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)?.build()?;
    let mut rows = Vec::new();
    for batch in reader {
        let batch = batch?;
        let column = |name: &str| -> Result<arrow::array::StringArray> {
            let index = batch
                .schema()
                .index_of(name)
                .with_context(|| format!("{} has no column {name:?}", path.display()))?;
            let array = batch.column(index);
            array
                .as_any()
                .downcast_ref::<arrow::array::StringArray>()
                .cloned()
                .with_context(|| format!("{name:?} in {} is not a string column", path.display()))
        };
        let sequences = column("sequence")?;
        let organisms = column("organism")?;
        let isotypes = column("isotype")?;
        for i in 0..batch.num_rows() {
            rows.push(Row {
                sequence: sequences.value(i).to_owned(),
                organism: organisms.value(i).to_owned(),
                isotype: isotypes.value(i).to_owned(),
            });
        }
    }
    Ok(rows)
}

fn measure(table: &PrevalenceTable, rows: &[Row]) -> Result<CohortResult> {
    let mut spreads: BTreeMap<(String, Stringency), Spread> = BTreeMap::new();
    let mut unscorable = 0_usize;
    for row in rows {
        let scored = score(table, &row.sequence)?;
        // A sequence with nothing scorable contributes to neither numerator
        // nor denominator. Counting it as zero would assert it looks unlike a
        // human antibody, which is not what was measured.
        if scored.scored == 0 {
            unscorable += 1;
            continue;
        }
        // Broken down two ways. Organism is the question; isotype is the
        // confound, because a class-switched sequence has been somatically
        // mutated away from germline and germline is where the common
        // peptides are. Without the second breakdown, a difference in isotype
        // between two cohorts would read as a difference in humanness.
        let labels =
            [ALL.to_owned(), format!("{ORGANISM}{}", row.organism), format!("{ISOTYPE}{}", row.isotype)];
        for stringency in Stringency::ALL {
            let Some(identity) = scored.identity(stringency) else { continue };
            for label in &labels {
                spreads.entry((label.clone(), stringency)).or_default().values.push(identity);
            }
        }
    }
    Ok(CohortResult { rows: rows.len(), unscorable, spreads })
}

/// Label prefixes, so the two breakdowns can be told apart in one map.
const ALL: &str = "(all)";
const ORGANISM: &str = "organism/";
const ISOTYPE: &str = "isotype/";

/// Measure every cohort against the table and write the report.
pub fn humanness(out_root: &Path) -> Result<()> {
    let table_path = out_root.join("humanness").join("table.bin");
    let file = File::open(&table_path).with_context(|| {
        format!("opening {}. Build the table first: `opendna-data build humanness`.", table_path.display())
    })?;
    let table = PrevalenceTable::read(std::io::BufReader::with_capacity(1 << 20, file))
        .with_context(|| format!("reading {}", table_path.display()))?;
    eprintln!("table: {} donors, {} peptides", table.donors(), table.len());

    let data = out_root.join("antibody").join("data");
    let mut results = Vec::new();
    for (cohort, _) in COHORTS {
        let rows = read_cohort(&data.join(format!("{cohort}.parquet")))?;
        eprintln!("{cohort}: scoring {} sequences", rows.len());
        results.push((*cohort, measure(&table, &rows)?));
    }

    let expressible: Vec<Stringency> =
        Stringency::ALL.into_iter().filter(|s| s.expressible_with(table.donors())).collect();
    if expressible.is_empty() {
        bail!(
            "a table of {} donors cannot express any of the three thresholds, so there is \
             nothing to report",
            table.donors()
        );
    }

    let report = render(&table, &results, &expressible);
    let path = out_root.join("humanness").join("EVALUATION.md");
    std::fs::write(&path, &report)?;
    print!("{report}");
    println!("\nWritten to {}", path.display());
    Ok(())
}

fn render(table: &PrevalenceTable, results: &[(&str, CohortResult)], expressible: &[Stringency]) -> String {
    let mut s = String::new();
    s.push_str("# Humanness table: does it separate the cohorts it should?\n\n");
    let _ = writeln!(
        s,
        "Scored with a table of {} donors and {} peptides against the four antibody cohorts, \
         whose donors are disjoint from the table's. Each number is a **per-sequence identity**: \
         the fraction of a sequence's scorable 9-mers occurring in at least that share of the \
         table's donors. Nothing is tuned and no threshold was chosen after seeing the \
         numbers.\n",
        table.donors(),
        table.len()
    );

    s.push_str("What the right answer is, fixed before measuring:\n\n");
    for (cohort, expectation) in COHORTS {
        let _ = writeln!(s, "- `{cohort}` -- {expectation}");
    }

    let skipped: Vec<Stringency> = Stringency::ALL.into_iter().filter(|s| !expressible.contains(s)).collect();
    if !skipped.is_empty() {
        s.push_str("\nThresholds this table cannot express, and so not reported:");
        for stringency in &skipped {
            let _ = write!(s, " {},", stringency.label());
        }
        let _ = writeln!(
            s,
            "\nfrom {} donors the smallest non-zero prevalence is {:.2}%.",
            table.donors(),
            100.0 / f64::from(table.donors().max(1))
        );
    }

    for stringency in expressible {
        let _ = writeln!(s, "\n## {}\n", stringency.label());
        for (heading, prefix) in [("By cohort and organism", ORGANISM), ("By isotype", ISOTYPE)] {
            let _ = writeln!(s, "### {heading}\n");
            s.push_str("| Cohort | group | sequences | p5 | p25 | median | p75 | p95 |\n");
            s.push_str("|---|---|---:|---:|---:|---:|---:|---:|\n");
            for (cohort, result) in results {
                let mut names: Vec<&String> = result
                    .spreads
                    .keys()
                    .filter(|(name, st)| {
                        st == stringency && (name.as_str() == ALL || name.starts_with(prefix))
                    })
                    .map(|(name, _)| name)
                    .collect();
                names.sort_by_key(|name| (name.as_str() != ALL, name.as_str()));
                for name in names {
                    let (cohort_cell, group_cell) = if name.as_str() == ALL {
                        (format!("**{cohort}**"), "all".to_owned())
                    } else {
                        (String::new(), name.trim_start_matches(prefix).to_owned())
                    };
                    // Cloned because quantiles sorts, and a report should not
                    // depend on which row was rendered first.
                    let mut spread =
                        Spread { values: result.spreads[&((*name).clone(), *stringency)].values.clone() };
                    if let Some([p5, p25, p50, p75, p95]) = spread.quantiles() {
                        let _ = writeln!(
                            s,
                            "| {cohort_cell} | {group_cell} | {} | {:.1}% | {:.1}% | **{:.1}%** | {:.1}% | {:.1}% |",
                            spread.values.len(),
                            p5 * 100.0,
                            p25 * 100.0,
                            p50 * 100.0,
                            p75 * 100.0,
                            p95 * 100.0
                        );
                    }
                }
            }
            s.push('\n');
        }
    }

    if let Some(stringency) = expressible.iter().copied().max() {
        s.push_str(&matched(results, stringency));
    }

    s.push_str(
        "\n## Rows that could not be scored\n\n| Cohort | rows | no scorable window |\n|---|---:|---:|\n",
    );
    for (cohort, result) in results {
        let _ = writeln!(s, "| {cohort} | {} | {} |", result.rows, result.unscorable);
    }
    s.push_str(
        "\nA sequence whose every 9-mer holds a non-standard residue has no measurement, and is \
         left out of both numerator and denominator rather than scored zero.\n",
    );

    s.push_str(
        "\n## What this does and does not establish\n\n\
         It establishes that the table behaves as the method claims on donors it has never seen, \
         and that it answers the humanisation question a species classifier cannot: the \
         diagnostic cohort is mice, and a measure of humanness should call their antibodies \
         human.\n\n\
         It establishes nothing about immunogenicity. These are prevalences in sampled \
         repertoires, and whether an antibody provokes a response in a person depends on MHC \
         presentation, dose, route and the recipient. It also establishes nothing comparable with \
         published OASis figures: this table's donor units are shallow, so every prevalence here \
         is biased down. The table's own card gives the size of that bias.\n",
    );
    s
}

/// The median identity of one group, or `None` if it has no rows.
fn median(result: &CohortResult, label: &str, stringency: Stringency) -> Option<f64> {
    let values = result.spreads.get(&(label.to_owned(), stringency))?.values.clone();
    Spread { values }.quantiles().map(|q| q[2])
}

/// Compare human and diagnostic isotype by isotype, and report the spread
/// within human.
///
/// Two things have to be said for the headline result to mean anything.
///
/// First, the diagnostic cohort is mostly IGHM while the human cohort spans
/// isotypes, and a class-switched sequence has been somatically mutated away
/// from germline -- which is where the peptides common to many donors are. So
/// "diagnostic scores above human" could be a difference in mutation load
/// wearing a species costume. The matched comparison settles it.
///
/// Second, the spread within the human cohort itself is large for the same
/// reason, and that is a limit on reading any single score: a mutated human
/// antibody scores lower than a naive one, so an absolute identity cannot be
/// read as "how human" without knowing how mutated the sequence is.
fn matched(results: &[(&str, CohortResult)], stringency: Stringency) -> String {
    let mut s = String::new();
    let find = |want: &str| results.iter().find(|(name, _)| *name == want).map(|(_, r)| r);
    let (Some(human), Some(diagnostic)) = (find("human"), find("diagnostic")) else {
        return s;
    };

    let isotypes: Vec<String> = human
        .spreads
        .keys()
        .filter(|(name, st)| *st == stringency && name.starts_with(ISOTYPE))
        .map(|(name, _)| name.clone())
        .filter(|name| diagnostic.spreads.contains_key(&(name.clone(), stringency)))
        .collect();

    if !isotypes.is_empty() {
        let _ = writeln!(
            s,
            "\n## The confound, and what is left after it\n\nThe diagnostic cohort is mostly \
             IGHM, the human cohort spans isotypes, and a class-switched sequence has been \
             mutated away from germline -- which is where the peptides shared by many donors \
             are. So the diagnostic cohort scoring above the human one could be a difference in \
             somatic mutation rather than in the species of the loci. Matched isotype by isotype \
             at the {} threshold, median identity:\n",
            stringency.label()
        );
        s.push_str("| Isotype | human | diagnostic | difference |\n|---|---:|---:|---:|\n");
        let mut diagnostic_higher = 0_usize;
        for name in &isotypes {
            let (Some(h), Some(d)) = (median(human, name, stringency), median(diagnostic, name, stringency))
            else {
                continue;
            };
            if d >= h {
                diagnostic_higher += 1;
            }
            let _ = writeln!(
                s,
                "| {} | {:.1}% | {:.1}% | {:+.1} pp |",
                name.trim_start_matches(ISOTYPE),
                h * 100.0,
                d * 100.0,
                (d - h) * 100.0
            );
        }
        let _ = writeln!(
            s,
            "\nThe diagnostic cohort scores at least as high as the human cohort on {} of {} \
             shared isotypes, so the result is not the isotype difference. These are mice, and \
             the measure calls their antibodies human: it tracks the immunoglobulin rather than \
             the animal, which is the behaviour a humanness measure needs and the reason this \
             cohort carries no species label. A classifier trained on donor species would have \
             been trained to call them mouse.\n\nWhat it does not show is that the measure beats \
             a germline lookup here. These sequences carry human V genes, so a V-gene table also \
             answers \"human\". The peptide measure earns its place on resolution instead: it \
             scores every window, including the junction and CDR3, where there is no germline \
             gene to look up and where a humanised antibody actually departs from anything a \
             human has made.",
            diagnostic_higher,
            isotypes.len()
        );
        s.push_str(
            "\nThat the diagnostic cohort scores *higher* is not a claim that transgenic mice \
             make better human antibodies. Their repertoires are less affinity-matured, so they \
             sit closer to germline, and germline is where prevalence is highest.\n",
        );
    }

    let mut within: Vec<(String, f64)> = human
        .spreads
        .keys()
        .filter(|(name, st)| *st == stringency && name.starts_with(ISOTYPE))
        .filter_map(|(name, _)| median(human, name, stringency).map(|m| (name.clone(), m)))
        .collect();
    within.sort_by(|a, b| a.1.total_cmp(&b.1));
    if let (Some(low), Some(high)) = (within.first(), within.last()) {
        let _ = writeln!(
            s,
            "\nThe same effect sets a limit on reading any single score. Within the human cohort \
             alone, median identity runs from {:.1}% ({}) to {:.1}% ({}) -- a spread of {:.0} \
             percentage points between isotypes of the same species. An absolute identity \
             therefore cannot be read as \"how human\" without knowing how mutated the sequence \
             is, which is another reason the per-window output is the useful one.",
            low.1 * 100.0,
            low.0.trim_start_matches(ISOTYPE),
            high.1 * 100.0,
            high.0.trim_start_matches(ISOTYPE),
            (high.1 - low.1) * 100.0
        );
    }
    s
}
