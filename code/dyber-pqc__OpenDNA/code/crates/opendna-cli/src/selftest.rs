//! Does this installation actually work on this machine?
//!
//! `opendna models` answers whether the right files are present and intact.
//! That is a different question from whether running them here produces the
//! right numbers. A head can verify perfectly and still give nonsense because
//! the machine resolved a different BLAS, because the GPU path is broken, or
//! because the embedding model in the cache is not the one the head expects.
//!
//! So this runs all four heads over a bundled sequence and compares the
//! results against values recorded from a named build. It deliberately goes
//! through the same functions the `predict` subcommand calls, formatting
//! included, so a pass means the path a user takes works and not that some
//! parallel test harness does.
//!
//! Two kinds of check, because they fail for different reasons:
//!
//! - **Invariants** must hold for any build of these heads. They are facts
//!   about ubiquitin and about the shape of the output, so they catch a
//!   transposed tensor, a truncated sequence, or a model that has stopped
//!   tracking structure at all.
//! - **Pinned values** are numbers measured from one specific build. They
//!   catch small numerical drift that an invariant would pass. They can only
//!   be checked against the build they were taken from, so against any other
//!   build they are reported as not checked rather than quietly skipped or
//!   wrongly failed.
//!
//! A retrain is therefore expected to leave the pinned checks unchecked until
//! someone records new values. That is the honest state: nobody has measured
//! what the new head does on this fixture yet.

use crate::exit::{Classify as _, Fault::Install};
use anyhow::{anyhow, bail, Context, Result};
use opendna_ml::esm2::ModelSize;
use opendna_ml::heads::Heads;
use std::path::Path;

/// Human ubiquitin (P0CG48 residues 1-76).
///
/// Chosen because it is the most thoroughly characterised small protein there
/// is: its fold, its hydrophobic core and its Ile44 patch are textbook, so the
/// invariants below are claims anyone can check against the literature rather
/// than against our own output.
const SEQUENCE: &str = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

/// The FASTA name, which the mutation fixture has to match.
const NAME: &str = "ubiquitin";

/// Buried hydrophobic core positions, 1-based.
///
/// Ile44 is left out on purpose. It is hydrophobic but sits on the surface —
/// it is the patch the proteasome and most ubiquitin-binding domains read — so
/// including it would weaken a claim that is otherwise unambiguous.
const CORE: [usize; 11] = [3, 5, 13, 15, 23, 43, 50, 56, 61, 67, 69];

/// Surface lysines, 1-based: the seven chain-forming sites.
const SURFACE: [usize; 7] = [6, 11, 27, 29, 33, 48, 63];

/// Mutations the stability check scores.
const MUTATIONS: [&str; 2] = ["L8P", "L8M"];

/// A number measured from one build of one head.
struct Pinned {
    /// Which head produced it.
    task: &'static str,
    /// The build id it was measured from. Checking against a different build
    /// would be comparing two models, not testing one.
    build: &'static str,
    /// What was measured, matching the label the run reports.
    what: &'static str,
    /// The recorded value.
    expected: f64,
    /// How far from it this machine may land.
    ///
    /// Wide enough to absorb a different BLAS, a different device, or a
    /// different floating-point summation order; narrow enough that the wrong
    /// model fails. The width is a judgement rather than a measurement — these
    /// values were recorded on one machine — so a failure that lands just
    /// outside is worth reporting rather than widening away.
    tolerance: f64,
}

/// Values measured on 2026-11-15, CPU, from the builds named here.
///
/// The positions are chosen rather than sampled. Disorder and burial are the
/// same architecture differing only in the squash on the output, so at some
/// positions one head's number is close enough to the other's that a pin there
/// would survive the two being swapped. These are positions where the two
/// disagree most, which is the case a pinned value has to catch and a file
/// checksum cannot.
///
/// They had to be chosen again when the disorder head was calibrated.
/// Calibration compresses the probabilities toward zero, which shrinks every
/// such gap: position 15's margin fell from over 0.25 to 0.027, and position 5
/// took its place at 0.129. A fixture's discriminating power is a property of
/// the model it was measured against, not of the fixture.
const PINNED: [Pinned; 9] = [
    Pinned {
        task: "disorder",
        build: "f7d99273e4b9b961",
        what: "disorder probability at 3",
        expected: 0.1594,
        tolerance: 0.01,
    },
    Pinned {
        task: "disorder",
        build: "f7d99273e4b9b961",
        what: "disorder probability at 5",
        expected: 0.1383,
        tolerance: 0.01,
    },
    Pinned {
        task: "disorder",
        build: "f7d99273e4b9b961",
        what: "disorder probability at 47",
        expected: 0.0705,
        tolerance: 0.01,
    },
    Pinned {
        task: "burial",
        build: "1853a82710875d5e",
        what: "accessibility at 5",
        expected: 0.0177,
        tolerance: 0.01,
    },
    Pinned {
        task: "burial",
        build: "1853a82710875d5e",
        what: "accessibility at 48",
        expected: 0.4376,
        tolerance: 0.01,
    },
    Pinned {
        task: "burial",
        build: "1853a82710875d5e",
        what: "accessibility at 47",
        expected: 0.6791,
        tolerance: 0.01,
    },
    Pinned {
        task: "stability",
        build: "0f87e0457b663d81",
        what: "ddg of L8P",
        expected: -0.779,
        tolerance: 0.05,
    },
    Pinned {
        task: "stability",
        build: "0f87e0457b663d81",
        what: "ddg of L8M",
        expected: 0.333,
        tolerance: 0.05,
    },
    // Above 100 because the head is an unbounded regressor on a percentage
    // scale, not a classifier. Recording what it actually says is the point.
    Pinned {
        task: "solubility",
        build: "5ecb3f9d62ea7fae",
        what: "solubility percent",
        expected: 100.1,
        tolerance: 1.0,
    },
];

/// What one check concluded.
enum Outcome {
    /// It held.
    Pass,
    /// It could not be made, with the reason.
    NotChecked(String),
    /// It did not hold, with what was seen.
    Fail(String),
}

/// One check and its result.
struct Check {
    what: String,
    outcome: Outcome,
}

impl Check {
    #[cfg(test)]
    fn pass(what: impl Into<String>) -> Self {
        Self { what: what.into(), outcome: Outcome::Pass }
    }

    /// Assert a claim, recording what was seen when it fails.
    fn holds(what: impl Into<String>, held: bool, saw: impl Into<String>) -> Self {
        let outcome = if held { Outcome::Pass } else { Outcome::Fail(saw.into()) };
        Self { what: what.into(), outcome }
    }

    fn not_checked(what: impl Into<String>, why: impl Into<String>) -> Self {
        Self { what: what.into(), outcome: Outcome::NotChecked(why.into()) }
    }
}

/// What the four heads said about the fixture.
struct Numbers {
    /// Per-residue disorder probability.
    probabilities: Vec<f64>,
    /// Per-residue relative solvent accessibility.
    accessibility: Vec<f64>,
    /// Whole-protein solubility percentage.
    percent: f64,
    /// One score per entry of `MUTATIONS`, in that order.
    ddg: Vec<f64>,
    /// Disorder from a second, identical run.
    repeated: Vec<f64>,
    /// The threshold the head says to call a residue disordered at.
    threshold: f64,
}

/// Every measurement the run produced, labelled the way `PINNED` names them.
struct Measured {
    values: Vec<(&'static str, String, f64)>,
}

impl Measured {
    fn get(&self, task: &str, what: &str) -> Option<f64> {
        self.values.iter().find(|(t, w, _)| *t == task && w == what).map(|(_, _, v)| *v)
    }
}

/// Pull a numeric column out of a rendered table, keyed by position.
fn column(table: &crate::Table, column: &str) -> Result<Vec<f64>> {
    let index = table
        .columns
        .iter()
        .position(|c| *c == column)
        .ok_or_else(|| anyhow!("the result table has no {column} column"))?;
    table
        .rows
        .iter()
        .map(|row| {
            row.get(index)
                .ok_or_else(|| anyhow!("a row is missing its {column} field"))
                .and_then(|v| v.parse::<f64>().with_context(|| format!("parsing {v:?} as a number")))
        })
        .collect()
}

/// The value at a 1-based position, or NaN if the run is short.
///
/// NaN rather than a default, because every comparison against it is false:
/// a missing value reports as a failure instead of passing by accident.
fn at(values: &[f64], position: usize) -> f64 {
    values.get(position - 1).copied().unwrap_or(f64::NAN)
}

/// Mean of the values at the given 1-based positions.
fn mean_at(values: &[f64], positions: &[usize]) -> f64 {
    let total: f64 = positions.iter().filter_map(|p| values.get(p - 1)).sum();
    #[allow(clippy::cast_precision_loss)]
    let count = positions.len() as f64;
    total / count
}

/// Claims that must hold for any build of these heads.
///
/// Each one is either a fact about the output's shape or a fact about
/// ubiquitin. A head that has drifted a little passes all of them; a head that
/// is the wrong model, reading the wrong tensor, or truncating the sequence
/// fails one.
fn invariants(n: &Numbers) -> Vec<Check> {
    let residues = SEQUENCE.len();

    // Shape. A head that silently drops or duplicates residues produces a
    // table that reads correctly and is wrong at every position after the
    // first, which is the hardest kind of error to notice downstream.
    let mut checks = vec![
        Check::holds(
            "disorder returns one value per residue",
            n.probabilities.len() == residues,
            format!("{} values for {residues} residues", n.probabilities.len()),
        ),
        Check::holds(
            "burial returns one value per residue",
            n.accessibility.len() == residues,
            format!("{} values for {residues} residues", n.accessibility.len()),
        ),
        Check::holds(
            "stability scores every mutation asked for",
            n.ddg.len() == MUTATIONS.len(),
            format!("{} scores for {} mutations", n.ddg.len(), MUTATIONS.len()),
        ),
    ];

    // Scale. Both per-residue heads squash their output, so a value outside
    // the range means the wrong tensor is being read, not a bad prediction.
    let off_scale = n.probabilities.iter().filter(|&&p| !(0.0..=1.0).contains(&p)).count();
    checks.push(Check::holds(
        "disorder probabilities lie in [0, 1]",
        off_scale == 0,
        format!("{off_scale} values outside the range"),
    ));
    let off_rsa = n.accessibility.iter().filter(|&&a| !(0.0..=1.25).contains(&a)).count();
    checks.push(Check::holds(
        "accessibility lies in [0, 1.25]",
        off_rsa == 0,
        format!("{off_rsa} values outside the range"),
    ));

    // Structure. Ubiquitin's core is buried and its chain-forming lysines are
    // not. A head that has stopped reading structure at all passes every check
    // above and fails this one.
    let core = mean_at(&n.accessibility, &CORE);
    let surface = mean_at(&n.accessibility, &SURFACE);
    checks.push(Check::holds(
        "the hydrophobic core reads as less accessible than the surface lysines",
        core + 0.1 < surface,
        format!("core {core:.3}, surface lysines {surface:.3}"),
    ));

    // Biology. Ubiquitin is the textbook folded protein: a compact beta-grasp
    // domain with no intrinsic disorder. A head that calls a substantial part
    // of it disordered has stopped predicting disorder, and no tolerance on
    // three individual positions would notice.
    let called = n.probabilities.iter().filter(|p| **p >= n.threshold).count();
    checks.push(Check::holds(
        "almost none of ubiquitin is called disordered",
        called <= residues / 10,
        format!("{called} of {residues} residues called at threshold {:.2}", n.threshold),
    ));

    // Chemistry. A proline in the middle of ubiquitin's first beta strand
    // breaks the backbone hydrogen bonding; L8M is conservative. Any head
    // worth running puts them in this order.
    let (l8p, l8m) = (at(&n.ddg, 1), at(&n.ddg, 2));
    checks.push(Check::holds(
        "L8P is predicted more destabilising than L8M",
        l8p < l8m,
        format!("L8P {l8p:+.3}, L8M {l8m:+.3}"),
    ));

    // Determinism. When this fails every number the installation produces is
    // unreproducible, and nothing else it reports matters.
    checks.push(Check::holds(
        "a second identical run gives identical values",
        n.repeated == n.probabilities,
        "the two runs differ".to_string(),
    ));
    checks
}

/// The recorded numbers, each against the build it was measured from.
fn pinned(n: &Numbers, build_of: &impl Fn(&str) -> Option<String>) -> Vec<Check> {
    let measured = Measured {
        values: vec![
            ("disorder", "disorder probability at 3".into(), at(&n.probabilities, 3)),
            ("disorder", "disorder probability at 5".into(), at(&n.probabilities, 5)),
            ("disorder", "disorder probability at 47".into(), at(&n.probabilities, 47)),
            ("burial", "accessibility at 5".into(), at(&n.accessibility, 5)),
            ("burial", "accessibility at 47".into(), at(&n.accessibility, 47)),
            ("burial", "accessibility at 48".into(), at(&n.accessibility, 48)),
            ("stability", "ddg of L8P".into(), at(&n.ddg, 1)),
            ("stability", "ddg of L8M".into(), at(&n.ddg, 2)),
            ("solubility", "solubility percent".into(), n.percent),
        ],
    };
    let mut checks = Vec::with_capacity(PINNED.len());
    for pin in &PINNED {
        let installed = build_of(pin.task);
        if installed.as_deref() != Some(pin.build) {
            checks.push(Check::not_checked(
                pin.what,
                format!(
                    "recorded for {} build {}, this install runs {}",
                    pin.task,
                    pin.build,
                    installed.as_deref().unwrap_or("an unidentified build")
                ),
            ));
            continue;
        }
        let Some(saw) = measured.get(pin.task, pin.what) else {
            checks.push(Check::not_checked(pin.what, "the run produced no such value"));
            continue;
        };
        checks.push(Check::holds(
            format!("{} is {:.4} ± {}", pin.what, pin.expected, pin.tolerance),
            (saw - pin.expected).abs() <= pin.tolerance,
            format!("{saw:.4}, off by {:.4}", saw - pin.expected),
        ));
    }
    checks
}

/// Run every check and report.
///
/// # Errors
/// If a head is missing or corrupt, if the embedding model cannot be loaded,
/// or if any check fails. Failure is an error rather than a warning because
/// the whole purpose is to be runnable in an installation script.
///
/// With `offline` set it will not reach the network, so on a staged bundle it
/// tests the same files the predictions will use.
///
/// `device` is passed in rather than chosen here, so the same checks can be
/// run on the CPU and on a GPU and the two compared. A head that reproduces
/// its recorded numbers on one device and not the other is worth knowing
/// about, and nothing else would report it.
///
/// `budget` comes in for the same reason, though the fixture is a single
/// sequence and so has nothing to be batched with: passing it means the
/// self-test exercises whatever a real run on this machine would use rather
/// than a path of its own.
pub fn run(
    dir: &Path,
    embedder: ModelSize,
    offline: bool,
    device: candle_core::Device,
    budget: usize,
) -> Result<(), crate::exit::Failure> {
    let found = crate::inventory::scan(dir, embedder);
    let unusable: Vec<String> =
        found.iter().filter(|h| !h.ready()).map(|h| format!("  {:11} {}", h.task, h.describe())).collect();
    if !unusable.is_empty() {
        return Err(crate::exit::Failure {
            fault: crate::exit::Fault::Install,
            error: anyhow!(
                "selftest needs all four heads installed and intact. In {}:\n{}",
                dir.display(),
                unusable.join("\n")
            ),
        });
    }
    let metadata_of = |task: &str| -> Option<opendna_ml::heads::HeadMetadata> {
        found.iter().find(|h| h.task == task).and_then(|h| h.metadata.clone())
    };
    let threshold = metadata_of("disorder")
        .and_then(|m| m.calibration.map(|c| c.threshold))
        .unwrap_or(crate::DISORDER_THRESHOLD);

    let build_of = |task: &str| -> Option<String> {
        found.iter().find(|h| h.task == task).and_then(|h| h.metadata.as_ref()).and_then(|m| m.build.clone())
    };

    let description = opendna_ml::device_description(&device);
    println!("Self-test on {description}");
    println!("Fixture: {NAME}, {} residues", SEQUENCE.len());

    let heads = Heads::load(dir, embedder, &device).blame(Install)?;
    let model = crate::staging::open_embedder(dir, embedder, &device, offline).blame(Install)?;
    println!(
        "Embedding model: {} ({})",
        embedder.repo(),
        if crate::staging::staged_embedder(dir).is_some() {
            "staged in the models directory"
        } else {
            "Hugging Face cache"
        }
    );
    let proteins = vec![(NAME.to_string(), SEQUENCE.to_string())];
    let mutations: Vec<(String, String)> =
        MUTATIONS.iter().map(|m| (NAME.to_string(), (*m).to_string())).collect();
    // Through the same resolution the CLI does, so the fixture exercises the
    // wild-type agreement check as well as the head.
    let resolved = crate::resolve_mutations(&proteins, &mutations).blame(Install)?;

    // The same calls `predict` makes, in the same order.
    let disorder_runs =
        crate::disorder(&proteins, &heads, &model, &device, budget, threshold).blame(Install)?;
    let disorder_table = crate::disorder_table(&disorder_runs, threshold, None);
    let burial_table = crate::burial(&proteins, &heads, &model, &device, budget).blame(Install)?;
    let solubility_table = crate::solubility(&proteins, &heads, &model, &device, budget).blame(Install)?;
    let stability_table = crate::stability(&proteins, &resolved, &heads, &model, &device).blame(Install)?;

    let numbers = Numbers {
        repeated: {
            // Nothing here is sampled, so a second run on this machine has to
            // agree with the first. Run it before anything is compared, so the
            // determinism check cannot be skipped by an earlier failure.
            let again = crate::disorder_table(
                &crate::disorder(&proteins, &heads, &model, &device, budget, threshold)?,
                threshold,
                None,
            );
            column(&again, "probability")?
        },
        probabilities: column(&disorder_table, "probability")?,
        accessibility: column(&burial_table, "accessibility")?,
        percent: at(&column(&solubility_table, "solubility_percent")?, 1),
        ddg: column(&stability_table, "ddg_kcal_per_mol")?,
        threshold,
    };

    let mut checks = invariants(&numbers);
    checks.extend(pinned(&numbers, &build_of));
    report(&checks).blame(crate::exit::Fault::Checks)
}

/// Print the checks and decide the exit status.
fn report(checks: &[Check]) -> Result<()> {
    println!();
    let mut failed = 0_usize;
    let mut unchecked = 0_usize;
    for check in checks {
        match &check.outcome {
            Outcome::Pass => println!("  {:<11} {}", "ok", check.what),
            Outcome::NotChecked(why) => {
                unchecked += 1;
                println!("  {:<11} {} — {why}", "not checked", check.what);
            }
            Outcome::Fail(saw) => {
                failed += 1;
                println!("  {:<11} {} — saw {saw}", "FAILED", check.what);
            }
        }
    }
    let passed = checks.len() - failed - unchecked;
    println!("\n{passed} passed, {failed} failed, {unchecked} not checked");
    if failed > 0 {
        bail!(
            "this installation does not reproduce its recorded results. Do not trust numbers from \
             it until that is explained: re-download the heads, and if the failures persist report \
             them with the machine, the device line above and this output."
        );
    }
    if unchecked > 0 {
        println!(
            "\nThe numeric checks above were recorded from particular builds and this install \
             runs others, so only the invariants were tested. Record values for these builds to \
             get the rest back."
        );
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// The fixture's positions have to name the residues the comments claim,
    /// or every structural check is about the wrong side chains.
    #[test]
    fn the_fixture_positions_name_the_residues_they_claim() {
        let residues: Vec<char> = SEQUENCE.chars().collect();
        assert_eq!(residues.len(), 76, "ubiquitin is 76 residues");
        for p in CORE {
            assert!(
                matches!(residues[p - 1], 'I' | 'L' | 'V' | 'F' | 'M'),
                "core position {p} is {}, not hydrophobic",
                residues[p - 1]
            );
        }
        for p in SURFACE {
            assert_eq!(residues[p - 1], 'K', "surface position {p} is not a lysine");
        }
        // The mutation labels have to agree with the sequence too, or the run
        // stops with a parse error rather than a result.
        assert_eq!(residues[7], 'L', "both fixture mutations are at L8");
    }

    #[test]
    fn a_pinned_value_exists_for_every_head() {
        for task in crate::inventory::TASKS {
            assert!(PINNED.iter().any(|p| p.task == task), "no pinned value for {task}");
        }
    }

    #[test]
    fn a_mean_over_positions_is_one_based() {
        let values = [0.0, 1.0, 2.0, 3.0];
        assert!((mean_at(&values, &[1, 2]) - 0.5).abs() < 1e-12);
        assert!((mean_at(&values, &[4]) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn a_failing_check_makes_the_report_an_error() {
        let ok = report(&[Check::pass("something")]);
        assert!(ok.is_ok());

        let bad = report(&[Check::holds("something", false, "a wrong number")]);
        let err = bad.unwrap_err().to_string();
        assert!(err.contains("does not reproduce"), "{err}");

        // Not checked is not a failure: a retrained head has not been measured
        // on this fixture, which is different from being wrong on it.
        let partial = report(&[Check::not_checked("something", "another build")]);
        assert!(partial.is_ok());
    }

    #[test]
    fn a_column_is_read_by_name_and_not_by_position() {
        let table = crate::Table {
            columns: vec!["name", "probability"],
            rows: vec![vec!["p".into(), "0.25".into()], vec!["p".into(), "0.75".into()]],
        };
        assert_eq!(column(&table, "probability").unwrap(), vec![0.25, 0.75]);
        assert!(column(&table, "absent").is_err());
    }
}
