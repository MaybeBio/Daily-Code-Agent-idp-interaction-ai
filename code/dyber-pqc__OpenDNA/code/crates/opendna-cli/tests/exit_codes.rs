//! The exit-code contract, checked against the real binary.
//!
//! A pipeline reads the status, not the message. The codes are documented in
//! the README and in `exit.rs`, and a documented contract nobody exercises is
//! a comment. So this runs the compiled `opendna` and asserts what each kind
//! of failure returns.
//!
//! It needs trained heads, which are not in the repository: they are hundreds
//! of kilobytes of weights published separately. With none installed the test
//! prints what it would have needed and passes, rather than failing on every
//! machine that has not downloaded them. That is a deliberate trade — a test
//! that fails for everybody gets disabled, and then it checks nothing at all.
//!
//! Every fixture is built in a temporary directory from a copy of the real
//! models, so nothing here can damage the installation it is testing.
//!
//! The cases split by cost. Most of the contract is about refusing before any
//! inference happens, which is fast even unoptimised. Two cases have to embed
//! a sequence, and candle in a debug build turns that from seconds into
//! minutes, so they are `#[ignore]`d and run deliberately:
//!
//! ```text
//! cargo test --release -p opendna-cli -- --ignored
//! ```

#![allow(clippy::unwrap_used, clippy::print_stdout)]

use std::path::{Path, PathBuf};
use std::process::Command;

/// Where the heads are, if they are anywhere.
fn models() -> Option<PathBuf> {
    let dir = std::env::var_os("OPENDNA_MODELS").map_or_else(|| PathBuf::from("models"), PathBuf::from);
    // The second candidate is relative to the workspace root, since cargo may
    // run us from either there or the crate directory.
    [dir.clone(), Path::new("../..").join(&dir)]
        .into_iter()
        .find(|candidate| candidate.join("disorder.safetensors").is_file())
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("opendna-exit-{name}-{}", std::process::id()));
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Copy a directory's files one level deep, which is all a models directory is.
fn copy_models(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        if entry.file_type().unwrap().is_file() {
            std::fs::copy(entry.path(), to.join(entry.file_name())).unwrap();
        }
    }
}

fn sha256(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

/// Held for the whole of any test that spawns a run which embeds something.
///
/// The Hugging Face cache cannot be read by two processes at once: on Windows
/// the second one is refused the weights file. Cargo runs tests in parallel,
/// so without this the two inference tests race and whichever loses reports a
/// broken installation. Not a workaround for a test bug -- it is the same
/// constraint real users have, and `staging::open_embedder` says so in its
/// error.
static EMBEDDER: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Take that lock, ignoring poisoning: a panic in the other test tells us
/// nothing about whether the cache is readable now.
fn one_at_a_time() -> std::sync::MutexGuard<'static, ()> {
    EMBEDDER.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// Run the binary and return (status code, stdout).
fn run(args: &[&str]) -> (i32, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_opendna")).args(args).output().unwrap();
    (output.status.code().unwrap_or(-1), String::from_utf8_lossy(&output.stdout).into_owned())
}

/// A stand-in for a folding engine: answers `--version`, and writes a
/// two-residue structure into whatever directory it is given with `-o`.
///
/// Not a mock of anything in this repository. It stands where `ESMFold` stands and is
/// driven through the same discovery, the same command line and the same
/// checks, which is the only way to exercise a contract whose other side is
/// a 2.6 GB download nobody can require to run a test.
fn stand_in_engine(dir: &Path) -> PathBuf {
    // One Ca per residue is all the matcher reads: MET then GLN, with the
    // two confidence values in the B-factor column.
    let met = "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 88.00           C";
    let gln = "ATOM      2  CA  GLN A   2       3.800   0.000   0.000  1.00 78.00           C";

    #[cfg(windows)]
    {
        let path = dir.join("stand-in-fold.cmd");
        std::fs::write(
            &path,
            format!(
                "@echo off\r\n\
                 if \"%1\"==\"--version\" (echo stand-in 1.0 & exit /b 0)\r\n\
                 > \"%4\\model.pdb\" echo {met}\r\n\
                 >>\"%4\\model.pdb\" echo {gln}\r\n\
                 >>\"%4\\model.pdb\" echo END\r\n\
                 exit /b 0\r\n"
            ),
        )
        .unwrap();
        path
    }

    #[cfg(not(windows))]
    {
        use std::os::unix::fs::PermissionsExt as _;
        let path = dir.join("stand-in-fold.sh");
        std::fs::write(
            &path,
            format!(
                "#!/bin/sh\n\
                 if [ \"$1\" = \"--version\" ]; then echo 'stand-in 1.0'; exit 0; fi\n\
                 printf '%s\\n%s\\nEND\\n' '{met}' '{gln}' > \"$4/model.pdb\"\n\
                 exit 0\n"
            ),
        )
        .unwrap();
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
        path
    }
}

/// Run the binary with the environment a fold test needs.
///
/// Both override variables are always set: pointing the unused one at a file
/// that is not there is what stops a machine that happens to have `ColabFold`
/// installed from turning "no engine" into "two engines".
fn run_fold(engine: Option<&Path>, args: &[&str]) -> (i32, String, String) {
    let nowhere = std::env::temp_dir().join("opendna-no-such-engine");
    let esmfold = engine.map_or(nowhere.clone(), Path::to_path_buf);
    let output = Command::new(env!("CARGO_BIN_EXE_opendna"))
        .args(args)
        .env("OPENDNA_ESMFOLD", &esmfold)
        .env("OPENDNA_COLABFOLD", &nowhere)
        .output()
        .unwrap();
    (
        output.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// Run the binary with no engines and no models, for commands that need none.
fn run_bare(args: &[&str]) -> (i32, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_opendna")).args(args).output().unwrap();
    let said =
        format!("{}{}", String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
    (output.status.code().unwrap_or(-1), said)
}

/// `humanness`'s share of the contract, which needs no heads and no engines.
///
/// The case worth having is the last one: a table too narrow to express the
/// threshold asked for. Nothing crashes and every lookup succeeds; the
/// resulting percentage is simply a number that measured nothing, and it is
/// refused for that reason alone.
#[test]
fn humanness_refuses_in_the_ways_the_contract_promises() {
    let base = scratch("humanness");
    let fasta = base.join("ab.fasta");
    std::fs::write(&fasta, ">h1\nACDEFGHIKLMNPQRSTVWY\n").unwrap();
    let missing = base.join("absent.bin");
    let garbage = base.join("garbage.bin");
    std::fs::write(&garbage, b"this is not a prevalence table, and is not 16 bytes of magic").unwrap();

    // A ten-donor table, written through the real writer so the format under
    // test is the one the builder produces. Ten puts the boundary exactly on
    // a threshold: 1/10 is 10%, so `medium` is expressible and `loose` (1%)
    // is not.
    let table = base.join("ten.bin");
    {
        use opendna_analysis::humanness::PrevalenceTable;
        let built = PrevalenceTable::build(std::iter::repeat_n(vec!["ACDEFGHIKLMNPQRSTVWY"], 10), 1);
        assert_eq!(built.donors(), 10);
        built.write(std::fs::File::create(&table).unwrap()).unwrap();
    }

    let f = fasta.to_str().unwrap();
    let t = table.to_str().unwrap();
    let mut failures = Vec::new();

    // No table: an installation problem, and there is nothing to fall back on
    // -- a score without a table would be a number with no denominator.
    let (code, said) = run_bare(&["humanness", "-i", f, "-t", missing.to_str().unwrap()]);
    check(&mut failures, 3, "no prevalence table", code, "opendna-data build humanness", &said);

    // A file that is not a table at all.
    let (code, said) = run_bare(&["humanness", "-i", f, "-t", garbage.to_str().unwrap()]);
    check(&mut failures, 3, "table is not a table", code, "expected magic", &said);

    // Unreadable input.
    let (code, said) = run_bare(&["humanness", "-i", base.join("none.fasta").to_str().unwrap(), "-t", t]);
    check(&mut failures, 4, "missing input", code, "reading", &said);

    // A sequence with no scorable window: refusing beats scoring it zero,
    // which would assert it looks nothing like a human antibody.
    let allx = base.join("allx.fasta");
    std::fs::write(&allx, format!(">x\n{}\n", "X".repeat(30))).unwrap();
    let (code, said) = run_bare(&["humanness", "-i", allx.to_str().unwrap(), "-t", t]);
    check(&mut failures, 4, "nothing scorable", code, "no scorable window", &said);

    // A threshold finer than ten donors can express. Checks ran and did not
    // hold: the install is fine and so is the input, and every lookup would
    // have succeeded -- the percentage would simply have measured nothing.
    let (code, said) = run_bare(&["humanness", "-i", f, "-t", t, "--stringency", "loose"]);
    check(&mut failures, 5, "threshold not expressible", code, "cannot be evaluated", &said);

    // And the cases that must work: the default threshold, which ten donors
    // can express exactly, and the strictest one.
    let (code, said) = run_bare(&["humanness", "-i", f, "-t", t]);
    check(&mut failures, 0, "the default threshold", code, "of donors", &said);
    let (code, said) = run_bare(&["humanness", "-i", f, "-t", t, "--stringency", "strict"]);
    check(&mut failures, 0, "the strictest threshold", code, "100.0%", &said);

    assert!(failures.is_empty(), "{}", failures.join("\n\n"));
}

/// Record a case that did not behave as the contract promises.
fn check(into: &mut Vec<String>, want: i32, what: &str, got: i32, needle: &str, said: &str) {
    if got != want {
        into.push(format!("{what}: wanted {want}, got {got}\n{said}"));
    } else if !said.contains(needle) {
        into.push(format!("{what}: exit {got} was right but it did not say {needle:?}\n{said}"));
    }
}

/// `fold`'s share of the contract, which needs no heads at all.
///
/// The case that matters most is the last one: an engine that ran, exited
/// cleanly and returned a structure of the wrong protein. That is not a
/// crash, and nothing in the engine's own output says anything is wrong; it
/// is caught only because the chain that came back is compared with the
/// sequence that went out.
#[test]
fn folding_refuses_in_the_ways_the_contract_promises() {
    let base = scratch("fold");
    let right = base.join("right.fasta");
    std::fs::write(&right, ">p1\nMQ\n").unwrap();
    let wrong = base.join("wrong.fasta");
    std::fs::write(&wrong, ">p1\nMQAAK\n").unwrap();
    let clash = base.join("clash.fasta");
    std::fs::write(&clash, ">p 1\nMQ\n>p|1\nMA\n").unwrap();
    let engine = stand_in_engine(&base);
    let out = |name: &str| base.join(name).to_str().unwrap().to_string();

    let mut failures = Vec::new();

    // No engine: an installation problem, and nothing to fall back on.
    let (code, _, said) = run_fold(None, &["fold", "-i", right.to_str().unwrap(), "-o", &out("a")]);
    check(&mut failures, 3, "no engine installed", code, "no folding engine is installed", &said);

    // Two sequences whose names become one filename. Refused before the
    // engine is started, which is the point: the engine is the slow part.
    let (code, _, said) = run_fold(Some(&engine), &["fold", "-i", clash.to_str().unwrap(), "-o", &out("b")]);
    check(&mut failures, 4, "names that collide", code, "would both be written to", &said);
    if said.contains("Running:") {
        failures.push("a name collision started the engine before refusing".to_string());
    }

    // The engine works and the structure is the sequence that was sent.
    let (code, printed, said) =
        run_fold(Some(&engine), &["fold", "-i", right.to_str().unwrap(), "-o", &out("c")]);
    check(&mut failures, 0, "a working fold", code, "structure(s)", &printed);
    if code == 0 {
        let written = base.join("c").join("p1.pdb");
        if !written.is_file() {
            failures.push(format!("nothing was written to {}", written.display()));
        }
        if !base.join("c").join("p1.manifest.json").is_file() {
            failures.push("no manifest was written beside the structure".to_string());
        }
        // Matched by its residues: the engine named its file `model.pdb`,
        // which has nothing to do with `p1`.
        if !said.contains("Running:") {
            failures.push("the command that ran should be printed".to_string());
        }
    }

    // The engine returned a structure of something else, and said nothing
    // was wrong. A check ran and did not hold.
    let (code, _, said) = run_fold(Some(&engine), &["fold", "-i", wrong.to_str().unwrap(), "-o", &out("d")]);
    check(
        &mut failures,
        5,
        "a structure of the wrong sequence",
        code,
        "did not return a structure matching",
        &said,
    );

    std::fs::remove_dir_all(&base).ok();
    assert!(failures.is_empty(), "the fold contract does not hold:\n  {}", failures.join("\n  "));
}

/// The failures that are refused before any inference starts.
// Long because it is a list. Splitting a table of cases across functions to
// satisfy a line count would make it harder to read and easier to forget a
// case in.
#[allow(clippy::too_many_lines)]
#[test]
fn every_kind_of_failure_returns_the_status_the_contract_promises() {
    let Some(installed) = models() else {
        println!("no heads installed; skipping. Set OPENDNA_MODELS or run `opendna stage`.");
        return;
    };
    let base = scratch("fixtures");
    let fasta = base.join("ubq.fasta");
    std::fs::write(
        &fasta,
        ">ubiquitin\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n",
    )
    .unwrap();

    // Each case is (expected status, what it is, arguments).
    let empty = base.join("empty");
    std::fs::create_dir_all(&empty).unwrap();

    // Corrupt: one byte appended, so the weights no longer hash to what the
    // metadata records. This is what a truncated download looks like.
    let corrupt = base.join("corrupt");
    copy_models(&installed, &corrupt);
    {
        let path = corrupt.join("disorder.safetensors");
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.push(b'x');
        std::fs::write(&path, bytes).unwrap();
    }

    // Swapped: burial's weights presented as the disorder head, with the
    // metadata rewritten to match so the file verifies. Nothing about the
    // files is wrong; the numbers are.
    let swapped = base.join("swapped");
    copy_models(&installed, &swapped);
    {
        let weights = std::fs::read(swapped.join("burial.safetensors")).unwrap();
        std::fs::write(swapped.join("disorder.safetensors"), &weights).unwrap();
        let path = swapped.join("disorder.metadata.json");
        let mut meta: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        meta["weights_sha256"] = sha256(&weights).into();
        std::fs::write(&path, serde_json::to_string_pretty(&meta).unwrap()).unwrap();
    }

    let bad_fasta = base.join("bad.fasta");
    std::fs::write(&bad_fasta, "not a fasta at all\n").unwrap();
    let wrong_mutation = base.join("wrong.csv");
    // Position 8 is a leucine, so a label calling it Q cannot be scored: a
    // list numbered against a different construct must stop the run.
    std::fs::write(&wrong_mutation, "name,mutation\nubiquitin,Q8P\n").unwrap();

    // A lockfile pinning an engine that is certainly not installed, so the
    // check has something to fail against. Engines are a different kind of
    // dependency from heads and reach the same exit code, which is the part
    // worth holding: a pipeline reads the number, not the noun.
    let engines_lock = base.join("engines.lock");
    std::fs::write(
        &engines_lock,
        "engine vina version = AutoDock Vina v1.2.5 sha256 = ".to_owned()
            + &"ab".repeat(32)
            + " identification = executable
",
    )
    .unwrap();
    let engines_lock = engines_lock.to_str().unwrap().to_owned();

    let checkpoint = base.join("run.jsonl");
    let checkpoint = checkpoint.to_str().unwrap().to_owned();
    let installed = installed.to_str().unwrap().to_owned();
    let (fasta, bad_fasta) = (fasta.to_str().unwrap(), bad_fasta.to_str().unwrap());
    let wrong_mutation = wrong_mutation.to_str().unwrap();
    let (empty, corrupt, swapped) =
        (empty.to_str().unwrap(), corrupt.to_str().unwrap(), swapped.to_str().unwrap());

    let cases: Vec<(i32, &str, Vec<&str>)> = vec![
        (2, "an unknown flag", vec!["predict", "--task", "disorder", "--nonsense"]),
        (
            2,
            "stability without --mutations",
            vec!["--models", &installed, "predict", "--task", "stability", "-i", fasta],
        ),
        (
            2,
            "the CAID format on the wrong task",
            vec!["--models", &installed, "predict", "--task", "burial", "--format", "caid", "-i", fasta],
        ),
        (3, "no heads installed", vec!["--models", empty, "predict", "--task", "disorder", "-i", fasta]),
        (3, "corrupt weights", vec!["--models", corrupt, "predict", "--task", "disorder", "-i", fasta]),
        (
            3,
            "--offline with nothing staged",
            vec!["--models", &installed, "--offline", "predict", "--task", "disorder", "-i", fasta],
        ),
        (
            4,
            "a file that is not FASTA",
            vec!["--models", &installed, "predict", "--task", "disorder", "-i", bad_fasta],
        ),
        (3, "an engine lockfile that no longer matches", vec!["engines", "--check", &engines_lock]),
        (
            2,
            "a chunk size of zero",
            vec![
                "--models",
                &installed,
                "predict",
                "--task",
                "disorder",
                "-i",
                fasta,
                "--checkpoint",
                &checkpoint,
                "--chunk",
                "0",
            ],
        ),
        (
            2,
            "a checkpoint with the CAID format",
            vec![
                "--models",
                &installed,
                "predict",
                "--task",
                "disorder",
                "--format",
                "caid",
                "-i",
                fasta,
                "--checkpoint",
                &checkpoint,
            ],
        ),
        (
            4,
            "a mutation the sequence contradicts",
            vec![
                "--models",
                &installed,
                "predict",
                "--task",
                "stability",
                "-i",
                fasta,
                "--mutations",
                wrong_mutation,
            ],
        ),
    ];
    let _ = swapped;

    let mut wrong = Vec::new();
    for (want, what, args) in &cases {
        let (got, stdout) = run(args);
        if got != *want {
            wrong.push(format!("{what}: wanted {want}, got {got}"));
        }
        // `predict` streams results to standard output, so a failing run must
        // leave it empty: a pipe collecting results must never receive an
        // error as though it were data. `selftest` is a report command whose
        // standard output is the report, so this does not apply to it.
        if got != 0 && args.contains(&"predict") && !stdout.is_empty() {
            wrong.push(format!("{what}: wrote {} bytes to stdout while failing", stdout.len()));
        }
    }
    // A refused run must not leave a checkpoint behind. A file created by a
    // run that never started would be resumed by the next one.
    if Path::new(&checkpoint).exists() {
        wrong.push("a refused run created its checkpoint file".to_string());
    }
    std::fs::remove_dir_all(&base).ok();
    assert!(wrong.is_empty(), "the exit-code contract does not hold:\n  {}", wrong.join("\n  "));
}

/// Resuming is only safe if continuing the wrong checkpoint is impossible.
///
/// Needs real inference, because the checkpoint under test has to be one the
/// binary itself wrote: hand-building the header here would be a second
/// statement of its shape, and the two would drift.
#[test]
#[ignore = "embeds a sequence; run with --release -- --ignored"]
fn a_checkpoint_from_another_run_is_refused_and_an_interrupted_one_is_continued() {
    let _one_at_a_time = one_at_a_time();
    let Some(installed) = models() else {
        println!("no heads installed; skipping.");
        return;
    };
    let installed = installed.to_str().unwrap().to_owned();
    let base = scratch("resume");
    let ubq = base.join("ubq.fasta");
    std::fs::write(
        &ubq,
        ">ubiquitin\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n",
    )
    .unwrap();
    // Two sequences, so a chunk size of one gives a run with two chunks.
    let pair = base.join("pair.fasta");
    std::fs::write(
        &pair,
        ">ubiquitin\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n\
         >fragment\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEG\n",
    )
    .unwrap();

    let reference = base.join("reference.csv");
    let checkpoint = base.join("run.jsonl");
    let predict = |input: &Path, output: &Path, checkpoint: &Path| -> (i32, String) {
        let output = Command::new(env!("CARGO_BIN_EXE_opendna"))
            .args([
                "--models",
                &installed,
                "predict",
                "--task",
                "disorder",
                "-i",
                input.to_str().unwrap(),
                "-o",
                output.to_str().unwrap(),
                "--checkpoint",
                checkpoint.to_str().unwrap(),
                "--chunk",
                "1",
                "--no-manifest",
            ])
            .output()
            .unwrap();
        (output.status.code().unwrap_or(-1), String::from_utf8_lossy(&output.stderr).into_owned())
    };

    let (code, stderr) = predict(&pair, &reference, &checkpoint);
    assert_eq!(code, 0, "the reference run must finish:\n{stderr}");

    // The same checkpoint against a different input. The rows in it describe
    // sequences this run is not predicting, so continuing would produce a
    // table that is half one run and half another.
    let (code, stderr) = predict(&ubq, &base.join("wrong.csv"), &checkpoint);
    assert_eq!(code, 4, "a checkpoint for another input must be refused:\n{stderr}");
    assert!(stderr.contains("input digest"), "the refusal must name the field:\n{stderr}");

    // An interrupted run, as a kill partway through a write leaves it: the
    // first chunk intact and the second one half-written.
    let interrupted = base.join("interrupted.jsonl");
    let full = std::fs::read_to_string(&checkpoint).unwrap();
    let mut lines = full.lines();
    let header = lines.next().unwrap();
    let first = lines.next().unwrap();
    let torn = &lines.next().unwrap()[..20];
    std::fs::write(&interrupted, format!("{header}\n{first}\n{torn}")).unwrap();

    let resumed = base.join("resumed.csv");
    let (code, stderr) = predict(&pair, &resumed, &interrupted);
    assert_eq!(code, 0, "an interrupted run must continue:\n{stderr}");
    assert!(stderr.contains("Resuming"), "it should say it is resuming:\n{stderr}");

    // The point of the whole exercise: being interrupted changed nothing.
    assert_eq!(
        std::fs::read_to_string(&resumed).unwrap(),
        std::fs::read_to_string(&reference).unwrap(),
        "a resumed run must produce the same bytes as an uninterrupted one"
    );
    std::fs::remove_dir_all(&base).ok();
}

/// The two cases that need real inference.
///
/// Slow unoptimised, so run with `--release -- --ignored`. They are still
/// worth having: a contract that only describes failures says nothing about
/// the success it is contrasted with, and exit 5 is the only code that proves
/// the checks themselves can fail the run.
#[test]
#[ignore = "embeds a sequence; run with --release -- --ignored"]
fn a_working_run_exits_zero_and_a_wrong_model_exits_five() {
    let _one_at_a_time = one_at_a_time();
    let Some(installed) = models() else {
        println!("no heads installed; skipping.");
        return;
    };
    let base = scratch("inference");
    let fasta = base.join("ubq.fasta");
    std::fs::write(
        &fasta,
        ">ubiquitin\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n",
    )
    .unwrap();

    let (code, out) = run(&[
        "--models",
        installed.to_str().unwrap(),
        "predict",
        "--task",
        "disorder",
        "-i",
        fasta.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "a good install must predict");
    assert!(out.contains("probability"), "the table should reach standard output");

    // Burial's weights as the disorder head, metadata rewritten so the file
    // verifies: the install is intact and the numbers are wrong, which is
    // exactly what exit 5 is for.
    let swapped = base.join("swapped");
    copy_models(&installed, &swapped);
    let weights = std::fs::read(swapped.join("burial.safetensors")).unwrap();
    std::fs::write(swapped.join("disorder.safetensors"), &weights).unwrap();
    let path = swapped.join("disorder.metadata.json");
    let mut meta: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    meta["weights_sha256"] = sha256(&weights).into();
    std::fs::write(&path, serde_json::to_string_pretty(&meta).unwrap()).unwrap();

    let (code, out) = run(&["--models", swapped.to_str().unwrap(), "selftest"]);
    assert_eq!(code, 5, "a model that does not reproduce must exit 5");
    assert!(out.contains("FAILED"), "the report should name what failed:\n{out}");
    std::fs::remove_dir_all(&base).ok();
}

/// A corrupt head must stop `models` too, not only a prediction: the whole
/// point of the command is to be the thing an install script runs first.
#[test]
fn models_fails_on_corrupt_weights_and_passes_on_a_good_install() {
    let Some(installed) = models() else {
        println!("no heads installed; skipping.");
        return;
    };
    let base = scratch("models");
    let good = base.join("good");
    copy_models(&installed, &good);
    assert_eq!(run(&["--models", good.to_str().unwrap(), "models"]).0, 0);

    let path = good.join("stability.safetensors");
    let mut bytes = std::fs::read(&path).unwrap();
    bytes.push(0);
    std::fs::write(&path, bytes).unwrap();
    assert_eq!(run(&["--models", good.to_str().unwrap(), "models"]).0, 3);
    std::fs::remove_dir_all(&base).ok();
}
