//! Checks this ESM-2 port against activations captured from the reference
//! `PyTorch` implementation (`transformers`' `EsmForMaskedLM`).
//!
//! The fixtures in `tests/data` were produced once by `scripts/gen_esm2_reference.py`
//! and are committed, so the comparison needs no Python. It does need the
//! `facebook/esm2_t6_8M_UR50D` weights, which are ~31 MB and are downloaded on
//! first run; the tests are therefore `#[ignore]`d by default:
//!
//! ```text
//! cargo test -p opendna-ml --test esm2_reference -- --ignored
//! ```

#![allow(clippy::unwrap_used, missing_docs)]

use opendna_ml::esm2::{Esm2, ModelSize};
use std::path::PathBuf;
use std::sync::{Mutex, MutexGuard, OnceLock, PoisonError};

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn reference() -> &'static serde_json::Value {
    static REF: OnceLock<serde_json::Value> = OnceLock::new();
    REF.get_or_init(|| {
        let text = std::fs::read_to_string(data_dir().join("reference.json")).unwrap();
        serde_json::from_str(&text).unwrap()
    })
}

fn load_f32(name: &str) -> Vec<f32> {
    let bytes = std::fs::read(data_dir().join(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_eq!(bytes.len() % 4, 0, "{name} is not a whole number of f32 values");
    bytes.as_chunks::<4>().0.iter().map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}

fn cases() -> Vec<(String, usize)> {
    reference()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| (c["sequence"].as_str().unwrap().to_string(), c["n_tokens"].as_u64().unwrap() as usize))
        .collect()
}

/// Serialises the tests that run the model.
///
/// `cargo test` runs tests on several threads, and they share one [`Esm2`]
/// through the `OnceLock` below. Driving one CUDA model from several threads
/// at once changes results in the last bits — about 2e-7 relative, which is
/// algorithm selection inside cuBLAS rather than a wrong answer, but it is
/// enough to break the exactness check in
/// `windowing_is_exact_for_short_sequences`. Serialising restores
/// reproducibility instead of hiding the variation behind a tolerance.
///
/// Verified: the whole file passes with `--test-threads=1` and that one test
/// fails without it, while three sequential `embed` calls of the same
/// sequence on the GPU are bit-identical.
fn exclusive() -> MutexGuard<'static, ()> {
    static LOCK: Mutex<()> = Mutex::new(());
    LOCK.lock().unwrap_or_else(PoisonError::into_inner)
}

fn model() -> &'static Esm2 {
    static MODEL: OnceLock<Esm2> = OnceLock::new();
    MODEL.get_or_init(|| {
        let device = opendna_ml::best_device().unwrap();
        Esm2::from_pretrained(ModelSize::Small8M, &device)
            .expect("ESM-2 8M weights; run with network access on first use")
    })
}

/// Largest absolute difference, and where it is.
fn max_diff(a: &[f32], b: &[f32]) -> (f32, usize) {
    assert_eq!(a.len(), b.len(), "length mismatch: {} vs {}", a.len(), b.len());
    a.iter().zip(b).enumerate().map(|(i, (x, y))| ((x - y).abs(), i)).fold((0.0, 0), |acc, d| {
        if d.0 > acc.0 {
            d
        } else {
            acc
        }
    })
}

/// Tolerance for f32 accumulation differences between two implementations.
///
/// The observed maxima on this fixture set are 1.4e-5 (hidden states) and
/// 6.3e-5 (logits) against activations reaching ~17, so this leaves roughly
/// eightfold headroom for platform differences while still catching real
/// defects: loading the rotary frequencies from the formula instead of from
/// the checkpoint, for instance, produces errors of 6.6e-3.
const TOLERANCE: f32 = 5e-4;

#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn hidden_states_match_reference() {
    let _serial = exclusive();
    let model = model();
    for (i, (sequence, n_tokens)) in cases().into_iter().enumerate() {
        let hidden = model.hidden_states(&[&sequence]).unwrap();
        assert_eq!(hidden.dims(), &[1, n_tokens, model.hidden_size()], "case {i} shape");
        let got = hidden.flatten_all().unwrap().to_vec1::<f32>().unwrap();
        let want = load_f32(&format!("case{i}_hidden.f32"));
        let (diff, at) = max_diff(&got, &want);
        assert!(
            diff < TOLERANCE,
            "case {i}: hidden states differ by {diff} at element {at} (got {}, want {})",
            got[at],
            want[at]
        );
    }
}

#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn logits_match_reference() {
    let _serial = exclusive();
    let model = model();
    for (i, (sequence, _)) in cases().into_iter().enumerate() {
        let got = model.logits(&[&sequence]).unwrap().flatten_all().unwrap().to_vec1::<f32>().unwrap();
        let want = load_f32(&format!("case{i}_logits.f32"));
        let (diff, at) = max_diff(&got, &want);
        assert!(diff < TOLERANCE, "case {i}: logits differ by {diff} at element {at}");
    }
}

#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn pooled_embeddings_match_reference() {
    let _serial = exclusive();
    let model = model();
    for (i, (sequence, _)) in cases().into_iter().enumerate() {
        let out = model.embed(&[&sequence]).unwrap();
        let want = load_f32(&format!("case{i}_pooled.f32"));
        let (diff, at) = max_diff(&out[0].pooled, &want);
        assert!(diff < TOLERANCE, "case {i}: pooled embedding differs by {diff} at element {at}");
        assert_eq!(out[0].length, sequence.len());
        assert_eq!(out[0].per_residue.len(), sequence.len() * model.hidden_size());
    }
}

/// Padding must not change a sequence's embedding: batching short and long
/// sequences together has to give the same answer as embedding them alone.
#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn padding_does_not_affect_results() {
    let _serial = exclusive();
    let model = model();
    let all = cases();
    let sequences: Vec<&str> = all.iter().map(|(s, _)| s.as_str()).collect();
    let batched = model.embed(&sequences).unwrap();
    for (i, sequence) in sequences.iter().enumerate() {
        let alone = model.embed(&[sequence]).unwrap();
        let (diff, at) = max_diff(&batched[i].pooled, &alone[0].pooled);
        assert!(diff < TOLERANCE, "case {i}: batching changed the embedding by {diff} at element {at}");
    }
}

#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn rejects_empty_and_overlong_sequences() {
    let _serial = exclusive();
    let model = model();
    assert!(matches!(model.embed(&[""]), Err(opendna_ml::MlError::EmptySequence { index: 0 })));
    let too_long = "A".repeat(model.max_residues() + 1);
    assert!(matches!(model.embed(&[&too_long]), Err(opendna_ml::MlError::SequenceTooLong { index: 0, .. })));
}

/// Windowing must be a no-op for sequences that already fit.
#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn windowing_is_exact_for_short_sequences() {
    let _serial = exclusive();
    let model = model();
    let (sequence, _) = cases().remove(0);
    let direct = model.embed(&[&sequence]).unwrap().remove(0);
    let windowed = model.embed_windowed(&sequence).unwrap();
    assert_eq!(windowed.length, direct.length);
    let (diff, at) = max_diff(&windowed.per_residue, &direct.per_residue);
    assert!(diff == 0.0, "windowing changed a short sequence by {diff} at {at}");
}

/// A sequence past the model's limit still produces one vector per residue,
/// and the vectors are real values rather than the zero-initialised buffer.
#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn windowing_covers_every_residue_of_a_long_sequence() {
    let _serial = exclusive();
    let model = model();
    let unit = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    let long: String = unit.repeat(20); // 1520 residues, past the 1024 limit
    assert!(long.len() > model.max_residues());

    let out = model.embed_windowed(&long).unwrap();
    assert_eq!(out.length, long.len());
    assert_eq!(out.per_residue.len(), long.len() * model.hidden_size());
    for i in [0, 1, long.len() / 2, long.len() - 1] {
        let v = out.residue(i).unwrap();
        assert!(v.iter().any(|x| *x != 0.0), "residue {i} was never written");
        assert!(v.iter().all(|x| x.is_finite()), "residue {i} has non-finite values");
    }
    // Direct embedding of the same sequence must be refused, which is why
    // windowing exists.
    assert!(matches!(model.embed(&[&long]), Err(opendna_ml::MlError::SequenceTooLong { .. })));
}

/// The masked distribution should favour the residue that is actually there:
/// a language model that cannot do this is not usable as a mutation predictor.
#[test]
#[ignore = "downloads ~31 MB of ESM-2 weights on first run"]
fn masked_predictions_recover_the_true_residue() {
    use opendna_ml::esm2::alphabet;
    let _serial = exclusive();
    let model = model();
    let sequence = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    let positions: Vec<usize> = (0..sequence.len()).collect();
    let rows = model.masked_log_probs(sequence, &positions).unwrap();
    assert_eq!(rows.len(), sequence.len());

    let mut recovered = 0;
    for (position, row) in rows.iter().enumerate() {
        assert_eq!(row.len(), alphabet::TOKENS.len());
        let best = row
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i as u32)
            .unwrap();
        let actual = alphabet::token_id(sequence.as_bytes()[position] as char);
        if best == actual {
            recovered += 1;
        }
        // Rows are log-probabilities, so they must be non-positive.
        assert!(row.iter().all(|v| *v <= 1e-4), "position {position} has a positive log-probability");
    }
    // ESM-2 8M is the weakest published checkpoint; on ubiquitin it should
    // still beat the 5% a uniform guess over 20 residues would give.
    assert!(
        recovered * 4 >= sequence.len(),
        "recovered only {recovered}/{} residues, which suggests the masking is wrong",
        sequence.len()
    );
}
