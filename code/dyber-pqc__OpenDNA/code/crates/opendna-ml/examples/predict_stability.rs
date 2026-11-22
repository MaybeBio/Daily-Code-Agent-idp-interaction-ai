//! Scores every substitution at one position of a sequence.
//!
//! The position is 1-based, matching how mutations are written (`L8P` is
//! position 8) and matching the desktop application's own API.
//!
//! ```text
//! cargo run --release -p opendna-ml --example predict_stability -- models MQIFVKTLTGK 8
//! ```
//!
//! Needs a trained `stability` head in the given directory; see `training/`.

#![allow(clippy::print_stdout, clippy::unwrap_used)]

use opendna_core::AminoAcid;
use opendna_ml::esm2::{Esm2, ModelSize};
use opendna_ml::heads::{Heads, StabilityHead};
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let dir = PathBuf::from(args.next().unwrap_or_else(|| "models".into()));
    let sequence = args.next().unwrap_or_else(|| {
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG".into()
    });
    let position: usize = args.next().map_or(Ok(8), |s| s.parse())?;
    let residues: Vec<char> = sequence.chars().collect();
    if position == 0 || position > residues.len() {
        println!(
            "Position {position} is outside a {}-residue sequence; positions are 1-based.",
            residues.len()
        );
        return Ok(());
    }
    let index = position - 1;

    let embedder = ModelSize::Large150M;
    let device = opendna_ml::best_device()?;
    println!("Loading on {}", opendna_ml::device_description(&device));

    let heads = Heads::load(&dir, embedder, &device)?;
    let Some(head) = heads.stability.as_ref() else {
        println!("No stability head in {}. Train one first; see training/README.md.", dir.display());
        return Ok(());
    };

    let model = Esm2::from_pretrained(embedder, &device)?;
    let embedding = model.embed(&[&sequence])?.remove(0);
    let log_probs = model.masked_log_probs(&sequence, &[index])?.remove(0);

    let wild_type = residues[index];
    let site = embedding.residue(index).unwrap();
    let wt_index = AminoAcid::STANDARD.iter().position(|a| a.one_letter() == wild_type).unwrap();
    let wt_logp = log_probs[opendna_ml::esm2::alphabet::token_id(wild_type) as usize];

    // One batch through the head, via the same helper the desktop
    // application uses, so the two cannot produce different numbers.
    let mut rows = Vec::new();
    let mut mutants = Vec::new();
    for (mutant_index, aa) in AminoAcid::STANDARD.iter().enumerate() {
        if mutant_index == wt_index {
            continue;
        }
        let logp = log_probs[opendna_ml::esm2::alphabet::token_id(aa.one_letter()) as usize];
        rows.push(StabilityHead::features(site, wt_index, mutant_index, logp - wt_logp));
        mutants.push(aa.one_letter());
    }
    let ddgs = head.predict(&rows, model.hidden_size(), &device)?;
    let mut scored: Vec<(char, f32)> = mutants.into_iter().zip(ddgs).collect();
    scored.sort_by(|a, b| a.1.total_cmp(&b.1));

    println!("\nPredicted ddG at position {position} (wild type {wild_type}), kcal/mol");
    println!("Negative is destabilising. Most to least destabilising:\n");
    for (residue, ddg) in &scored {
        println!("  {wild_type}{position}{residue}  {ddg:+.2}");
    }
    Ok(())
}
