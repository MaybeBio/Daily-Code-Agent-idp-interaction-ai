//! `stability`: single-substitution ddG from Tsuboyama et al. 2023.
//!
//! Canonical convention for every row in the output:
//!
//! ```text
//! ddg = dG_unfold(mutant) - dG_unfold(wild type)   [kcal/mol]
//! ```
//!
//! so a positive value means the substitution is stabilising. Tsuboyama et
//! al. already report in this convention (every row they flag
//! `Stabilizing_mut` has `ddG_ML > 0`, which the pipeline checks); the conversion layer
//! exists so that further sources with other conventions can be added
//! without ambiguity.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, Read};

use anyhow::{bail, Context, Result};

use crate::card::{CardInput, CardLicense};
use crate::cluster::ClusterParams;
use crate::context::{write_dataset, BuildContext, FinishedDataset};
use crate::parquet_out::{ColumnData, Table};
use crate::sources::{self, Source};
use crate::split::{assign_groups, Fractions, Split, SplitSummary};
use crate::{is_standard_aa, is_standard_protein, sequence_id};

pub const NAME: &str = "stability";

/// CSV inside the Zenodo archive holding Tsuboyama's Datasets 2 and 3.
pub const TSUBOYAMA_ENTRY: &str = "Processed_K50_dG_datasets/Tsuboyama2023_Dataset2_Dataset3_20230416.csv";

/// Replicate measurements of the same substitution further apart than this
/// are treated as conflicting and dropped.
pub const MAX_DUPLICATE_SPREAD_KCAL: f64 = 1.0;

const KJ_PER_KCAL: f64 = 4.184;

/// Sign convention used by a source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignConvention {
    /// dG_unfold(mut) - dG_unfold(wt); positive = stabilising (Tsuboyama 2023).
    StabilizingPositive,
    /// dG_unfold(wt) - dG_unfold(mut), equivalently ddG of folding; positive = destabilising
    /// (the convention of FoldX, Rosetta `cartesian_ddg` and many older databases).
    DestabilizingPositive,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyUnit {
    KcalPerMol,
    KjPerMol,
}

/// Converts a source ddG to the canonical convention (kcal/mol, positive = stabilising).
pub fn to_canonical_ddg(value: f64, sign: SignConvention, unit: EnergyUnit) -> f64 {
    let kcal = match unit {
        EnergyUnit::KcalPerMol => value,
        EnergyUnit::KjPerMol => value / KJ_PER_KCAL,
    };
    match sign {
        SignConvention::StabilizingPositive => kcal,
        SignConvention::DestabilizingPositive => -kcal,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Substitution {
    pub wt: u8,
    /// 1-based position in the wild-type sequence.
    pub pos: u32,
    pub mt: u8,
}

impl Substitution {
    pub fn label(&self) -> String {
        format!("{}{}{}", self.wt as char, self.pos, self.mt as char)
    }

    /// Applies the substitution to `wt_seq`, checking that the wild-type residue matches.
    pub fn apply(&self, wt_seq: &str) -> Option<String> {
        let i = (self.pos as usize).checked_sub(1)?;
        let bytes = wt_seq.as_bytes();
        if bytes.get(i) != Some(&self.wt) {
            return None;
        }
        let mut out = bytes.to_vec();
        out[i] = self.mt;
        String::from_utf8(out).ok()
    }
}

/// Parses a single substitution such as `D12Q` (standard amino acids only, 1-based).
pub fn parse_substitution(s: &str) -> Option<Substitution> {
    let b = s.as_bytes();
    if b.len() < 3 {
        return None;
    }
    let (wt, mt) = (b[0], b[b.len() - 1]);
    if !is_standard_aa(wt) || !is_standard_aa(mt) {
        return None;
    }
    let digits = &s[1..s.len() - 1];
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) || digits.starts_with('0') {
        return None;
    }
    Some(Substitution { wt, pos: digits.parse().ok()?, mt })
}

#[derive(Debug, Clone, PartialEq)]
pub struct Measurement {
    pub protein_id: String,
    pub wt_sequence: String,
    pub sub: Substitution,
    /// Canonical ddG (kcal/mol, positive = stabilising).
    pub ddg: f64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadStats {
    pub rows: u64,
    pub wild_type_rows: u64,
    pub not_single_substitution: u64,
    pub outside_ml_subset: u64,
    pub synonymous: u64,
    pub missing_wild_type: u64,
    pub nonstandard_wild_type: u64,
    pub sequence_mismatch: u64,
    pub kept: u64,
}

fn col(headers: &csv::StringRecord, name: &str) -> Result<usize> {
    headers.iter().position(|h| h == name).with_context(|| format!("Tsuboyama CSV has no column {name:?}"))
}

/// Reads Tsuboyama Dataset 2/3 rows and keeps single substitutions with a
/// value in the ML-grade subset (Dataset 3: numeric `ddG_ML` and `dG_ML`).
pub fn read_tsuboyama(r: impl Read) -> Result<(Vec<Measurement>, ReadStats)> {
    let mut rdr = csv::ReaderBuilder::new().has_headers(true).from_reader(r);
    let h = rdr.headers()?.clone();
    let (c_seq, c_mut, c_wt, c_ddg, c_dg, c_stab) = (
        col(&h, "aa_seq")?,
        col(&h, "mut_type")?,
        col(&h, "WT_name")?,
        col(&h, "ddG_ML")?,
        col(&h, "dG_ML")?,
        col(&h, "Stabilizing_mut")?,
    );
    let mut stats = ReadStats::default();
    let mut wild_types: HashMap<String, String> = HashMap::new();
    let mut conflicting_wt: Vec<String> = Vec::new();
    // (wt name, substitution, mutant sequence, ddg)
    let mut pending: Vec<(String, Substitution, String, f64)> = Vec::new();

    for rec in rdr.records() {
        let rec = rec?;
        stats.rows += 1;
        let mut_type = &rec[c_mut];
        let wt_name = &rec[c_wt];
        if mut_type == "wt" {
            stats.wild_type_rows += 1;
            let seq = rec[c_seq].to_owned();
            match wild_types.get(wt_name) {
                Some(prev) if *prev != seq => conflicting_wt.push(wt_name.to_owned()),
                Some(_) => {}
                None => {
                    wild_types.insert(wt_name.to_owned(), seq);
                }
            }
            continue;
        }
        let Some(sub) = parse_substitution(mut_type) else {
            stats.not_single_substitution += 1;
            continue;
        };
        let (Ok(ddg), Ok(_dg)) = (rec[c_ddg].parse::<f64>(), rec[c_dg].parse::<f64>()) else {
            stats.outside_ml_subset += 1;
            continue;
        };
        if !ddg.is_finite() {
            stats.outside_ml_subset += 1;
            continue;
        }
        // Consistency check of the sign convention against the source's own flag.
        // `Stabilizing_mut` is a thresholded call, so only the True direction is a hard
        // implication (True => ddG > 0); False rows may have small positive ddG.
        if &rec[c_stab] == "True" && ddg <= 0.0 {
            bail!("row {mut_type} of {wt_name}: Stabilizing_mut=True but ddG_ML={ddg}; sign convention changed?");
        }
        if sub.wt == sub.mt {
            stats.synonymous += 1;
            continue;
        }
        pending.push((wt_name.to_owned(), sub, rec[c_seq].to_owned(), ddg));
    }
    for name in &conflicting_wt {
        wild_types.remove(name);
    }

    let mut out = Vec::with_capacity(pending.len());
    for (wt_name, sub, mut_seq, ddg) in pending {
        let Some(wt_seq) = wild_types.get(&wt_name) else {
            stats.missing_wild_type += 1;
            continue;
        };
        if !is_standard_protein(wt_seq) {
            stats.nonstandard_wild_type += 1;
            continue;
        }
        if sub.apply(wt_seq).as_deref() != Some(mut_seq.as_str()) {
            stats.sequence_mismatch += 1;
            continue;
        }
        out.push(Measurement {
            protein_id: wt_name,
            wt_sequence: wt_seq.clone(),
            sub,
            ddg: to_canonical_ddg(ddg, SignConvention::StabilizingPositive, EnergyUnit::KcalPerMol),
        });
    }
    stats.kept = out.len() as u64;
    Ok((out, stats))
}

#[derive(Debug, Clone, PartialEq)]
pub struct Merged {
    pub protein_id: String,
    pub wt_sequence: String,
    pub sub: Substitution,
    pub ddg: f64,
    pub n_measurements: u32,
    pub spread: f64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MergeStats {
    pub duplicate_groups: u64,
    pub conflicting_groups: u64,
    pub conflicting_measurements: u64,
}

/// Collapses measurements of the same substitution in the same wild-type
/// sequence. Consistent replicates (spread <= `max_spread`) are averaged;
/// conflicting ones are dropped entirely rather than resolved arbitrarily.
pub fn merge_duplicates(ms: Vec<Measurement>, max_spread: f64) -> (Vec<Merged>, MergeStats) {
    let mut groups: BTreeMap<(String, u32, u8), Vec<Measurement>> = BTreeMap::new();
    for m in ms {
        groups.entry((m.wt_sequence.clone(), m.sub.pos, m.sub.mt)).or_default().push(m);
    }
    let mut stats = MergeStats::default();
    let mut out = Vec::with_capacity(groups.len());
    for (_, g) in groups {
        let n = g.len();
        let min = g.iter().map(|m| m.ddg).fold(f64::INFINITY, f64::min);
        let max = g.iter().map(|m| m.ddg).fold(f64::NEG_INFINITY, f64::max);
        if n > 1 {
            stats.duplicate_groups += 1;
            if max - min > max_spread {
                stats.conflicting_groups += 1;
                stats.conflicting_measurements += n as u64;
                continue;
            }
        }
        let mean = g.iter().map(|m| m.ddg).sum::<f64>() / n as f64;
        let first = g.iter().min_by(|a, b| a.protein_id.cmp(&b.protein_id)).expect("group is non-empty");
        out.push(Merged {
            protein_id: first.protein_id.clone(),
            wt_sequence: first.wt_sequence.clone(),
            sub: first.sub,
            ddg: mean,
            n_measurements: n as u32,
            spread: max - min,
        });
    }
    (out, stats)
}

pub fn to_table(rows: &[&Merged], cluster_of: &BTreeMap<String, String>) -> Table {
    let mut t = Table::new();
    t.push(
        "protein_id",
        "Tsuboyama 2023 wild-type name (WT_name); for merged duplicates, the lexicographically smallest",
        ColumnData::Utf8(rows.iter().map(|r| r.protein_id.clone()).collect()),
    );
    t.push(
        "wt_sequence",
        "Wild-type domain sequence (aa_seq, without the display linkers)",
        ColumnData::Utf8(rows.iter().map(|r| r.wt_sequence.clone()).collect()),
    );
    t.push(
        "mutation",
        "Substitution as <wt><1-based position><mut>, e.g. D12Q",
        ColumnData::Utf8(rows.iter().map(|r| r.sub.label()).collect()),
    );
    t.push(
        "position",
        "1-based position in wt_sequence",
        ColumnData::Int32(rows.iter().map(|r| r.sub.pos as i32).collect()),
    );
    t.push(
        "wt_aa",
        "Wild-type residue (one-letter code)",
        ColumnData::Utf8(rows.iter().map(|r| (r.sub.wt as char).to_string()).collect()),
    );
    t.push(
        "mut_aa",
        "Mutant residue (one-letter code)",
        ColumnData::Utf8(rows.iter().map(|r| (r.sub.mt as char).to_string()).collect()),
    );
    t.push(
        "ddg",
        "dG_unfold(mutant) - dG_unfold(wild type) in kcal/mol; positive = stabilising. Mean over consistent replicates",
        ColumnData::Float64(rows.iter().map(|r| r.ddg).collect()),
    );
    t.push(
        "n_measurements",
        "Number of source measurements averaged into ddg",
        ColumnData::Int32(rows.iter().map(|r| r.n_measurements as i32).collect()),
    );
    t.push(
        "ddg_spread",
        "max - min of the averaged measurements in kcal/mol (0 for a single measurement)",
        ColumnData::Float64(rows.iter().map(|r| r.spread).collect()),
    );
    t.push(
        "source",
        "Upstream source id",
        ColumnData::Utf8(rows.iter().map(|_| "tsuboyama2023".to_owned()).collect()),
    );
    t.push(
        "cluster_id",
        "MMseqs2 30% identity cluster of the wild-type sequence (representative sequence hash); the split unit",
        ColumnData::Utf8(
            rows.iter()
                .map(|r| cluster_of[&sequence_id(&r.wt_sequence)].clone())
                .collect(),
        ),
    );
    t
}

pub const SOURCES: &[Source] = &[sources::TSUBOYAMA_2023];

pub fn build(ctx: &mut BuildContext) -> Result<()> {
    let zip_path = ctx.cache.fetch(&sources::TSUBOYAMA_2023)?;
    let mut archive = zip::ZipArchive::new(BufReader::new(File::open(&zip_path)?))?;
    let entry = archive
        .by_name(TSUBOYAMA_ENTRY)
        .with_context(|| format!("{TSUBOYAMA_ENTRY} not found in {}", zip_path.display()))?;
    eprintln!("reading {TSUBOYAMA_ENTRY}");
    let (measurements, rs) = read_tsuboyama(BufReader::with_capacity(1 << 20, entry))?;
    let (merged, ms) = merge_duplicates(measurements, MAX_DUPLICATE_SPREAD_KCAL);

    let seqs: BTreeMap<String, String> =
        merged.iter().map(|m| (sequence_id(&m.wt_sequence), m.wt_sequence.clone())).collect();
    let params = ClusterParams::ID30;
    let clustering = ctx.cluster(&seqs, params)?;
    let cluster_of: BTreeMap<String, String> = clustering.member_to_rep.clone();

    let mut weights: BTreeMap<String, u64> = BTreeMap::new();
    for m in &merged {
        *weights.entry(cluster_of[&sequence_id(&m.wt_sequence)].clone()).or_default() += 1;
    }
    let fractions = Fractions::new(0.8, 0.1, 0.1)?;
    let assignment = assign_groups(&weights, &BTreeMap::new(), fractions, ctx.seed)?;

    let mut sorted: Vec<&Merged> = merged.iter().collect();
    sorted.sort_by(|a, b| (&a.protein_id, a.sub.pos, a.sub.mt).cmp(&(&b.protein_id, b.sub.pos, b.sub.mt)));
    let split_of = |m: &Merged| assignment[&cluster_of[&sequence_id(&m.wt_sequence)]];
    let summary = SplitSummary::from_rows(
        sorted.iter().map(|m| cluster_of[&sequence_id(&m.wt_sequence)].as_str()),
        &assignment,
    );
    let tables: Vec<(Split, Table)> = Split::ALL
        .iter()
        .map(|s| {
            let rows: Vec<&Merged> = sorted.iter().copied().filter(|m| split_of(m) == *s).collect();
            (*s, to_table(&rows, &cluster_of))
        })
        .collect();
    let n_proteins = seqs.len();
    let statistics = vec![
        ("CSV rows read".into(), rs.rows.to_string()),
        ("wild-type rows".into(), rs.wild_type_rows.to_string()),
        (
            "dropped: not a single substitution (double mutants, insertions, deletions)".into(),
            rs.not_single_substitution.to_string(),
        ),
        (
            "dropped: outside Tsuboyama's ML-grade subset (ddG_ML or dG_ML not numeric)".into(),
            rs.outside_ml_subset.to_string(),
        ),
        ("dropped: synonymous".into(), rs.synonymous.to_string()),
        ("dropped: wild type missing or ambiguous".into(), rs.missing_wild_type.to_string()),
        ("dropped: non-standard residue in wild type".into(), rs.nonstandard_wild_type.to_string()),
        ("dropped: mutant sequence inconsistent with substitution".into(), rs.sequence_mismatch.to_string()),
        ("measurements kept".into(), rs.kept.to_string()),
        ("duplicate substitution groups".into(), ms.duplicate_groups.to_string()),
        (
            format!("dropped: conflicting duplicate groups (spread > {MAX_DUPLICATE_SPREAD_KCAL} kcal/mol)"),
            format!("{} groups / {} measurements", ms.conflicting_groups, ms.conflicting_measurements),
        ),
        ("final rows".into(), merged.len().to_string()),
        ("distinct wild-type sequences".into(), n_proteins.to_string()),
        ("MMseqs2 clusters".into(), clustering.n_clusters().to_string()),
    ];
    let construction = format!(
        "1. Read `{TSUBOYAMA_ENTRY}` from the Zenodo archive.\n\
         2. Keep rows whose `mut_type` is a single substitution of standard amino acids. \
            Double mutants, insertions and deletions are dropped.\n\
         3. Keep rows in Tsuboyama's ML-grade subset (Dataset 3): numeric `ddG_ML` and `dG_ML`. \
            Rows whose dG was outside the assay's dynamic range (`<-1`, `>5`) or was flagged unreliable (`-`) are dropped.\n\
         4. Check the sign convention against the source's own `Stabilizing_mut` flag (the build aborts on a contradiction), \
            and convert to kcal/mol with positive = stabilising (a no-op for this source).\n\
         5. Verify that each mutant sequence equals the wild-type sequence with the substitution applied; drop rows that do not.\n\
         6. Merge measurements of the same substitution in the same wild-type sequence: average when their spread is at most \
            {MAX_DUPLICATE_SPREAD_KCAL} kcal/mol, otherwise drop all of them.\n\
         7. Cluster the distinct wild-type sequences with `{}`.\n\
         8. Assign whole clusters to train/validation/test (target 80/10/10 of rows) with a seeded deterministic procedure (seed {}). \
            All mutations of one domain, and of all its homologs at 30% identity, are in the same split.",
        params.describe(),
        ctx.seed
    );
    let card = CardInput {
        pretty_name: "OpenDNA Stability (Tsuboyama 2023 single substitutions)",
        license: CardLicense::Known("cc-by-4.0"),
        license_note: "The only included source is licensed CC BY 4.0, and this derived dataset is released under CC BY 4.0. \
                       Attribution: Tsuboyama et al. 2023 (citation below).",
        tags: &["biology", "protein", "protein-stability", "ddg", "mutation-effect"],
        summary: format!(
            "Single-substitution folding stability changes (ddG, kcal/mol) for {n_proteins} small natural and designed \
             protein domains, measured by cDNA display proteolysis (Tsuboyama et al. 2023). Prepared for training a ddG \
             regression head on frozen protein language model embeddings. Positive ddG means the substitution is stabilising. \
             Splits are by 30% sequence identity cluster, so test domains have no close homolog in train."
        ),
        sources: SOURCES,
        excluded: sources::EXCLUDED_SOURCES
            .iter()
            .filter(|e| e.dataset == NAME)
            .collect(),
        construction,
        schema: &tables[0].1,
        row_unit: "mutations",
        group_unit: "clusters",
        splits: &summary,
        cohorts: &[],
        statistics,
        citations: TSUBOYAMA_BIBTEX,
        limitations: &[
            "ddG is inferred from protease resistance through a kinetic model, not measured calorimetrically. The assay's \
             dynamic range is roughly -1 to 5 kcal/mol for dG; mutations that push dG outside it are excluded, which \
             truncates the distribution of strongly destabilising mutations.",
            "Wild-type domains are 32-74 residues long and many are de novo designs. Performance on large, multi-domain or membrane proteins \
             cannot be inferred from this data.",
            "Only single substitutions are included.",
            "FireProtDB was requested as a second source but is excluded for licensing reasons (see above), so all values \
             come from one assay in one laboratory; assay-specific biases are not averaged out.",
            "Clustering uses 30% identity over 80% bidirectional coverage. Homology below that threshold, or shared \
             structural motifs, can still connect train and test.",
            "The duplicate-conflict threshold (1 kcal/mol) is a pipeline choice, not a property of the data.",
        ],
        lock: ctx.cache.lock(),
    };
    write_dataset(
        ctx,
        FinishedDataset {
            name: NAME,
            tables: &tables,
            extra: &[],
            card,
            fractions: Some(fractions),
            cluster_params: Some(params),
        },
    )?;
    Ok(())
}

const TSUBOYAMA_BIBTEX: &str = r#"@article{tsuboyama2023megascale,
  title   = {Mega-scale experimental analysis of protein folding stability in biology and design},
  author  = {Tsuboyama, Kotaro and Dauparas, Justas and Chen, Jonathan and Laine, Elodie and
             Mohseni Behbahani, Yasser and Weinstein, Jonathan J. and Mangan, Niall M. and
             Ovchinnikov, Sergey and Rocklin, Gabriel J.},
  journal = {Nature},
  volume  = {620},
  pages   = {434--444},
  year    = {2023},
  doi     = {10.1038/s41586-023-06328-6}
}

@dataset{tsuboyama2023zenodo,
  title     = {Mega-scale experimental analysis of protein folding stability in biology and design},
  author    = {Tsuboyama, Kotaro and others},
  publisher = {Zenodo},
  year      = {2023},
  doi       = {10.5281/zenodo.7992926}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = include_str!("../../tests/fixtures/tsuboyama_mini.csv");

    #[test]
    fn parses_substitutions() {
        assert_eq!(parse_substitution("D12Q"), Some(Substitution { wt: b'D', pos: 12, mt: b'Q' }));
        for bad in ["wt", "insG1", "delD1", "T9Q:N23Q", "D0Q", "D012Q", "X5A", "D5", "d5q"] {
            assert_eq!(parse_substitution(bad), None, "{bad}");
        }
    }

    #[test]
    fn applies_substitution_with_wt_check() {
        let s = parse_substitution("K2A").unwrap();
        assert_eq!(s.apply("MKV").as_deref(), Some("MAV"));
        assert_eq!(parse_substitution("V2A").unwrap().apply("MKV"), None);
        assert_eq!(parse_substitution("K9A").unwrap().apply("MKV"), None);
    }

    #[test]
    fn normalises_sign_and_units() {
        use EnergyUnit::*;
        use SignConvention::*;
        assert_eq!(to_canonical_ddg(1.5, StabilizingPositive, KcalPerMol), 1.5);
        assert_eq!(to_canonical_ddg(1.5, DestabilizingPositive, KcalPerMol), -1.5);
        assert!((to_canonical_ddg(4.184, StabilizingPositive, KjPerMol) - 1.0).abs() < 1e-12);
        assert!((to_canonical_ddg(8.368, DestabilizingPositive, KjPerMol) + 2.0).abs() < 1e-12);
    }

    #[test]
    fn reads_fixture_with_expected_filtering() {
        let (ms, st) = read_tsuboyama(FIXTURE.as_bytes()).unwrap();
        assert_eq!(st.rows, 11);
        assert_eq!(st.wild_type_rows, 3);
        assert_eq!(st.not_single_substitution, 2); // double + insertion
        assert_eq!(st.outside_ml_subset, 2); // ddG_ML "-" and dG_ML ">5"
        assert_eq!(st.sequence_mismatch, 1); // mutant aa_seq does not match
        assert_eq!(st.missing_wild_type, 0);
        assert_eq!(st.kept, 3);
        let labels: Vec<String> = ms.iter().map(|m| m.sub.label()).collect();
        assert_eq!(labels, ["K2A", "V3G", "K2A"]);
        assert!(ms.iter().any(|m| m.ddg > 0.0));
    }

    #[test]
    fn contradicting_sign_flag_aborts() {
        let bad = FIXTURE.replace("0.40,True", "-0.40,True");
        assert!(read_tsuboyama(bad.as_bytes()).is_err());
    }

    fn m(id: &str, seq: &str, sub: &str, ddg: f64) -> Measurement {
        Measurement {
            protein_id: id.into(),
            wt_sequence: seq.into(),
            sub: parse_substitution(sub).unwrap(),
            ddg,
        }
    }

    #[test]
    fn merges_consistent_and_drops_conflicting_duplicates() {
        let input = vec![
            m("b", "MKV", "K2A", 1.0),
            m("a", "MKV", "K2A", 1.4), // same sequence + substitution, different name
            m("a", "MKV", "V3G", -2.0),
            m("a", "MKV", "V3G", 0.5), // conflicting: spread 2.5
            m("c", "MKL", "K2A", 0.3),
        ];
        let (out, st) = merge_duplicates(input, 1.0);
        assert_eq!(st.duplicate_groups, 2);
        assert_eq!(st.conflicting_groups, 1);
        assert_eq!(st.conflicting_measurements, 2);
        assert_eq!(out.len(), 2);
        let k2a = out.iter().find(|r| r.wt_sequence == "MKV").unwrap();
        assert_eq!(k2a.protein_id, "a");
        assert_eq!(k2a.n_measurements, 2);
        assert!((k2a.ddg - 1.2).abs() < 1e-12);
        assert!((k2a.spread - 0.4).abs() < 1e-12);
    }
}
