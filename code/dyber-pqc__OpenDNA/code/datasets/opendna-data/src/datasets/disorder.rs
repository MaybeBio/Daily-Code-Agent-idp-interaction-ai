//! `disorder`: per-residue intrinsic disorder labels from DisProt, with the
//! CAID3 Disorder-PDB proteins (and everything clustered with them) held out
//! as the test split.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use anyhow::{bail, Context, Result};
use serde::Deserialize;

use crate::card::{CardInput, CardLicense};
use crate::cluster::ClusterParams;
use crate::context::{write_dataset, BuildContext, FinishedDataset};
use crate::is_standard_protein;
use crate::parquet_out::{ColumnData, Table};
use crate::sources::{self, Source};
use crate::split::{assign_groups, Fractions, Split, SplitSummary};

pub const NAME: &str = "disorder";

/// IDPO term for "disorder" in the DisProt Structural state namespace.
pub const IDPO_DISORDER: &str = "IDPO:0000002";
/// Structural-state terms that describe partially folded states; counted but not labelled.
pub const IDPO_MOLTEN: &[&str] = &["IDPO:0000003", "IDPO:0000004"];

const CAID_PREFIX: &str = "caid_";

#[derive(Debug, Deserialize)]
pub struct DisprotRelease {
    pub data: Vec<DisprotEntry>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct DisprotEntry {
    pub disprot_id: String,
    pub acc: String,
    pub sequence: String,
    pub length: usize,
    #[serde(default)]
    pub organism: Option<String>,
    #[serde(default)]
    pub ncbi_taxon_id: Option<i64>,
    #[serde(default)]
    pub regions: Vec<DisprotRegion>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct DisprotRegion {
    #[serde(default)]
    pub term_id: Option<String>,
    pub start: i64,
    pub end: i64,
}

pub fn read_disprot(r: impl Read) -> Result<Vec<DisprotEntry>> {
    let rel: DisprotRelease = serde_json::from_reader(r).context("parsing DisProt JSON")?;
    let ids: BTreeSet<&str> = rel.data.iter().map(|e| e.disprot_id.as_str()).collect();
    if ids.len() != rel.data.len() {
        bail!("DisProt release contains duplicate disprot_id values");
    }
    Ok(rel.data)
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LabelStats {
    pub disorder_regions: u64,
    pub out_of_range_regions: u64,
    pub molten_globule_regions: u64,
}

/// Per-residue labels: 1 inside any `IDPO:0000002` region, 0 elsewhere.
///
/// DisProt coordinates are 1-based and inclusive. Regions that fall outside
/// the sequence are skipped and counted rather than clipped.
pub fn disorder_labels(len: usize, regions: &[DisprotRegion], stats: &mut LabelStats) -> Vec<i8> {
    let mut labels = vec![0i8; len];
    for r in regions {
        let term = r.term_id.as_deref().unwrap_or("");
        if IDPO_MOLTEN.contains(&term) {
            stats.molten_globule_regions += 1;
            continue;
        }
        if term != IDPO_DISORDER {
            continue;
        }
        if r.start < 1 || r.end < r.start || r.end as usize > len {
            stats.out_of_range_regions += 1;
            continue;
        }
        stats.disorder_regions += 1;
        for l in &mut labels[(r.start as usize - 1)..(r.end as usize)] {
            *l = 1;
        }
    }
    labels
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CaidEntry {
    pub id: String,
    pub sequence: String,
}

/// Parses a CAID reference file: three lines per protein (`>id`, sequence,
/// labels in `{0,1,-}`). The labels are validated for shape but not used.
pub fn read_caid_reference(r: impl BufRead) -> Result<Vec<CaidEntry>> {
    let lines: Vec<String> = r
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .filter(|l| !l.trim().is_empty())
        .collect();
    if !lines.len().is_multiple_of(3) {
        bail!("CAID reference has {} non-empty lines, expected a multiple of 3", lines.len());
    }
    let mut out = Vec::with_capacity(lines.len() / 3);
    for chunk in lines.chunks(3) {
        let Some(id) = chunk[0].strip_prefix('>') else {
            bail!("expected a CAID header line, got {:?}", chunk[0]);
        };
        let seq = chunk[1].trim();
        let labels = chunk[2].trim();
        if seq.len() != labels.len() {
            bail!("CAID entry {id}: sequence and label lengths differ");
        }
        if !labels.bytes().all(|c| matches!(c, b'0' | b'1' | b'-')) {
            bail!("CAID entry {id}: unexpected label characters");
        }
        out.push(CaidEntry {
            id: id.split_whitespace().next().unwrap_or("").to_owned(),
            sequence: seq.to_owned(),
        });
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct Protein {
    pub entry: DisprotEntry,
    pub labels: Vec<i8>,
    pub caid3: bool,
}

#[derive(Debug, Clone, Default)]
pub struct CleanStats {
    pub entries: u64,
    pub nonstandard: u64,
    pub length_mismatch: u64,
    pub without_disorder: u64,
    pub labels: LabelStats,
}

pub fn clean(entries: Vec<DisprotEntry>, caid_ids: &BTreeSet<String>) -> (Vec<Protein>, CleanStats) {
    let mut st = CleanStats::default();
    let mut out = Vec::new();
    for e in entries {
        st.entries += 1;
        if !is_standard_protein(&e.sequence) {
            st.nonstandard += 1;
            continue;
        }
        if e.sequence.len() != e.length {
            st.length_mismatch += 1;
            continue;
        }
        let labels = disorder_labels(e.sequence.len(), &e.regions, &mut st.labels);
        if !labels.contains(&1) {
            st.without_disorder += 1;
        }
        let caid3 = caid_ids.contains(&e.disprot_id);
        out.push(Protein { entry: e, labels, caid3 });
    }
    out.sort_by(|a, b| a.entry.disprot_id.cmp(&b.entry.disprot_id));
    (out, st)
}

/// Clusters that must go to test: any cluster containing a CAID protein,
/// identified either by its DisProt id or by its CAID reference sequence.
pub fn caid_clusters(
    member_to_rep: &BTreeMap<String, String>,
    caid_ids: &BTreeSet<String>,
) -> BTreeSet<String> {
    member_to_rep
        .iter()
        .filter(|(m, _)| m.starts_with(CAID_PREFIX) || caid_ids.contains(*m))
        .map(|(_, rep)| rep.clone())
        .collect()
}

pub fn to_table(rows: &[&Protein], cluster_of: &BTreeMap<String, String>) -> Table {
    let mut t = Table::new();
    t.push(
        "disprot_id",
        "DisProt accession",
        ColumnData::Utf8(rows.iter().map(|p| p.entry.disprot_id.clone()).collect()),
    );
    t.push(
        "uniprot_acc",
        "UniProt accession of the DisProt entry",
        ColumnData::Utf8(rows.iter().map(|p| p.entry.acc.clone()).collect()),
    );
    t.push(
        "organism",
        "Source organism as given by DisProt",
        ColumnData::OptUtf8(rows.iter().map(|p| p.entry.organism.clone()).collect()),
    );
    t.push(
        "ncbi_taxon_id",
        "NCBI taxonomy id",
        ColumnData::OptInt64(rows.iter().map(|p| p.entry.ncbi_taxon_id).collect()),
    );
    t.push(
        "sequence",
        "Full UniProt sequence used by DisProt",
        ColumnData::Utf8(rows.iter().map(|p| p.entry.sequence.clone()).collect()),
    );
    t.push(
        "length",
        "Sequence length",
        ColumnData::Int32(rows.iter().map(|p| p.entry.sequence.len() as i32).collect()),
    );
    t.push(
        "disorder",
        "Per-residue label: 1 = inside a DisProt IDPO:0000002 (disorder) region, 0 = not annotated as disordered (not evidence of order)",
        ColumnData::ListInt8(rows.iter().map(|p| p.labels.clone()).collect()),
    );
    t.push(
        "n_disordered",
        "Number of residues labelled 1",
        ColumnData::Int32(rows.iter().map(|p| p.labels.iter().filter(|l| **l == 1).count() as i32).collect()),
    );
    t.push(
        "caid3_target",
        "True if the protein is in the CAID3 Disorder-PDB reference set",
        ColumnData::Bool(rows.iter().map(|p| p.caid3).collect()),
    );
    t.push(
        "cluster_id",
        "MMseqs2 30% identity cluster (representative id); the split unit",
        ColumnData::Utf8(rows.iter().map(|p| cluster_of[&p.entry.disprot_id].clone()).collect()),
    );
    t
}

pub const SOURCES: &[Source] = &[sources::DISPROT_2025_12, sources::CAID3_DISORDER_PDB];

pub fn build(ctx: &mut BuildContext) -> Result<()> {
    let dp_path = ctx.cache.fetch(&sources::DISPROT_2025_12)?;
    let caid_path = ctx.cache.fetch(&sources::CAID3_DISORDER_PDB)?;
    let entries = read_disprot(BufReader::with_capacity(1 << 20, File::open(&dp_path)?))?;
    let caid = read_caid_reference(BufReader::new(File::open(&caid_path)?))?;
    let caid_ids: BTreeSet<String> = caid.iter().map(|c| c.id.clone()).collect();

    let (proteins, st) = clean(entries, &caid_ids);
    let present: BTreeSet<&str> =
        proteins.iter().filter(|p| p.caid3).map(|p| p.entry.disprot_id.as_str()).collect();

    let mut seqs: BTreeMap<String, String> =
        proteins.iter().map(|p| (p.entry.disprot_id.clone(), p.entry.sequence.clone())).collect();
    // CAID sequences are clustered alongside DisProt so that homologs of CAID
    // targets are caught even when the DisProt sequence differs from the CAID one.
    for c in &caid {
        if is_standard_protein(&c.sequence) {
            seqs.insert(format!("{CAID_PREFIX}{}", c.id), c.sequence.clone());
        }
    }
    let params = ClusterParams::ID30;
    let clustering = ctx.cluster(&seqs, params)?;
    let cluster_of = clustering.member_to_rep.clone();

    let test_clusters = caid_clusters(&cluster_of, &caid_ids);
    let mut weights: BTreeMap<String, u64> = cluster_of.values().map(|r| (r.clone(), 0)).collect();
    for p in &proteins {
        *weights.get_mut(&cluster_of[&p.entry.disprot_id]).expect("clustered") +=
            p.entry.sequence.len() as u64;
    }
    let fixed: BTreeMap<String, Split> = test_clusters.iter().map(|c| (c.clone(), Split::Test)).collect();
    let fractions = Fractions::new(0.9, 0.1, 0.0)?;
    let assignment = assign_groups(&weights, &fixed, fractions, ctx.seed)?;

    let split_of = |p: &Protein| assignment[&cluster_of[&p.entry.disprot_id]];
    let homolog_in_test = proteins.iter().filter(|p| !p.caid3 && split_of(p) == Split::Test).count();
    let summary = SplitSummary::from_rows(
        proteins.iter().map(|p| cluster_of[&p.entry.disprot_id].as_str()),
        &assignment,
    );
    let tables: Vec<(Split, Table)> = Split::ALL
        .iter()
        .map(|s| {
            let rows: Vec<&Protein> = proteins.iter().filter(|p| split_of(p) == *s).collect();
            (*s, to_table(&rows, &cluster_of))
        })
        .collect();
    let residues = |s: Split| -> u64 {
        proteins.iter().filter(|p| split_of(p) == s).map(|p| p.labels.len() as u64).sum()
    };
    let disordered: u64 = proteins.iter().map(|p| p.labels.iter().filter(|l| **l == 1).count() as u64).sum();
    let total_res: u64 = proteins.iter().map(|p| p.labels.len() as u64).sum();
    let statistics = vec![
        ("DisProt 2025_12 entries".into(), st.entries.to_string()),
        ("dropped: non-standard residues in sequence".into(), st.nonstandard.to_string()),
        ("dropped: sequence length disagrees with `length`".into(), st.length_mismatch.to_string()),
        ("IDPO:0000002 regions used".into(), st.labels.disorder_regions.to_string()),
        ("regions skipped: outside the sequence".into(), st.labels.out_of_range_regions.to_string()),
        ("regions not labelled: (pre-)molten globule".into(), st.labels.molten_globule_regions.to_string()),
        ("proteins with no IDPO:0000002 residue (kept, all zeros)".into(), st.without_disorder.to_string()),
        ("final proteins".into(), proteins.len().to_string()),
        ("residues / labelled disordered".into(), format!("{total_res} / {disordered}")),
        ("CAID3 Disorder-PDB proteins".into(), caid_ids.len().to_string()),
        ("CAID3 proteins present in this DisProt release".into(), present.len().to_string()),
        (
            "non-CAID proteins moved to test (clustered with a CAID3 protein)".into(),
            homolog_in_test.to_string(),
        ),
        ("MMseqs2 clusters (including CAID-only)".into(), clustering.n_clusters().to_string()),
        (
            "residues train / validation / test".into(),
            format!(
                "{} / {} / {}",
                residues(Split::Train),
                residues(Split::Validation),
                residues(Split::Test)
            ),
        ),
    ];
    let construction = format!(
        "1. Download DisProt release 2025_12 through the DisProt API with ambiguous and obsolete annotations excluded.\n\
         2. Drop entries whose sequence contains non-standard residues or disagrees with the recorded length.\n\
         3. Label each residue 1 if it lies inside any region with term `{IDPO_DISORDER}` (Structural state: disorder), \
            using DisProt's 1-based inclusive coordinates, and 0 otherwise. Molten-globule and pre-molten-globule regions are \
            not labelled. Function, transition and binding annotations are not used.\n\
         4. Read the CAID3 Disorder-PDB reference ({} proteins) from the CAID repository at a pinned commit.\n\
         5. Cluster DisProt sequences together with the CAID3 reference sequences using `{}`.\n\
         6. Every cluster containing a CAID3 protein goes to test, so no train or validation protein shares a 30%-identity \
            cluster with a CAID3 target. Non-CAID proteins in those clusters are also in test; filter on `caid3_target` to \
            evaluate on the exact CAID3 subset.\n\
         7. The remaining clusters are split 90/10 (by residue count) into train/validation with seed {}.\n\
         8. Test labels are DisProt's, not CAID's; CAID scores against a different negative definition (residues observed \
            in PDB), so numbers on this test split are not directly comparable to published CAID results.",
        caid_ids.len(),
        params.describe(),
        ctx.seed
    );
    let card = CardInput {
        pretty_name: "OpenDNA Disorder (DisProt 2025_12 per-residue labels)",
        license: CardLicense::Known("cc-by-4.0"),
        license_note: "DisProt is distributed under CC BY 4.0, and this derived dataset is released under CC BY 4.0. \
                       The CAID3 reference is used only to choose which proteins go to test; none of its labels are \
                       redistributed here.",
        tags: &["biology", "protein", "intrinsic-disorder", "token-classification"],
        summary: format!(
            "Per-residue intrinsic disorder labels for {} proteins from DisProt release 2025_12, prepared for training a \
             disorder head on frozen protein language model embeddings. The test split holds out the CAID3 Disorder-PDB \
             proteins and every protein in the same 30% identity cluster.",
            proteins.len()
        ),
        sources: SOURCES,
        excluded: vec![],
        construction,
        schema: &tables[0].1,
        row_unit: "proteins",
        group_unit: "clusters",
        splits: &summary,
        cohorts: &[],
        statistics,
        citations: DISORDER_BIBTEX,
        limitations: &[
            "A label of 0 means \"not annotated as disordered in DisProt\", not \"experimentally ordered\". DisProt annotates \
             regions with experimental evidence of disorder; the rest of each sequence is mostly unexamined, so negatives \
             are noisy and the positive rate is underestimated.",
            "Sequences are full UniProt sequences, some longer than 30,000 residues. ESM-2 accepts at most 1,024 residues \
             per pass; training code must crop or tile.",
            "Coverage is biased toward well-studied proteins and organisms (human, model organisms, viruses).",
            "Clustering uses 30% identity over 80% bidirectional coverage; a short disordered region shared between \
             otherwise unrelated proteins will not be detected by this criterion.",
            "Only CAID3 Disorder-PDB is held out. Older CAID rounds are not; their targets may be in train.",
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

const DISORDER_BIBTEX: &str = r#"@article{aspromonte2024disprot,
  title   = {DisProt in 2024: improving function annotation of intrinsically disordered proteins},
  author  = {Aspromonte, Maria Cristina and others},
  journal = {Nucleic Acids Research},
  volume  = {52},
  number  = {D1},
  pages   = {D434--D441},
  year    = {2024},
  doi     = {10.1093/nar/gkad928}
}

@misc{caid_repository,
  title        = {CAID: Critical Assessment of protein Intrinsic Disorder, assessment software and reference data},
  author       = {{BioComputingUP}},
  howpublished = {\url{https://github.com/BioComputingUP/CAID}},
  note         = {commit 7bab7e8880d7f949e480fdb377c5a5c8546ae336}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    const DISPROT: &str = include_str!("../../tests/fixtures/disprot_mini.json");
    const CAID: &str = include_str!("../../tests/fixtures/caid_mini.fasta");

    fn region(term: &str, start: i64, end: i64) -> DisprotRegion {
        DisprotRegion { term_id: Some(term.into()), start, end }
    }

    #[test]
    fn labels_are_one_based_inclusive_and_skip_bad_regions() {
        let mut st = LabelStats::default();
        let regions = [
            region(IDPO_DISORDER, 2, 3),
            region(IDPO_DISORDER, 3, 4), // overlapping is fine
            region(IDPO_DISORDER, 5, 9), // past the end: skipped
            region(IDPO_DISORDER, 0, 1), // start below 1: skipped
            region("IDPO:0000004", 1, 1),
            region("GO:0005515", 1, 6), // function annotation: ignored
        ];
        let l = disorder_labels(6, &regions, &mut st);
        assert_eq!(l, vec![0, 1, 1, 1, 0, 0]);
        assert_eq!(st.disorder_regions, 2);
        assert_eq!(st.out_of_range_regions, 2);
        assert_eq!(st.molten_globule_regions, 1);
    }

    #[test]
    fn reads_fixture_and_cleans() {
        let entries = read_disprot(DISPROT.as_bytes()).unwrap();
        assert_eq!(entries.len(), 4);
        let caid = read_caid_reference(CAID.as_bytes()).unwrap();
        assert_eq!(caid.len(), 2);
        let ids: BTreeSet<String> = caid.iter().map(|c| c.id.clone()).collect();
        let (ps, st) = clean(entries, &ids);
        assert_eq!(st.nonstandard, 1);
        assert_eq!(st.length_mismatch, 1);
        assert_eq!(ps.len(), 2);
        assert_eq!(ps[0].entry.disprot_id, "DP00001");
        assert_eq!(ps[0].labels, vec![1, 1, 1, 0, 0, 0, 0, 0]);
        assert!(ps[1].caid3);
    }

    #[test]
    fn caid_reference_shape_is_validated() {
        assert!(read_caid_reference(">a\nMK\n01\n>b\nM\n".as_bytes()).is_err());
        assert!(read_caid_reference(">a\nMK\n0\n".as_bytes()).is_err());
        assert!(read_caid_reference(">a\nMK\n0x\n".as_bytes()).is_err());
    }

    #[test]
    fn duplicate_ids_are_rejected() {
        let dup = r#"{"data":[{"disprot_id":"DP1","acc":"P1","sequence":"M","length":1},
                              {"disprot_id":"DP1","acc":"P2","sequence":"K","length":1}]}"#;
        assert!(read_disprot(dup.as_bytes()).is_err());
    }

    #[test]
    fn caid_members_pull_their_clusters_into_test() {
        let mut m = BTreeMap::new();
        m.insert("DP1".to_string(), "DP1".to_string());
        m.insert("DP2".to_string(), "DP1".to_string());
        m.insert("DP3".to_string(), "DP3".to_string());
        m.insert("caid_DP9".to_string(), "DP3".to_string());
        m.insert("DP4".to_string(), "DP4".to_string());
        let ids: BTreeSet<String> = ["DP2".to_string()].into();
        let t = caid_clusters(&m, &ids);
        assert_eq!(t, ["DP1".to_string(), "DP3".to_string()].into());
    }
}
