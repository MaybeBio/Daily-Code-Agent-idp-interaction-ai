//! `solubility`: eSOL cell-free solubility of E. coli proteins, with sequences
//! from UniProt.
//!
//! License review (2026-11-14):
//! * eSOL, as archived by NBDC: CC BY-SA 2.1 JP (stated in the archive README).
//!   Redistribution allowed with attribution; derivatives must use the same license.
//! * UniProtKB: CC BY 4.0.
//! * NESG / PSI:Biology: excluded, see [`crate::sources::EXCLUDED_SOURCES`].
//!
//! Because of the share-alike term, the derived dataset is CC BY-SA 2.1 JP.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufReader, Read};

use anyhow::{Context, Result};

use crate::card::{CardInput, CardLicense};
use crate::cluster::ClusterParams;
use crate::context::{write_dataset, BuildContext, FinishedDataset};
use crate::parquet_out::{ColumnData, Table};
use crate::sources::{self, Source};
use crate::split::{assign_groups, Fractions, Split, SplitSummary};
use crate::{is_standard_protein, sequence_id};

pub const NAME: &str = "solubility";
pub const ESOL_ENTRY: &str = "esol.csv";

/// Replicate solubility values further apart than this (percentage points) are dropped.
pub const MAX_DUPLICATE_SPREAD_PCT: f64 = 10.0;

#[derive(Debug, Clone, PartialEq)]
pub struct EsolRow {
    pub jw_id: String,
    pub b_number: String,
    pub gene_name: String,
    pub solubility_pct: f64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EsolStats {
    pub rows: u64,
    pub no_solubility: u64,
    pub bad_b_number: u64,
}

fn is_b_number(s: &str) -> bool {
    s.len() == 5 && s.starts_with('b') && s[1..].bytes().all(|c| c.is_ascii_digit())
}

fn header_index(h: &csv::StringRecord, name: &str) -> Result<usize> {
    h.iter().position(|x| x.trim() == name).with_context(|| format!("missing column {name:?}"))
}

/// Reads `esol.csv`. Rows without a solubility value (proteins that were not
/// quantified) are skipped and counted.
pub fn read_esol(r: impl Read) -> Result<(Vec<EsolRow>, EsolStats)> {
    let mut rdr = csv::ReaderBuilder::new().flexible(false).from_reader(r);
    let h = rdr.headers()?.clone();
    let (cj, cb, cg, cs) = (
        header_index(&h, "JW_ID")?,
        header_index(&h, "B number")?,
        header_index(&h, "Gene name K-12")?,
        header_index(&h, "Solubility (%)")?,
    );
    let mut st = EsolStats::default();
    let mut out = Vec::new();
    for rec in rdr.byte_records() {
        let rec = rec?;
        st.rows += 1;
        let field = |i: usize| String::from_utf8_lossy(&rec[i]).trim().to_owned();
        let sol = field(cs);
        let Ok(v) = sol.parse::<f64>() else {
            st.no_solubility += 1;
            continue;
        };
        let b = field(cb);
        if !is_b_number(&b) {
            st.bad_b_number += 1;
            continue;
        }
        out.push(EsolRow { jw_id: field(cj), b_number: b, gene_name: field(cg), solubility_pct: v });
    }
    Ok((out, st))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UniprotEntry {
    pub accession: String,
    pub reviewed: bool,
    pub sequence: String,
}

/// Reads the UniProt TSV and indexes entries by ordered locus name (b-number).
pub fn read_uniprot_tsv(r: impl Read) -> Result<BTreeMap<String, Vec<UniprotEntry>>> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(r);
    let h = rdr.headers()?.clone();
    let (ca, cr, co, cs) = (
        header_index(&h, "Entry")?,
        header_index(&h, "Reviewed")?,
        header_index(&h, "Gene Names (ordered locus)")?,
        header_index(&h, "Sequence")?,
    );
    let mut by_b: BTreeMap<String, Vec<UniprotEntry>> = BTreeMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let e = UniprotEntry {
            accession: rec[ca].to_owned(),
            reviewed: &rec[cr] == "reviewed",
            sequence: rec[cs].to_owned(),
        };
        for tok in rec[co].split(|c: char| c.is_whitespace() || c == ';' || c == ',') {
            if is_b_number(tok) {
                by_b.entry(tok.to_owned()).or_default().push(e.clone());
            }
        }
    }
    Ok(by_b)
}

/// Picks the UniProt entry for a b-number: the single candidate, or the single
/// reviewed candidate when there are several. Anything else is ambiguous.
pub fn resolve(cands: &[UniprotEntry]) -> Option<&UniprotEntry> {
    match cands {
        [one] => Some(one),
        _ => {
            let reviewed: Vec<_> = cands.iter().filter(|e| e.reviewed).collect();
            (reviewed.len() == 1).then(|| reviewed[0])
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub b_numbers: Vec<String>,
    pub jw_ids: Vec<String>,
    pub gene_names: Vec<String>,
    pub accession: String,
    pub sequence: String,
    pub solubility_pct: f64,
    pub n_measurements: u32,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct JoinStats {
    pub unmapped: u64,
    pub ambiguous: u64,
    pub nonstandard: u64,
    pub duplicate_groups: u64,
    pub conflicting_groups: u64,
    pub conflicting_rows: u64,
    pub above_100: u64,
}

/// Joins eSOL rows to UniProt sequences and merges rows that share a sequence.
pub fn join(
    rows: &[EsolRow],
    uniprot: &BTreeMap<String, Vec<UniprotEntry>>,
    max_spread: f64,
) -> (Vec<Record>, JoinStats) {
    let mut st = JoinStats::default();
    let mut by_seq: BTreeMap<String, Vec<(&EsolRow, &UniprotEntry)>> = BTreeMap::new();
    for r in rows {
        let Some(cands) = uniprot.get(&r.b_number) else {
            st.unmapped += 1;
            continue;
        };
        let Some(e) = resolve(cands) else {
            st.ambiguous += 1;
            continue;
        };
        if !is_standard_protein(&e.sequence) {
            st.nonstandard += 1;
            continue;
        }
        by_seq.entry(e.sequence.clone()).or_default().push((r, e));
    }
    let mut out = Vec::new();
    for (seq, group) in by_seq {
        let vals: Vec<f64> = group.iter().map(|(r, _)| r.solubility_pct).collect();
        let (min, max) =
            vals.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(*v), b.max(*v)));
        if group.len() > 1 {
            st.duplicate_groups += 1;
            if max - min > max_spread {
                st.conflicting_groups += 1;
                st.conflicting_rows += group.len() as u64;
                continue;
            }
        }
        let mean = vals.iter().sum::<f64>() / vals.len() as f64;
        if mean > 100.0 {
            st.above_100 += 1;
        }
        let uniq = |f: &dyn Fn(&EsolRow) -> &str| -> Vec<String> {
            group.iter().map(|(r, _)| f(r).to_owned()).collect::<BTreeSet<_>>().into_iter().collect()
        };
        let accession =
            group.iter().map(|(_, e)| e.accession.as_str()).min().expect("group is non-empty").to_owned();
        out.push(Record {
            b_numbers: uniq(&|r| &r.b_number),
            jw_ids: uniq(&|r| &r.jw_id),
            gene_names: uniq(&|r| &r.gene_name),
            accession,
            sequence: seq,
            solubility_pct: mean,
            n_measurements: group.len() as u32,
        });
    }
    out.sort_by(|a, b| a.b_numbers.cmp(&b.b_numbers));
    (out, st)
}

pub fn to_table(rows: &[&Record], cluster_of: &BTreeMap<String, String>) -> Table {
    let mut t = Table::new();
    t.push(
        "b_number",
        "Blattner locus tag(s) in E. coli K-12 MG1655; ';'-joined when several genes encode the same sequence",
        ColumnData::Utf8(rows.iter().map(|r| r.b_numbers.join(";")).collect()),
    );
    t.push(
        "jw_id",
        "ASKA/GenoBase JW id(s) of the expressed clone",
        ColumnData::Utf8(rows.iter().map(|r| r.jw_ids.join(";")).collect()),
    );
    t.push(
        "gene_name",
        "Gene name(s) as given by eSOL",
        ColumnData::Utf8(rows.iter().map(|r| r.gene_names.join(";")).collect()),
    );
    t.push(
        "uniprot_acc",
        "UniProtKB accession the sequence was taken from",
        ColumnData::Utf8(rows.iter().map(|r| r.accession.clone()).collect()),
    );
    t.push(
        "sequence",
        "UniProtKB canonical sequence (not the expressed construct; see limitations)",
        ColumnData::Utf8(rows.iter().map(|r| r.sequence.clone()).collect()),
    );
    t.push(
        "solubility",
        "eSOL solubility in percent: supernatant / total protein after cell-free (PURE) synthesis without chaperones. Not clipped; about 7% of upstream values exceed 100 (maximum 147)",
        ColumnData::Float64(rows.iter().map(|r| r.solubility_pct).collect()),
    );
    t.push(
        "n_measurements",
        "Number of eSOL rows averaged into solubility",
        ColumnData::Int32(rows.iter().map(|r| r.n_measurements as i32).collect()),
    );
    t.push(
        "cluster_id",
        "MMseqs2 30% identity cluster (representative sequence hash); the split unit",
        ColumnData::Utf8(rows.iter().map(|r| cluster_of[&sequence_id(&r.sequence)].clone()).collect()),
    );
    t
}

pub const SOURCES: &[Source] = &[sources::ESOL_NBDC, sources::UNIPROT_ECOLI_K12];

pub fn build(ctx: &mut BuildContext) -> Result<()> {
    let zip_path = ctx.cache.fetch(&sources::ESOL_NBDC)?;
    let up_path = ctx.cache.fetch(&sources::UNIPROT_ECOLI_K12)?;
    let mut archive = zip::ZipArchive::new(BufReader::new(File::open(&zip_path)?))?;
    let (rows, es) = read_esol(archive.by_name(ESOL_ENTRY)?)?;
    let uniprot = read_uniprot_tsv(BufReader::new(File::open(&up_path)?))?;
    let (records, js) = join(&rows, &uniprot, MAX_DUPLICATE_SPREAD_PCT);

    let seqs: BTreeMap<String, String> =
        records.iter().map(|r| (sequence_id(&r.sequence), r.sequence.clone())).collect();
    let params = ClusterParams::ID30;
    let clustering = ctx.cluster(&seqs, params)?;
    let cluster_of = clustering.member_to_rep.clone();
    let mut weights: BTreeMap<String, u64> = BTreeMap::new();
    for r in &records {
        *weights.entry(cluster_of[&sequence_id(&r.sequence)].clone()).or_default() += 1;
    }
    let fractions = Fractions::new(0.8, 0.1, 0.1)?;
    let assignment = assign_groups(&weights, &BTreeMap::new(), fractions, ctx.seed)?;
    let split_of = |r: &Record| assignment[&cluster_of[&sequence_id(&r.sequence)]];
    let summary = SplitSummary::from_rows(
        records.iter().map(|r| cluster_of[&sequence_id(&r.sequence)].as_str()),
        &assignment,
    );
    let tables: Vec<(Split, Table)> = Split::ALL
        .iter()
        .map(|s| {
            let rows: Vec<&Record> = records.iter().filter(|r| split_of(r) == *s).collect();
            (*s, to_table(&rows, &cluster_of))
        })
        .collect();
    let statistics = vec![
        ("eSOL rows".into(), es.rows.to_string()),
        ("dropped: no solubility value (not quantified)".into(), es.no_solubility.to_string()),
        ("dropped: malformed b-number".into(), es.bad_b_number.to_string()),
        ("dropped: b-number not in UniProt UP000000625".into(), js.unmapped.to_string()),
        ("dropped: b-number maps to several UniProt entries".into(), js.ambiguous.to_string()),
        ("dropped: non-standard residue".into(), js.nonstandard.to_string()),
        ("sequences with more than one eSOL row".into(), js.duplicate_groups.to_string()),
        (
            format!("dropped: conflicting duplicates (spread > {MAX_DUPLICATE_SPREAD_PCT} points)"),
            format!("{} groups / {} rows", js.conflicting_groups, js.conflicting_rows),
        ),
        ("final rows with solubility > 100%".into(), js.above_100.to_string()),
        ("final rows".into(), records.len().to_string()),
        ("MMseqs2 clusters".into(), clustering.n_clusters().to_string()),
    ];
    let construction = format!(
        "1. Read `esol.csv` from the NBDC LSDB Archive copy of eSOL and keep rows with a numeric `Solubility (%)`.\n\
         2. Map each b-number to UniProtKB through the ordered-locus names of reference proteome UP000000625. A b-number \
            with several candidate entries is resolved to the single reviewed entry, or dropped if that is not unique.\n\
         3. Merge rows that map to the same sequence: average when their values are within {MAX_DUPLICATE_SPREAD_PCT} \
            percentage points, otherwise drop them all.\n\
         4. Cluster sequences with `{}` and assign whole clusters to train/validation/test (80/10/10 of rows, seed {}).",
        params.describe(),
        ctx.seed
    );
    let card = CardInput {
        pretty_name: "OpenDNA Solubility (eSOL, E. coli cell-free)",
        license: CardLicense::Other {
            name: "cc-by-sa-2.1-jp",
            link: "https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en",
        },
        license_note: "eSOL is licensed CC BY-SA 2.1 JP by its maintainers (NBDC archive README). Its share-alike term \
                       requires derivatives to carry the same license, so this dataset is released under CC BY-SA 2.1 JP. \
                       Required attribution: \"eSOL (c) Hideki Taguchi (Tokyo Institute of Technology) licensed under CC \
                       Attribution-Share Alike 2.1 Japan\". Sequences are from UniProtKB (CC BY 4.0).",
        tags: &["biology", "protein", "protein-solubility", "escherichia-coli"],
        summary: format!(
            "Solubility of {} E. coli K-12 proteins expressed in a chaperone-free reconstituted cell-free system (eSOL, \
             Niwa et al. 2009), with sequences from UniProtKB. Prepared for training a solubility regression head on \
             frozen protein language model embeddings. Splits are by 30% sequence identity cluster.",
            records.len()
        ),
        sources: SOURCES,
        excluded: sources::EXCLUDED_SOURCES
            .iter()
            .filter(|e| e.dataset == NAME)
            .collect(),
        construction,
        schema: &tables[0].1,
        row_unit: "proteins",
        group_unit: "clusters",
        splits: &summary,
        cohorts: &[],
        statistics,
        citations: SOLUBILITY_BIBTEX,
        limitations: &[
            "One organism (E. coli K-12) and one expression system (PURE, no chaperones). Solubility in vivo, with tags, or \
             in other hosts can differ substantially.",
            "Sequences are UniProt canonical sequences. The expressed constructs came from the ASKA clone library and may \
             differ (tags, start codon choice).",
            "NESG and PSI:Biology data were requested but are excluded because no public, licensed release exists; the \
             dataset is therefore small.",
            "About 7% of eSOL values exceed 100% (up to 147%), reflecting measurement noise in the supernatant/total ratio. \n             They are left unclipped; many published solubility models clip to 100%, so do that explicitly if comparing.",
            "The share-alike license applies to any redistributed derivative of this dataset.",
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

const SOLUBILITY_BIBTEX: &str = r#"@article{niwa2009bimodal,
  title   = {Bimodal protein solubility distribution revealed by an aggregation analysis of the entire ensemble of Escherichia coli proteins},
  author  = {Niwa, Tatsuya and Ying, Bei-Wen and Saito, Katsuyo and Jin, WenZhen and Takada, Shoji and Ueda, Takuya and Taguchi, Hideki},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {106},
  number  = {11},
  pages   = {4201--4206},
  year    = {2009},
  doi     = {10.1073/pnas.0811922106}
}

@article{niwa2012chaperone,
  title   = {Global analysis of chaperone effects using a reconstituted cell-free translation system},
  author  = {Niwa, Tatsuya and Kanamori, Takashi and Ueda, Takuya and Taguchi, Hideki},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {109},
  number  = {23},
  pages   = {8937--8942},
  year    = {2012},
  doi     = {10.1073/pnas.1201380109}
}

@misc{esol_nbdc,
  title  = {eSOL, LSDB Archive},
  author = {Taguchi, Hideki},
  doi    = {10.18908/lsdba.nbdc00440-000}
}

@article{uniprot2025,
  title   = {UniProt: the Universal Protein Knowledgebase in 2025},
  author  = {{The UniProt Consortium}},
  journal = {Nucleic Acids Research},
  volume  = {53},
  number  = {D1},
  pages   = {D609--D617},
  year    = {2025},
  doi     = {10.1093/nar/gkae1010}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    const ESOL: &str = include_str!("../../tests/fixtures/esol_mini.csv");
    const UNIPROT: &str = include_str!("../../tests/fixtures/uniprot_mini.tsv");

    #[test]
    fn reads_esol_and_skips_unquantified() {
        let (rows, st) = read_esol(ESOL.as_bytes()).unwrap();
        assert_eq!(st.rows, 6);
        assert_eq!(st.no_solubility, 1);
        assert_eq!(rows.len(), 5);
        assert_eq!(rows[0].b_number, "b0003");
        assert_eq!(rows[0].solubility_pct, 32.0);
    }

    #[test]
    fn indexes_uniprot_by_locus_and_resolves_ambiguity() {
        let up = read_uniprot_tsv(UNIPROT.as_bytes()).unwrap();
        assert_eq!(up["b0003"].len(), 1);
        // b0005 has a reviewed and an unreviewed entry: the reviewed one wins.
        assert_eq!(resolve(&up["b0005"]).unwrap().accession, "P0AAA5");
        // b0007 has two unreviewed entries: ambiguous.
        assert!(resolve(&up["b0007"]).is_none());
    }

    #[test]
    fn joins_and_merges_by_sequence() {
        let (rows, _) = read_esol(ESOL.as_bytes()).unwrap();
        let up = read_uniprot_tsv(UNIPROT.as_bytes()).unwrap();
        let (recs, st) = join(&rows, &up, 10.0);
        assert_eq!(st.unmapped, 1); // b0009
        assert_eq!(st.ambiguous, 1); // b0007
        assert_eq!(st.duplicate_groups, 1);
        assert_eq!(st.conflicting_groups, 0);
        assert_eq!(recs.len(), 2);
        // b0003 and b0004 encode the same sequence in the fixture and are merged.
        let merged = &recs[0];
        assert_eq!(merged.b_numbers, ["b0003", "b0004"]);
        assert_eq!(merged.n_measurements, 2);
        assert!((merged.solubility_pct - 35.0).abs() < 1e-12);
        // With a tighter tolerance the pair conflicts and is dropped.
        let (recs, st) = join(&rows, &up, 5.0);
        assert_eq!(st.conflicting_groups, 1);
        assert_eq!(recs.len(), 1);
    }
}
