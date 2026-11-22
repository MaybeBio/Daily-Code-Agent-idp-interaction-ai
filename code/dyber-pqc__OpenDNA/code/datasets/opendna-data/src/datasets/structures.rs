//! `structures`: CATH domains (S100 representatives) with their CATH labels,
//! joined with wwPDB entry metadata, split by CATH topology.
//!
//! Leakage control: two domains end up in the same split if they share a CATH
//! topology (C.A.T) **or** an MMseqs2 30% identity cluster. Topologies and
//! clusters are nodes of a graph, each domain is an edge between its topology
//! and its cluster, and connected components are the split unit. A component
//! can be large when sequence clusters bridge topologies; its size is reported.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{bail, Context, Result};

use crate::card::{CardInput, CardLicense};
use crate::cluster::ClusterParams;
use crate::context::{write_dataset, BuildContext, FinishedDataset};
use crate::parquet_out::{ColumnData, Table};
use crate::sources::{self, Source};
use crate::split::{assign_groups, Fractions, Split, SplitSummary};
use crate::unionfind::UnionFind;
use crate::{fasta, is_standard_protein};

pub const NAME: &str = "structures";

#[derive(Debug, Clone, PartialEq)]
pub struct CathDomain {
    pub domain_id: String,
    pub class: u32,
    pub architecture: u32,
    pub topology: u32,
    pub homology: u32,
    pub length: u32,
    /// Angstrom; `None` for CATH's sentinel values (999 = NMR, 1000 = obsolete/unknown).
    pub resolution: Option<f64>,
}

impl CathDomain {
    pub fn topology_code(&self) -> String {
        format!("{}.{}.{}", self.class, self.architecture, self.topology)
    }

    pub fn superfamily_code(&self) -> String {
        format!("{}.{}", self.topology_code(), self.homology)
    }

    pub fn pdb_id(&self) -> &str {
        &self.domain_id[..4]
    }

    pub fn chain_id(&self) -> &str {
        &self.domain_id[4..5]
    }
}

/// Parses a CATH domain list (CLF 2.0): 12 whitespace-separated columns.
pub fn read_domain_list(r: impl BufRead) -> Result<Vec<CathDomain>> {
    let mut out = Vec::new();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() != 12 || f[0].len() != 7 {
            bail!("CATH domain list line {} is malformed: {line:?}", i + 1);
        }
        let num = |k: usize| -> Result<u32> {
            f[k].parse().with_context(|| format!("line {}: column {} is not an integer", i + 1, k + 1))
        };
        let res: f64 = f[11].parse().with_context(|| format!("line {}: bad resolution", i + 1))?;
        out.push(CathDomain {
            domain_id: f[0].to_owned(),
            class: num(1)?,
            architecture: num(2)?,
            topology: num(3)?,
            homology: num(4)?,
            length: num(10)?,
            resolution: (res < 999.0).then_some(res),
        });
    }
    Ok(out)
}

/// Domain id and segment string from a CATH FASTA header such as
/// `cath|4_4_0|107lA00/1-162`.
pub fn parse_cath_header(h: &str) -> Option<(&str, &str)> {
    let last = h.split('|').nth(2)?;
    let (id, segs) = last.split_once('/')?;
    (id.len() == 7).then_some((id, segs))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PdbEntry {
    pub deposition_date: Option<String>,
    pub method: String,
}

/// Converts the wwPDB index date `MM/DD/YY` to ISO 8601. Two-digit years
/// below 70 are 20xx (the PDB was founded in 1971).
pub fn iso_date(mmddyy: &str) -> Option<String> {
    let p: Vec<&str> = mmddyy.split('/').collect();
    let [m, d, y] = p.as_slice() else {
        return None;
    };
    let (m, d, y): (u32, u32, u32) = (m.parse().ok()?, d.parse().ok()?, y.parse().ok()?);
    if !(1..=12).contains(&m) || !(1..=31).contains(&d) || y > 99 {
        return None;
    }
    let year = if y < 70 { 2000 + y } else { 1900 + y };
    Some(format!("{year:04}-{m:02}-{d:02}"))
}

/// Parses `entries.idx` (two header lines, then tab-separated rows).
pub fn read_entries_idx(r: impl BufRead) -> Result<HashMap<String, PdbEntry>> {
    let mut out = HashMap::new();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        if i < 2 || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            bail!("entries.idx line {} has {} fields", i + 1, f.len());
        }
        out.insert(
            f[0].trim().to_ascii_lowercase(),
            PdbEntry { deposition_date: iso_date(f[2].trim()), method: f[7].trim().to_owned() },
        );
    }
    Ok(out)
}

#[derive(Debug, Clone)]
pub struct Domain {
    pub cath: CathDomain,
    pub segments: String,
    pub sequence: String,
    pub pdb: Option<PdbEntry>,
}

/// Connected components over topology and sequence-cluster nodes.
/// Returns domain id -> component name.
pub fn topology_components(
    domains: &[(String, String, String)], // (domain id, topology code, cluster rep)
) -> BTreeMap<String, String> {
    let mut uf = UnionFind::new();
    for (_, top, clu) in domains {
        uf.union(&format!("T:{top}"), &format!("C:{clu}"));
    }
    let comp = uf.components();
    domains.iter().map(|(id, top, _)| (id.clone(), comp[&format!("T:{top}")].clone())).collect()
}

fn to_table(
    rows: &[&Domain],
    cluster_of: &BTreeMap<String, String>,
    comp: &BTreeMap<String, String>,
) -> Table {
    let mut t = Table::new();
    let id = |d: &&Domain| d.cath.domain_id.clone();
    t.push(
        "domain_id",
        "CATH domain id (PDB id, chain, domain number)",
        ColumnData::Utf8(rows.iter().map(id).collect()),
    );
    t.push(
        "pdb_id",
        "PDB entry id (lowercase)",
        ColumnData::Utf8(rows.iter().map(|d| d.cath.pdb_id().to_owned()).collect()),
    );
    t.push(
        "chain_id",
        "PDB chain id",
        ColumnData::Utf8(rows.iter().map(|d| d.cath.chain_id().to_owned()).collect()),
    );
    t.push(
        "segments",
        "Domain boundaries in PDB residue numbering as given in the CATH FASTA header",
        ColumnData::Utf8(rows.iter().map(|d| d.segments.clone()).collect()),
    );
    t.push(
        "sequence",
        "Domain sequence from ATOM records (CATH)",
        ColumnData::Utf8(rows.iter().map(|d| d.sequence.clone()).collect()),
    );
    t.push(
        "length",
        "Sequence length",
        ColumnData::Int32(rows.iter().map(|d| d.sequence.len() as i32).collect()),
    );
    t.push(
        "cath_class",
        "CATH class (C)",
        ColumnData::Int32(rows.iter().map(|d| d.cath.class as i32).collect()),
    );
    t.push(
        "cath_architecture",
        "CATH architecture code (C.A)",
        ColumnData::Utf8(rows.iter().map(|d| format!("{}.{}", d.cath.class, d.cath.architecture)).collect()),
    );
    t.push(
        "cath_topology",
        "CATH topology code (C.A.T)",
        ColumnData::Utf8(rows.iter().map(|d| d.cath.topology_code()).collect()),
    );
    t.push(
        "cath_superfamily",
        "CATH homologous superfamily code (C.A.T.H)",
        ColumnData::Utf8(rows.iter().map(|d| d.cath.superfamily_code()).collect()),
    );
    t.push(
        "resolution",
        "Resolution in Angstrom from CATH; null for NMR and other non-diffraction methods",
        ColumnData::OptFloat64(rows.iter().map(|d| d.cath.resolution).collect()),
    );
    t.push(
        "experimental_method",
        "Experimental method from the wwPDB entries index",
        ColumnData::OptUtf8(rows.iter().map(|d| d.pdb.as_ref().map(|p| p.method.clone())).collect()),
    );
    t.push(
        "deposition_date",
        "Deposition date (ISO 8601) from the wwPDB entries index",
        ColumnData::OptUtf8(
            rows.iter().map(|d| d.pdb.as_ref().and_then(|p| p.deposition_date.clone())).collect(),
        ),
    );
    t.push(
        "cluster_id",
        "MMseqs2 30% identity cluster (representative domain id)",
        ColumnData::Utf8(rows.iter().map(|d| cluster_of[&d.cath.domain_id].clone()).collect()),
    );
    t.push(
        "split_group",
        "Connected component of topologies and clusters; the split unit",
        ColumnData::Utf8(rows.iter().map(|d| comp[&d.cath.domain_id].clone()).collect()),
    );
    t
}

pub const SOURCES: &[Source] =
    &[sources::CATH_DOMAIN_LIST, sources::CATH_DOMAIN_SEQS, sources::PDB_ENTRIES_IDX];

pub fn build(ctx: &mut BuildContext) -> Result<()> {
    let list_path = ctx.cache.fetch(&sources::CATH_DOMAIN_LIST)?;
    let seq_path = ctx.cache.fetch(&sources::CATH_DOMAIN_SEQS)?;
    let idx_path = ctx.cache.fetch(&sources::PDB_ENTRIES_IDX)?;
    let cath = read_domain_list(BufReader::new(File::open(&list_path)?))?;
    let fa = fasta::read(BufReader::new(File::open(&seq_path)?))?;
    let pdb = read_entries_idx(BufReader::new(File::open(&idx_path)?))?;

    let mut seq_by_id: HashMap<String, (String, String)> = HashMap::new();
    for r in fa {
        let Some((id, segs)) = parse_cath_header(&r.header) else {
            bail!("unexpected CATH FASTA header {:?}", r.header);
        };
        seq_by_id.insert(id.to_owned(), (segs.to_owned(), r.sequence));
    }
    let (mut no_seq, mut nonstd, mut no_pdb) = (0u64, 0u64, 0u64);
    let mut domains = Vec::new();
    for c in cath {
        let Some((segments, sequence)) = seq_by_id.remove(&c.domain_id) else {
            no_seq += 1;
            continue;
        };
        if !is_standard_protein(&sequence) {
            nonstd += 1;
            continue;
        }
        let entry = pdb.get(c.pdb_id()).cloned();
        if entry.is_none() {
            no_pdb += 1;
        }
        domains.push(Domain { cath: c, segments, sequence, pdb: entry });
    }
    domains.sort_by(|a, b| a.cath.domain_id.cmp(&b.cath.domain_id));

    let seqs: BTreeMap<String, String> =
        domains.iter().map(|d| (d.cath.domain_id.clone(), d.sequence.clone())).collect();
    let params = ClusterParams::ID30;
    let clustering = ctx.cluster(&seqs, params)?;
    let cluster_of = clustering.member_to_rep.clone();
    let triples: Vec<(String, String, String)> = domains
        .iter()
        .map(|d| (d.cath.domain_id.clone(), d.cath.topology_code(), cluster_of[&d.cath.domain_id].clone()))
        .collect();
    let comp = topology_components(&triples);
    let mut weights: BTreeMap<String, u64> = BTreeMap::new();
    for g in comp.values() {
        *weights.entry(g.clone()).or_default() += 1;
    }
    let largest = weights.values().copied().max().unwrap_or(0);
    let fractions = Fractions::new(0.8, 0.1, 0.1)?;
    let assignment = assign_groups(&weights, &BTreeMap::new(), fractions, ctx.seed)?;
    let split_of = |d: &Domain| assignment[&comp[&d.cath.domain_id]];
    let summary =
        SplitSummary::from_rows(domains.iter().map(|d| comp[&d.cath.domain_id].as_str()), &assignment);
    let tables: Vec<(Split, Table)> = Split::ALL
        .iter()
        .map(|s| {
            let rows: Vec<&Domain> = domains.iter().filter(|d| split_of(d) == *s).collect();
            (*s, to_table(&rows, &cluster_of, &comp))
        })
        .collect();
    let n_topologies = triples.iter().map(|t| t.1.as_str()).collect::<std::collections::BTreeSet<_>>().len();
    let statistics = vec![
        ("dropped: no CATH sequence".into(), no_seq.to_string()),
        ("dropped: non-standard residue".into(), nonstd.to_string()),
        ("kept without wwPDB metadata (entry obsolete or missing)".into(), no_pdb.to_string()),
        ("final domains".into(), domains.len().to_string()),
        ("CATH topologies".into(), n_topologies.to_string()),
        ("MMseqs2 clusters".into(), clustering.n_clusters().to_string()),
        ("split groups (components)".into(), weights.len().to_string()),
        ("largest split group (domains)".into(), largest.to_string()),
    ];
    let construction = format!(
        "1. Read the CATH v4.4.0 S100 domain list and the matching ATOM-record domain sequences.\n\
         2. Drop domains with no sequence or with non-standard residues. Join wwPDB entry metadata (method, deposition date) by PDB id.\n\
         3. Cluster domain sequences with `{}`.\n\
         4. Build a graph whose nodes are CATH topologies and sequence clusters, with one edge per domain; connected \
            components are the split unit. No topology and no 30% identity cluster spans two splits.\n\
         5. Assign components to train/validation/test (80/10/10 of domains, seed {}).",
        params.describe(),
        ctx.seed
    );
    let card = CardInput {
        pretty_name: "OpenDNA Structures (CATH v4.4.0 S100 domains)",
        license: CardLicense::Known("cc-by-4.0"),
        license_note: "CATH is CC BY 4.0 and wwPDB data is CC0 1.0; the derived dataset is released under CC BY 4.0.",
        tags: &["biology", "protein", "protein-structure", "cath"],
        summary: format!(
            "{} CATH v4.4.0 S100 domain sequences with CATH class, architecture, topology and superfamily labels and wwPDB \
             entry metadata. Splits are by CATH topology, merged with 30% identity sequence clusters, so no fold or close \
             homolog appears in more than one split. Coordinates are not included.",
            domains.len()
        ),
        sources: SOURCES,
        excluded: vec![],
        construction,
        schema: &tables[0].1,
        row_unit: "domains",
        group_unit: "components",
        splits: &summary,
        cohorts: &[],
        statistics,
        citations: STRUCTURES_BIBTEX,
        limitations: &[
            "Coordinates are not included; this release carries sequences and labels only.",
            "S100 representatives remove identical sequences only, so near-identical domains remain within a split.",
            "Splitting by topology means test folds are unseen in training. That is a harder setting than a random split \
             and the realised split fractions can deviate from 80/10/10 when a component is large.",
            "Domain sequences come from ATOM records, so unobserved residues are absent.",
            "Domains whose CATH sequence contains `X` (unknown or modified residues) are dropped; this removes about 13% \n             of S100 domains and may under-represent structures with many modified residues.",
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

const STRUCTURES_BIBTEX: &str = r#"@article{sillitoe2021cath,
  title   = {CATH: increased structural coverage of functional space},
  author  = {Sillitoe, Ian and Bordin, Nicola and Dawson, Natalie and others},
  journal = {Nucleic Acids Research},
  volume  = {49},
  number  = {D1},
  pages   = {D266--D273},
  year    = {2021},
  doi     = {10.1093/nar/gkaa1079}
}

@article{wwpdb2019,
  title   = {Protein Data Bank: the single global archive for 3D macromolecular structure data},
  author  = {{wwPDB consortium}},
  journal = {Nucleic Acids Research},
  volume  = {47},
  number  = {D1},
  pages   = {D520--D528},
  year    = {2019},
  doi     = {10.1093/nar/gky949}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    const LIST: &str = include_str!("../../tests/fixtures/cath_domain_list_mini.txt");
    const IDX: &str = include_str!("../../tests/fixtures/entries_mini.idx");

    #[test]
    fn parses_domain_list() {
        let d = read_domain_list(LIST.as_bytes()).unwrap();
        assert_eq!(d.len(), 3);
        assert_eq!(d[0].domain_id, "1oaiA00");
        assert_eq!(d[0].topology_code(), "1.10.8");
        assert_eq!(d[0].superfamily_code(), "1.10.8.10");
        assert_eq!(d[0].resolution, Some(1.0));
        assert_eq!(d[1].resolution, None); // 999.000 = NMR
        assert_eq!(d[2].pdb_id(), "3frh");
        assert_eq!(d[2].chain_id(), "A");
        assert!(read_domain_list("1oaiA00 1 10\n".as_bytes()).is_err());
    }

    #[test]
    fn parses_cath_fasta_header() {
        assert_eq!(parse_cath_header("cath|4_4_0|107lA00/1-162"), Some(("107lA00", "1-162")));
        assert_eq!(parse_cath_header("cath|4_4_0|1abcB02/10-50,60-80"), Some(("1abcB02", "10-50,60-80")));
        assert_eq!(parse_cath_header("nonsense"), None);
    }

    #[test]
    fn parses_entries_idx_and_dates() {
        let e = read_entries_idx(IDX.as_bytes()).unwrap();
        assert_eq!(e["101m"].method, "X-RAY DIFFRACTION");
        assert_eq!(e["101m"].deposition_date.as_deref(), Some("1997-12-13"));
        assert_eq!(iso_date("01/02/05").as_deref(), Some("2005-01-02"));
        assert_eq!(iso_date("13/02/05"), None);
        assert_eq!(iso_date("garbage"), None);
    }

    #[test]
    fn components_merge_topologies_bridged_by_a_cluster() {
        let t = |a: &str, b: &str, c: &str| (a.to_string(), b.to_string(), c.to_string());
        let comp = topology_components(&[
            t("d1", "1.10.8", "c1"),
            t("d2", "1.10.8", "c2"),
            t("d3", "2.40.50", "c2"), // shares cluster c2 with d2: merges the two topologies
            t("d4", "3.30.70", "c9"),
        ]);
        assert_eq!(comp["d1"], comp["d3"]);
        assert_ne!(comp["d1"], comp["d4"]);
    }
}
