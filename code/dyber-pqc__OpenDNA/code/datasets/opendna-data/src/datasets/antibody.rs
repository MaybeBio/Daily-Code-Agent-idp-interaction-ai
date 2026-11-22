//! `antibody`: heavy-chain variable regions from human and mouse repertoires.
//!
//! License review (2026-11-17):
//! * Observed Antibody Space (OAS): CC BY 4.0, stated by the maintainers on
//!   <https://opig.stats.ox.ac.uk/webapps/oas/> -- "the data contained within
//!   OAS is available under a CC-BY 4.0 license". Redistribution is permitted
//!   with attribution, so the derived dataset is CC BY 4.0.
//! * SAbDab: not used. It holds structures, and this is a sequence dataset.
//!   Its terms were not relied on and so were not reviewed.
//!
//! # What the label is, and what it is not
//!
//! Each row carries `is_human`: whether the sequence came out of a human
//! repertoire rather than a mouse one. That is a fact recorded by the study
//! that deposited it, not a prediction and not a measurement of anything
//! about the molecule.
//!
//! It is emphatically not an immunogenicity label. A sequence that resembles
//! a human repertoire is not thereby safe in humans, and the reverse is not
//! true either. What a head trained on this can answer is "does this look
//! like something a human B cell made", which is the question asked of a
//! humanised candidate -- and no more than that.
//!
//! # Why there is no train/test split
//!
//! This dataset was built to be split by group, so that neither a donor nor a
//! clonal lineage could span train and test. It cannot be. The mouse cohort is
//! one indivisible group, and two measurements on the real data say so:
//!
//! * Clustering the amino-acid V(D)J region at 90% identity put **all 37
//!   labelled donors in a single component**: one group over 176,543 rows. The
//!   region is mostly germline-constant, so 90% identity over its whole length
//!   groups V-gene families rather than lineages, and joining donors through
//!   shared clusters then merges everything transitively.
//! * Under the strictest reading of a shared lineage -- an identical
//!   rearrangement, meaning the same V gene, the same J gene and a byte-identical
//!   CDR3 -- the human donors mostly stay apart (13 alone, two groups of three),
//!   but **all 18 mouse donors join into one group**, across three inbred
//!   strains and nine studies.
//!
//! The second result is biology, not a loose threshold. C57BL/6 and BALB/c
//! animals are genetically near-identical, and their repertoires contain the
//! same public rearrangements; one laboratory mouse is not an independent
//! sample of a species the way one person is. So there is no grouped split
//! that puts mouse rows on both sides, and a split that leaves the entire
//! mouse cohort in `train` is not a split.
//!
//! What this dataset is, therefore, is four **cohorts** for evaluation, with
//! the grouping recorded in a `donor_group` column so anyone who wants to
//! split the human side can, knowing what the mouse side does. It is not a
//! training set for a species classifier -- which, separately, a lookup table
//! on the V gene call already answers at 98.5% held-out accuracy against a
//! 68.9% majority baseline, so there was never a model worth training here.
//!
//! # Three traps in the upstream metadata
//!
//! **The species field names the animal, not the immunoglobulin.** OAS
//! contains `HIS-Mouse` (a mouse carrying an engrafted human immune system)
//! and `Kymouse` (a transgenic mouse carrying human immunoglobulin loci).
//! Both are mice whose antibody V genes are human. Labelling them non-human
//! would teach a head that human V genes are mouse, on millions of sequences.
//! They are excluded from the labelled cohorts and kept as a separate
//! diagnostic set, where calling them human is the correct behaviour: that is
//! the humanisation question rather than a species-identification one.
//!
//! **One donor can appear in several studies.** Fifty individual identifiers
//! in OAS occur under more than one study. Some are the same person -- CAP256
//! is one HIV donor in Bhiman 2015, Doria-Rose 2015 and Johnson 2018 -- and
//! some are study-local placeholders, where `Donor-1` in two papers is two
//! people. Keying a donor by name alone merges those placeholders, which
//! costs a little data; keying by study would put one real donor in two
//! cohorts or two folds, which inflates every number computed from them.
//! [`donor_key`] takes the first trade, and says so.
//!
//! **A unit is two gzip members, not one.** OAS writes the metadata line as
//! its own gzip member and the CSV body as a second one. A single-member
//! decoder -- `flate2::read::GzDecoder`, the obvious choice -- stops after the
//! first, so it yields the metadata and no sequences at all, with no error.
//! Everything here reads units through [`open_unit`], which uses
//! `MultiGzDecoder`, and [`read_unit`] refuses a unit whose body is empty
//! rather than returning nothing and letting a build quietly produce a
//! dataset of zero rows.
//!
//! **Only mouse has enough independent studies to hold one out.** Human has
//! 69 studies and mouse 15, but rat, rabbit and bactrian camel have exactly
//! one study each, and rhesus two. A test score for those species would be
//! inseparable from one laboratory's protocol, so they are not part of the
//! labelled cohorts; they are kept as an out-of-distribution set whose
//! single-study origin the card states.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, BufReader, Read};

use anyhow::{bail, Context, Result};
use serde::Deserialize;

use crate::card::{CardInput, CardLicense};
use crate::context::{write_dataset, BuildContext, FinishedDataset};
use crate::is_standard_protein;
use crate::parquet_out::{ColumnData, Table};
use crate::sources::{self, Source};
use crate::split::SplitSummary;
use crate::unionfind::UnionFind;

pub const NAME: &str = "antibody";

/// Shortest and longest heavy-chain variable region kept, in residues.
///
/// A rearranged IGH V(D)J region runs about 110 to 130 residues. The band is
/// wide enough to keep unusual but real sequences and narrow enough to drop
/// the truncated reads that a sequencing run leaves behind; what falls outside
/// it is counted rather than discarded quietly.
pub const MIN_RESIDUES: usize = 90;
pub const MAX_RESIDUES: usize = 200;

/// What a unit's own metadata header says about it.
///
/// OAS puts this as a JSON object in the first line of every `csv.gz`, so the
/// organism and donor travel with the sequences rather than being recorded
/// separately here and going stale.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
pub struct UnitMetadata {
    #[serde(rename = "Species")]
    pub species: String,
    #[serde(rename = "Subject")]
    pub subject: String,
    #[serde(rename = "Chain")]
    pub chain: String,
    #[serde(rename = "Author")]
    pub author: String,
    /// Reported, never depended on. The `loose_string` deserializer below
    /// records why these three are read leniently: OAS does not type them
    /// consistently, and nothing here depends on their shape. Named in prose
    /// rather than linked, because it is private and a link to it would only
    /// resolve for someone building the docs with private items.
    #[serde(rename = "Isotype", default, deserialize_with = "loose_string")]
    pub isotype: String,
    #[serde(rename = "Disease", default, deserialize_with = "loose_string")]
    pub disease: String,
    #[serde(rename = "BType", default, deserialize_with = "loose_string")]
    pub b_type: String,
}

/// Read a metadata field that upstream does not type consistently.
///
/// OAS writes `Isotype` as `"Bulk"` in Greiff 2017 and as `["Bulk"]` in Li
/// 2017, and `Age` as `"no"` in one and `5` in another. The three fields this
/// is used on are reported in the card and never decide anything, so the
/// honest handling is to render whatever is there and carry on -- failing a
/// whole build over the shape of a field nothing depends on would be a
/// strictness that costs data and buys nothing.
///
/// The fields that do decide something -- `Species`, `Subject`, `Chain` --
/// stay plain strings, so a change in any of those is still an error.
fn loose_string<'de, D: serde::Deserializer<'de>>(d: D) -> std::result::Result<String, D::Error> {
    use serde::Deserialize as _;
    Ok(match serde_json::Value::deserialize(d)? {
        serde_json::Value::Null => String::new(),
        serde_json::Value::String(s) => s,
        serde_json::Value::Array(items) => items
            .iter()
            .map(|i| i.as_str().map_or_else(|| i.to_string(), str::to_owned))
            .collect::<Vec<_>>()
            .join(";"),
        other => other.to_string(),
    })
}

/// Which side of the question a unit is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Cohort {
    /// Human repertoire. Its own labelled cohort.
    Human,
    /// Mouse repertoire. Its own labelled cohort.
    Mouse,
    /// A mouse whose immunoglobulin loci are human: `HIS-Mouse`, `Kymouse`.
    /// Never labelled; the humanisation diagnostic, where a humanness measure
    /// answering "human" is right.
    HumanLociInMouse,
    /// Another species. Never labelled: one study each for most of them, so a
    /// score on one cannot be separated from one laboratory's protocol.
    OutOfDistribution,
}

impl Cohort {
    /// The split file a cohort's rows are written to.
    pub fn table(self) -> &'static str {
        match self {
            Self::Human | Self::Mouse => "labelled",
            Self::HumanLociInMouse => "diagnostic",
            Self::OutOfDistribution => "out_of_distribution",
        }
    }

    /// The label, where there is one.
    ///
    /// `None` for the two held-out cohorts, rather than a default. A
    /// transgenic mouse recorded as `is_human = false` is the exact mistake
    /// this enum exists to prevent, and a column that cannot hold "no answer"
    /// invites it.
    pub fn is_human(self) -> Option<bool> {
        match self {
            Self::Human => Some(true),
            Self::Mouse => Some(false),
            Self::HumanLociInMouse | Self::OutOfDistribution => None,
        }
    }
}

/// Read a cohort off the `Species` string.
///
/// Matched exactly against the values OAS actually uses, listed here in full,
/// so a new species appearing upstream is an error rather than a silent
/// assignment to whichever arm the match arms happen to end with.
///
/// # Errors
/// If the species is not one of the values this has been checked against.
pub fn cohort_of(species: &str) -> Result<Cohort> {
    Ok(match species {
        "human" => Cohort::Human,
        // Strains, all of them ordinary mice.
        "mouse" | "mouse_C57BL/6" | "mouse_BALB/c" | "mouse_RAG2-GFP/129Sve" | "mouse_Swiss-Webster" => {
            Cohort::Mouse
        }
        // Mice with human immunoglobulin. See the module documentation.
        "HIS-Mouse" | "Kymouse" => Cohort::HumanLociInMouse,
        "rat" | "rabbit" | "rhesus" | "bactrian_camel" => Cohort::OutOfDistribution,
        other => bail!(
            "unknown OAS species {other:?}. Add it to `cohort_of` deliberately: a mouse \
             carrying human immunoglobulin loci must not be labelled non-human, and a species \
             with one study behind it must not be trained on."
        ),
    })
}

/// The identifier a split may never be drawn through.
///
/// A named subject is global, because the same donor is deposited under the
/// same name by different studies. An absent name is a missing field rather
/// than a person, so those units fall back to the study, which keeps the
/// boundary safe without inventing an identity.
pub fn donor_key(metadata: &UnitMetadata, study: &str) -> String {
    let subject = metadata.subject.trim();
    if matches!(subject, "" | "no" | "No" | "NA" | "n/a" | "none" | "None" | "unknown") {
        format!("study:{study}")
    } else {
        format!("person:{subject}")
    }
}

/// One kept sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AntibodySequence {
    /// The amino-acid V(D)J region, as OAS aligned it.
    pub sequence: String,
    /// IMGT V gene call, as given. Kept for reporting, never as a feature.
    pub v_call: String,
    /// IMGT J gene call, as given.
    pub j_call: String,
    /// CDR3, amino acids.
    pub cdr3: String,
    /// How many identical reads OAS collapsed into this row.
    pub redundancy: u32,
}

/// Why rows were dropped, so the card can account for every one.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadStats {
    pub rows: u64,
    pub not_productive: u64,
    pub stop_codon: u64,
    pub out_of_frame: u64,
    pub wrong_locus: u64,
    pub empty_sequence: u64,
    pub nonstandard_residue: u64,
    pub out_of_length_band: u64,
    pub duplicate_within_unit: u64,
    pub kept: u64,
}

fn column(header: &csv::StringRecord, name: &str) -> Result<usize> {
    header
        .iter()
        .position(|field| field.trim() == name)
        .with_context(|| format!("OAS unit has no column {name:?}"))
}

/// Read one OAS `csv.gz` data unit: its metadata header and its usable rows.
///
/// The first line is a single CSV field holding a JSON object; the AIRR
/// header is the second line. Rows are filtered on the flags OAS already
/// computed -- productive, in frame, no stop codon -- rather than on anything
/// re-derived here, and every rejection is counted.
///
/// # Errors
/// If the stream is not an OAS unit: no metadata line, unparsable JSON, or a
/// missing AIRR column.
pub fn read_unit(reader: impl Read) -> Result<(UnitMetadata, Vec<AntibodySequence>, ReadStats)> {
    let mut buffered = BufReader::new(reader);
    let mut first = String::new();
    if buffered.read_line(&mut first)? == 0 {
        bail!("empty OAS unit: expected a metadata line");
    }
    let field = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_reader(first.as_bytes())
        .records()
        .next()
        .transpose()?
        .and_then(|record| record.get(0).map(str::to_owned))
        .context("OAS metadata line is not a CSV record")?;
    let metadata: UnitMetadata =
        serde_json::from_str(&field).context("OAS metadata line is not the expected JSON")?;

    let mut rdr = csv::ReaderBuilder::new().flexible(true).from_reader(buffered);
    let header = rdr.headers()?.clone();
    let (c_seq, c_prod, c_stop, c_frame, c_locus) = (
        column(&header, "sequence_alignment_aa")?,
        column(&header, "productive")?,
        column(&header, "stop_codon")?,
        column(&header, "vj_in_frame")?,
        column(&header, "locus")?,
    );
    let (c_v, c_j, c_cdr3, c_red) = (
        column(&header, "v_call")?,
        column(&header, "j_call")?,
        column(&header, "cdr3_aa")?,
        column(&header, "Redundancy")?,
    );

    let mut stats = ReadStats::default();
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();
    for record in rdr.records() {
        let record = record?;
        stats.rows += 1;
        let get = |i: usize| record.get(i).unwrap_or("").trim();

        // OAS writes T/F. Anything else is treated as not asserted, which is
        // the conservative reading for a flag that gates inclusion.
        if get(c_prod) != "T" {
            stats.not_productive += 1;
            continue;
        }
        if get(c_stop) != "F" {
            stats.stop_codon += 1;
            continue;
        }
        if get(c_frame) != "T" {
            stats.out_of_frame += 1;
            continue;
        }
        if get(c_locus) != "H" {
            stats.wrong_locus += 1;
            continue;
        }
        let sequence = get(c_seq).to_ascii_uppercase();
        if sequence.is_empty() {
            stats.empty_sequence += 1;
            continue;
        }
        // Rejects alignment gaps and ambiguity codes together, which is what
        // we want: an embedder should not be fed a dot.
        if !is_standard_protein(&sequence) {
            stats.nonstandard_residue += 1;
            continue;
        }
        if !(MIN_RESIDUES..=MAX_RESIDUES).contains(&sequence.len()) {
            stats.out_of_length_band += 1;
            continue;
        }
        if !seen.insert(sequence.clone()) {
            stats.duplicate_within_unit += 1;
            continue;
        }
        stats.kept += 1;
        out.push(AntibodySequence {
            sequence,
            v_call: get(c_v).to_owned(),
            j_call: get(c_j).to_owned(),
            cdr3: get(c_cdr3).to_owned(),
            redundancy: get(c_red).parse().unwrap_or(1),
        });
    }
    // A unit that parsed its metadata and then found no rows at all is the
    // signature of reading only the first gzip member. It is never a real
    // unit: OAS does not publish empty ones.
    if stats.rows == 0 {
        bail!(
            "OAS unit for {} / {} has metadata but no rows. If it was opened with a              single-member gzip decoder, that is the cause: a unit is two concatenated              members and the body is the second.",
            metadata.species,
            metadata.subject
        );
    }
    Ok((metadata, out, stats))
}

/// Open a cached OAS unit for reading.
///
/// The only way units are opened, so the two-member gzip layout is handled in
/// exactly one place.
///
/// # Errors
/// If the file cannot be opened.
pub fn open_unit(path: &std::path::Path) -> Result<impl Read> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(flate2::read::MultiGzDecoder::new(std::io::BufReader::new(file)))
}

/// Add one unit's counts into a running total.
pub fn accumulate(total: &mut ReadStats, one: &ReadStats) {
    total.rows += one.rows;
    total.not_productive += one.not_productive;
    total.stop_codon += one.stop_codon;
    total.out_of_frame += one.out_of_frame;
    total.wrong_locus += one.wrong_locus;
    total.empty_sequence += one.empty_sequence;
    total.nonstandard_residue += one.nonstandard_residue;
    total.out_of_length_band += one.out_of_length_band;
    total.duplicate_within_unit += one.duplicate_within_unit;
    total.kept += one.kept;
}

/// Every V gene seen in a cohort, with how often, for the card.
///
/// Reported because it is the first thing to check when a species classifier
/// looks too good: if one cohort's sequences are separable by their V call
/// alone, the head has learnt the germline catalogue rather than anything
/// about the molecule, and a reader deserves to see the overlap.
pub fn v_gene_counts<'a>(v_calls: impl IntoIterator<Item = &'a str>) -> BTreeMap<String, u64> {
    let mut counts = BTreeMap::new();
    for call in v_calls {
        // IGHV1-52*01 -> IGHV1-52: the allele is finer than anything here
        // needs and splits the counts across near-identical genes.
        let gene = call.split('*').next().unwrap_or(call);
        if !gene.is_empty() {
            *counts.entry(gene.to_owned()).or_default() += 1;
        }
    }
    counts
}

/// Where a data unit is fetched from. The unit path from the catalogue is
/// appended to this.
pub const OAS_UNIT_BASE: &str = "https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/";

/// The committed list of data units this dataset is built from.
///
/// Compiled in rather than read from disk, so the build cannot be run against
/// a catalogue that is not the one in git. OAS has no stable URL for its unit
/// index -- the index comes out of a POST to a search form -- so there is no
/// upstream catalogue to pin in `sources.lock` the way the CATH release pins
/// the burial dataset's file list. This file is the catalogue, and reviewing
/// it means reading a diff.
pub const CATALOGUE: &str = include_str!("../../oas_units.tsv");

/// Most sequences kept from any one donor.
///
/// A cap rather than a sample of the whole: repertoire depth varies by three
/// orders of magnitude between units, and without a cap two or three deeply
/// sequenced donors would supply most of the dataset and the head would learn
/// them rather than the species. Which sequences are kept is decided by the
/// hash of the sequence, so it does not depend on file order.
pub const MAX_PER_DONOR: usize = 5_000;

/// One row of the committed catalogue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CatalogueRow {
    pub group: String,
    pub purpose: String,
    pub organism: String,
    pub study: String,
    /// The donor key this unit was selected under, as [`donor_key`] computes it.
    ///
    /// Recorded in the file so that which donors a catalogue covers can be
    /// read without downloading a gigabyte, and checked against the unit's own
    /// metadata at build time so the recorded value cannot drift from the
    /// truth. Two catalogues are kept disjoint by comparing this column.
    pub donor: String,
    pub isotype: String,
    pub unit: String,
    /// SHA-256 of the unit, when it has been pinned.
    ///
    /// `None` for a unit that has been chosen but not yet downloaded and
    /// hashed. The build records what it actually got either way, and counts
    /// the unpinned ones so the card cannot imply they were verified.
    pub sha256: Option<String>,
}

impl CatalogueRow {
    /// The cache filename for this unit: its path, flattened.
    pub fn file_name(&self) -> String {
        self.unit.replace("/csv/", "__").replace('/', "__")
    }

    pub fn url(&self) -> String {
        format!("{OAS_UNIT_BASE}{}", self.unit)
    }
}

fn catalogue_column(header: &csv::StringRecord, name: &str) -> Result<usize> {
    header
        .iter()
        .position(|field| field.trim() == name)
        .with_context(|| format!("oas_units.tsv has no column {name:?}"))
}

/// Parse the committed catalogue.
///
/// # Errors
/// If a column is missing or a row is short.
pub fn read_catalogue(text: &str) -> Result<Vec<CatalogueRow>> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').comment(Some(b'#')).from_reader(text.as_bytes());
    let header = rdr.headers()?.clone();
    let (c_group, c_purpose, c_org, c_study, c_donor, c_iso, c_unit, c_sha) = (
        catalogue_column(&header, "group")?,
        catalogue_column(&header, "purpose")?,
        catalogue_column(&header, "organism")?,
        catalogue_column(&header, "study")?,
        catalogue_column(&header, "donor")?,
        catalogue_column(&header, "isotype")?,
        catalogue_column(&header, "unit")?,
        catalogue_column(&header, "sha256")?,
    );
    let mut rows = Vec::new();
    for record in rdr.records() {
        let record = record?;
        let get = |i: usize| record.get(i).unwrap_or("").trim().to_owned();
        let sha = get(c_sha);
        rows.push(CatalogueRow {
            group: get(c_group),
            purpose: get(c_purpose),
            organism: get(c_org),
            study: get(c_study),
            donor: get(c_donor),
            isotype: get(c_iso),
            unit: get(c_unit),
            // "-" is how the file spells "not pinned", so a blank field and a
            // missing one read the same way.
            sha256: (!sha.is_empty() && sha != "-").then_some(sha),
        });
    }
    if rows.is_empty() {
        bail!("oas_units.tsv lists no units");
    }
    Ok(rows)
}

/// One sequence, with everything the split and the card need.
#[derive(Debug, Clone)]
pub struct Record {
    pub sequence: String,
    pub cohort: Cohort,
    pub organism: String,
    pub study: String,
    pub donor: String,
    pub isotype: String,
    pub v_call: String,
    pub j_call: String,
    pub cdr3: String,
    pub redundancy: u32,
}

/// Keep at most [`MAX_PER_DONOR`] from each donor, by sequence hash.
pub fn cap_per_donor(mut records: Vec<Record>, cap: usize) -> (Vec<Record>, u64) {
    let mut by_donor: BTreeMap<String, Vec<Record>> = BTreeMap::new();
    for r in records.drain(..) {
        by_donor.entry(r.donor.clone()).or_default().push(r);
    }
    let mut kept = Vec::new();
    let mut dropped = 0_u64;
    for (_, mut group) in by_donor {
        if group.len() > cap {
            // By hash of the sequence: stable, and unrelated to the order the
            // rows happened to be read in.
            group.sort_by_cached_key(|r| crate::sha256_hex(r.sequence.as_bytes()));
            dropped += (group.len() - cap) as u64;
            group.truncate(cap);
        }
        kept.extend(group);
    }
    kept.sort_by(|a, b| (&a.donor, &a.sequence).cmp(&(&b.donor, &b.sequence)));
    (kept, dropped)
}

/// Build a labelled cohort table. `is_human` is a plain boolean here because
/// every row in it has an answer.
pub fn to_labelled_table(rows: &[&Record], group_of: &BTreeMap<String, String>) -> Table {
    let mut t = base_columns(rows);
    t.push(
        "is_human",
        "True if the sequence came from a human repertoire, false if from a mouse one. The study that deposited it recorded this; it is not a measurement of the molecule and not an immunogenicity label",
        ColumnData::Bool(
            rows.iter()
                .map(|r| r.cohort.is_human().expect("a labelled row has a label"))
                .collect(),
        ),
    );
    t.push(
        "donor_group",
        "Donors joined where they share an identical rearrangement (same V gene, same J gene, byte-identical CDR3). The smallest unit that can be held out without leaking a lineage. All mouse donors fall in one group, which is why this dataset has no split",
        ColumnData::Utf8(rows.iter().map(|r| group_of[&r.donor].clone()).collect()),
    );
    t
}

/// Build an unlabelled cohort table. Deliberately has no `is_human` column:
/// these rows have no answer, and a column that could hold a null would get
/// averaged.
pub fn to_held_out_table(rows: &[&Record]) -> Table {
    base_columns(rows)
}

/// Donors joined where they share an identical rearrangement.
///
/// "Identical" means the same V gene call, the same J gene call and a
/// byte-identical CDR3 -- one rearrangement, seen in two animals. This is the
/// strictest honest reading of a shared lineage, chosen over sequence
/// clustering because clustering the whole V(D)J region at any workable
/// identity threshold groups germline families instead: at 90% it puts every
/// labelled donor in one component.
///
/// Rows with an empty CDR3 join nothing: an absent junction is not a shared
/// one.
pub fn group_donors_by_shared_rearrangement(rows: &[&Record]) -> BTreeMap<String, String> {
    let mut uf = UnionFind::new();
    for r in rows {
        // Every donor is a node even if it shares nothing, or a donor with no
        // usable CDR3 would be missing from the result.
        uf.union(&format!("D:{}", r.donor), &format!("D:{}", r.donor));
        if !r.cdr3.is_empty() {
            uf.union(&format!("D:{}", r.donor), &format!("R:{}|{}|{}", r.v_call, r.j_call, r.cdr3));
        }
    }
    let component_of = uf.components();
    rows.iter().map(|r| (r.donor.clone(), component_of[&format!("D:{}", r.donor)].clone())).collect()
}

fn base_columns(rows: &[&Record]) -> Table {
    let mut t = Table::new();
    t.push(
        "sequence",
        "Heavy-chain variable region, amino acids, as OAS aligned it (`sequence_alignment_aa`)",
        ColumnData::Utf8(rows.iter().map(|r| r.sequence.clone()).collect()),
    );
    t.push(
        "organism",
        "OAS `Species` verbatim, including the strain. Names the animal, not the immunoglobulin: HIS-Mouse and Kymouse are mice with human V genes",
        ColumnData::Utf8(rows.iter().map(|r| r.organism.clone()).collect()),
    );
    t.push(
        "study",
        "OAS study the unit belongs to",
        ColumnData::Utf8(rows.iter().map(|r| r.study.clone()).collect()),
    );
    t.push(
        "donor",
        "Donor key. `person:<subject>` when the unit names a subject, which is global because one donor is deposited by several studies; `study:<study>` when it does not",
        ColumnData::Utf8(rows.iter().map(|r| r.donor.clone()).collect()),
    );
    t.push(
        "isotype",
        "OAS `Isotype` as rendered to a string",
        ColumnData::Utf8(rows.iter().map(|r| r.isotype.clone()).collect()),
    );
    t.push(
        "v_call",
        "IMGT V gene call from OAS. Reported so the germline overlap between cohorts can be checked; not a feature",
        ColumnData::Utf8(rows.iter().map(|r| r.v_call.clone()).collect()),
    );
    t.push(
        "j_call",
        "IMGT J gene call from OAS",
        ColumnData::Utf8(rows.iter().map(|r| r.j_call.clone()).collect()),
    );
    t.push(
        "cdr3",
        "CDR3 amino acids as OAS delimited them",
        ColumnData::Utf8(rows.iter().map(|r| r.cdr3.clone()).collect()),
    );
    t.push(
        "redundancy",
        "Identical reads OAS collapsed into this row",
        ColumnData::Int32(rows.iter().map(|r| i32::try_from(r.redundancy).unwrap_or(i32::MAX)).collect()),
    );
    t
}

pub const SOURCES: &[Source] = &[sources::OAS];

pub fn build(ctx: &mut BuildContext) -> Result<()> {
    let catalogue = read_catalogue(CATALOGUE)?;
    let catalogue_digest = crate::sha256_hex(CATALOGUE.as_bytes());
    eprintln!("{} units in oas_units.tsv (sha256 {catalogue_digest})", catalogue.len());

    let mut records: Vec<Record> = Vec::new();
    let mut stats = ReadStats::default();
    let mut unpinned = 0_u64;
    let mut digests: BTreeMap<String, String> = BTreeMap::new();
    for row in &catalogue {
        let (path, sha) = ctx.cache.fetch_bulk("oas", &row.file_name(), &row.url())?;
        match &row.sha256 {
            Some(pinned) if pinned != &sha => bail!(
                "{} changed upstream: oas_units.tsv pins sha256 {pinned}, the file is {sha}.\n\
                 Delete the cached file to re-download, or update the catalogue deliberately.",
                row.unit
            ),
            Some(_) => {}
            None => unpinned += 1,
        }
        digests.insert(row.unit.clone(), sha);

        let (metadata, sequences, unit_stats) =
            read_unit(open_unit(&path)?).with_context(|| format!("reading {}", path.display()))?;
        accumulate(&mut stats, &unit_stats);

        // The catalogue and the file must agree about what the file is. A
        // mismatch means the wrong unit was pinned, which no later check
        // would catch: the sequences would simply be labelled as whatever the
        // catalogue claimed.
        if metadata.species != row.organism {
            bail!("{} is {} but oas_units.tsv says {}", row.unit, metadata.species, row.organism);
        }
        let cohort = cohort_of(&metadata.species)?;
        let donor = donor_key(&metadata, &row.study);
        // The donor column is what keeps this catalogue disjoint from the
        // humanness table's, and that check is made without downloading
        // either file's units. So the column has to be the truth: compare it
        // with what the unit's own metadata says, every build.
        if donor != row.donor {
            bail!(
                "{} is donor {donor} but oas_units.tsv says {}. The donor column is what keeps \
                 this catalogue disjoint from oas_human_donors.tsv, so it cannot be stale.",
                row.unit,
                row.donor
            );
        }
        records.extend(sequences.into_iter().map(|s| Record {
            sequence: s.sequence,
            cohort,
            organism: metadata.species.clone(),
            study: row.study.clone(),
            donor: donor.clone(),
            isotype: metadata.isotype.clone(),
            v_call: s.v_call,
            j_call: s.j_call,
            cdr3: s.cdr3,
            redundancy: s.redundancy,
        }));
    }

    // A sequence seen under two donors is not evidence about either, and
    // keeping it would put the same row on both sides of a split.
    let before = records.len() as u64;
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for (i, r) in records.iter().enumerate() {
        seen.entry(r.sequence.clone()).or_insert(i);
    }
    let shared: BTreeMap<String, ()> = {
        let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
        for r in &records {
            *counts.entry(r.sequence.as_str()).or_default() += 1;
        }
        counts.into_iter().filter(|(_, n)| *n > 1).map(|(s, _)| (s.to_owned(), ())).collect()
    };
    records.retain(|r| !shared.contains_key(&r.sequence));
    let shared_dropped = before - records.len() as u64;

    let (records, over_cap) = cap_per_donor(records, MAX_PER_DONOR);
    let labelled: Vec<&Record> = records.iter().filter(|r| r.cohort.is_human().is_some()).collect();
    if labelled.is_empty() {
        bail!("no labelled sequences: the catalogue has no human or mouse units");
    }

    // Group donors by shared rearrangement, and report what that grouping
    // does. This is where the dataset's shape was decided: the grouping is
    // computed every build so the card's claim about it is measured rather
    // than remembered. See the module documentation.
    let group_of = group_donors_by_shared_rearrangement(&labelled);
    let mut weights: BTreeMap<&str, u64> = BTreeMap::new();
    for r in &labelled {
        *weights.entry(group_of[&r.donor].as_str()).or_default() += 1;
    }
    let largest = weights.values().copied().max().unwrap_or(0);
    let share = largest as f64 / labelled.len() as f64;
    eprintln!(
        "{} donor groups over {} labelled rows; largest holds {largest} ({:.1}%)",
        weights.len(),
        labelled.len(),
        share * 100.0
    );

    // Which cohorts the groups straddle, because that is the fact that makes
    // a split impossible rather than merely lopsided. If some group held both
    // human and mouse rows, no split could separate the two classes at all;
    // as it stands the mouse rows are one group, which is nearly as bad and
    // has to be stated the same way.
    let mut groups_per_cohort: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for r in &labelled {
        let cohort = if r.cohort == Cohort::Human { "human" } else { "mouse" };
        groups_per_cohort.entry(cohort).or_default().insert(group_of[&r.donor].as_str());
    }
    let mixed = weights
        .keys()
        .filter(|g| groups_per_cohort.values().filter(|set| set.contains(**g)).count() > 1)
        .count();
    if mixed > 0 {
        bail!(
            "{mixed} donor group(s) hold both human and mouse rows, which means an identical \
             rearrangement was recorded in both species. That is either a mislabelled unit or a \
             V/J naming collision, and it has to be understood before these cohorts are used."
        );
    }

    let cohort_rows =
        |cohort: Cohort| -> Vec<&Record> { records.iter().filter(|r| r.cohort == cohort).collect() };
    let humans = cohort_rows(Cohort::Human);
    let mice = cohort_rows(Cohort::Mouse);
    let diagnostic = cohort_rows(Cohort::HumanLociInMouse);
    let out_of_distribution = cohort_rows(Cohort::OutOfDistribution);
    // Every table named, none of them a split. `human` and `mouse` carry
    // `is_human`; the other two deliberately do not.
    let extra = vec![
        ("human", to_labelled_table(&humans, &group_of)),
        ("mouse", to_labelled_table(&mice, &group_of)),
        ("diagnostic", to_held_out_table(&diagnostic)),
        ("out_of_distribution", to_held_out_table(&out_of_distribution)),
    ];
    let cohorts: Vec<(&'static str, u64)> =
        extra.iter().map(|(name, table)| (*name, table.n_rows() as u64)).collect();
    let summary = SplitSummary::default();

    // How much of the germline catalogue the two cohorts share. The first
    // thing to check when a species classifier looks too good: if the V calls
    // do not overlap at all, the task is germline lookup.
    let v_of = |cohort: Cohort| -> BTreeMap<String, u64> {
        v_gene_counts(labelled.iter().copied().filter(|r| r.cohort == cohort).map(|r| r.v_call.as_str()))
    };
    let human_v = v_of(Cohort::Human);
    let mouse_v = v_of(Cohort::Mouse);
    let shared_v = human_v.keys().filter(|g| mouse_v.contains_key(*g)).count();

    let donors = |rows: &[&Record]| rows.iter().map(|r| r.donor.as_str()).collect::<BTreeSet<_>>().len();
    let studies = |rows: &[&Record]| rows.iter().map(|r| r.study.as_str()).collect::<BTreeSet<_>>().len();

    let statistics = vec![
        ("OAS unit catalogue (oas_units.tsv) sha256".into(), catalogue_digest.clone()),
        ("data units read".into(), catalogue.len().to_string()),
        ("data units not pinned in the catalogue".into(), unpinned.to_string()),
        ("rows read from units".into(), stats.rows.to_string()),
        ("dropped: not productive".into(), stats.not_productive.to_string()),
        ("dropped: stop codon".into(), stats.stop_codon.to_string()),
        ("dropped: V-J out of frame".into(), stats.out_of_frame.to_string()),
        ("dropped: not a heavy chain".into(), stats.wrong_locus.to_string()),
        ("dropped: empty amino-acid alignment".into(), stats.empty_sequence.to_string()),
        ("dropped: non-standard residue or alignment gap".into(), stats.nonstandard_residue.to_string()),
        (
            format!("dropped: outside {MIN_RESIDUES}-{MAX_RESIDUES} residues"),
            stats.out_of_length_band.to_string(),
        ),
        ("dropped: duplicate within its unit".into(), stats.duplicate_within_unit.to_string()),
        ("dropped: identical sequence under more than one donor".into(), shared_dropped.to_string()),
        (format!("dropped: over the {MAX_PER_DONOR}-per-donor cap"), over_cap.to_string()),
        ("labelled rows (human and mouse cohorts)".into(), labelled.len().to_string()),
        (
            "  human".into(),
            format!("{} rows, {} donors, {} studies", humans.len(), donors(&humans), studies(&humans)),
        ),
        (
            "  mouse".into(),
            format!("{} rows, {} donors, {} studies", mice.len(), donors(&mice), studies(&mice)),
        ),
        (
            "diagnostic rows (mice with human Ig loci, unlabelled)".into(),
            format!("{} rows, {} donors", diagnostic.len(), donors(&diagnostic)),
        ),
        (
            "out-of-distribution rows (other species, unlabelled)".into(),
            format!("{} rows, {} donors", out_of_distribution.len(), donors(&out_of_distribution)),
        ),
        ("donor groups (donors joined by a shared rearrangement)".into(), weights.len().to_string()),
        ("  among human donors".into(), groups_per_cohort["human"].len().to_string()),
        ("  among mouse donors".into(), groups_per_cohort["mouse"].len().to_string()),
        ("largest donor group, share of labelled rows".into(), format!("{:.2}%", share * 100.0)),
        ("V genes in human rows".into(), human_v.len().to_string()),
        ("V genes in mouse rows".into(), mouse_v.len().to_string()),
        ("V genes in both".into(), shared_v.to_string()),
    ];

    let construction = format!(
        "1. Read the {} data units listed in `oas_units.tsv`, which is committed to the repository \
            and compiled into the builder. OAS has no stable URL for its unit index, so this file \
            is the catalogue; its own SHA-256 is in the statistics above and each unit's is in \
            `manifest.json`. A unit whose digest does not match the catalogue stops the build.\n\
         2. Keep rows OAS already marked productive, in frame and free of stop codons, on the heavy \
            chain, whose amino-acid alignment is {MIN_RESIDUES}-{MAX_RESIDUES} standard residues.\n\
         3. Drop any sequence seen under more than one donor, and cap each donor at \
            {MAX_PER_DONOR} sequences chosen by sequence hash.\n\
         4. Write four cohorts. `human` and `mouse` carry `is_human`; HIS-Mouse and Kymouse carry \
            human immunoglobulin loci and go to `diagnostic` with no label, where answering \
            \"human\" is correct; rat, rabbit, rhesus and bactrian camel go to \
            `out_of_distribution`, also unlabelled.\n\
         5. Join donors wherever they share an identical rearrangement -- the same V gene, the \
            same J gene and a byte-identical CDR3 -- and record the resulting group in \
            `donor_group`. That is the smallest unit that can be held out without carrying a \
            lineage across the boundary.\n\n\
         **There is no train/validation/test split, because this data does not admit one.** All \
         {} mouse donors fall into a single donor group, spanning three inbred strains and nine \
         studies, because C57BL/6 and BALB/c animals are genetically near-identical and their \
         repertoires share public rearrangements. One laboratory mouse is not an independent \
         sample of its species the way one person is. The largest group therefore holds {:.1}% \
         of the labelled rows, and any grouped split leaves the whole mouse cohort on one side \
         -- which is not a split.\n\n\
         Sequence clustering does not rescue it: clustering the amino-acid V(D)J region at 90% \
         identity and joining donors through shared clusters put all {} labelled donors into one \
         component. That region is mostly germline-constant, so the threshold groups V-gene \
         families rather than lineages. The shared-rearrangement grouping is the strict \
         alternative, and `donor_group` records it: the human side can be split on it ({} \
         groups), the mouse side cannot ({}).",
        catalogue.len(),
        donors(&mice),
        share * 100.0,
        donors(&labelled),
        groups_per_cohort["human"].len(),
        groups_per_cohort["mouse"].len()
    );

    let card = CardInput {
        pretty_name: "OpenDNA Antibody (human and mouse heavy-chain cohorts, OAS)",
        license: CardLicense::Known("cc-by-4.0"),
        license_note: "OAS is licensed CC BY 4.0 by its maintainers, stated on \
             <https://opig.stats.ox.ac.uk/webapps/oas/>: \"the data contained within OAS is \
             available under a CC-BY 4.0 license\". Redistribution is permitted with attribution, \
             so this dataset is CC BY 4.0. Cite Olsen, Boyles and Deane (2021) and the study each \
             row's `study` column names.",
        tags: &["biology", "protein", "antibody", "immunoglobulin", "repertoire"],
        summary: format!(
            "Heavy-chain variable regions from human and mouse antibody repertoires in the \
             Observed Antibody Space, labelled by the species of the donor: {} rows from {} \
             donors across {} studies. Two unlabelled cohorts travel with them -- mice carrying \
             human immunoglobulin loci, where a humanness measure should answer human, and four \
             species with one or two studies each.\n\nThese are evaluation cohorts, not a \
             training split. There is no train/validation/test division, because the mouse donors \
             cannot be divided: all {} of them share identical rearrangements and form one \
             group. Construction says why in full, and `donor_group` carries the grouping so the \
             human side can still be split by anyone who needs to.",
            labelled.len(),
            donors(&labelled),
            studies(&labelled),
            donors(&mice)
        ),
        sources: SOURCES,
        excluded: vec![&sources::SABDAB_EXCLUDED],
        construction,
        schema: &extra[0].1,
        row_unit: "sequences",
        group_unit: "donor groups",
        splits: &summary,
        cohorts: &cohorts,
        statistics,
        citations: ANTIBODY_CITATIONS,
        limitations: ANTIBODY_LIMITATIONS,
        lock: ctx.cache.lock(),
    };
    write_dataset(
        ctx,
        FinishedDataset {
            name: NAME,
            tables: &[],
            extra: &extra,
            card,
            // Nothing was split and nothing was clustered, so the manifest
            // records neither rather than a fraction that divided nothing.
            fractions: None,
            cluster_params: None,
        },
    )?;
    Ok(())
}

// The card renderer wraps this in a bibtex fence, so this must be bibtex and
// nothing else. The sentence that used to follow it -- that each row's `study`
// column names the paper that produced it -- is in the licence note, where it
// is not inside a code block.
const ANTIBODY_CITATIONS: &str = "\
@article{olsen2022oas,
  title   = {Observed Antibody Space: A diverse database of cleaned, annotated,
             and translated unpaired and paired antibody sequences},
  author  = {Olsen, Tobias H. and Boyles, Fergus and Deane, Charlotte M.},
  journal = {Protein Science},
  volume  = {31},
  number  = {1},
  pages   = {141--146},
  year    = {2022},
  doi     = {10.1002/pro.4205}
}";

const ANTIBODY_LIMITATIONS: &[&str] = &[
    "The label is the species of the donor, recorded by the depositing study. It is not a \
     measurement of the molecule, and it is not an immunogenicity label: a sequence resembling a \
     human repertoire is not thereby safe in humans, and the converse does not hold either.",
    "Human and mouse V genes differ almost completely -- the card reports how many appear in \
     both cohorts -- so separating the two species is largely a germline lookup. Measured on \
     held-out rows, a lookup table on the V gene call alone reaches 98.5% against a 68.9% \
     majority baseline. An accuracy on this task is therefore not evidence that a model learned \
     anything about antibodies.",
    "Heavy chains only. Light chains have their own germline repertoire and are a separate task.",
    "The out-of-distribution cohort is one study per species for rat, rabbit and bactrian camel, \
     and two for rhesus. A score on those species cannot be separated from one laboratory's \
     protocol, which is why those rows carry no label.",
    "These are cohorts, not folds. There is no train/test split and none can be made: all 18 \
     mouse donors share identical rearrangements and form a single group, so any grouped split \
     puts the whole mouse cohort on one side. The `donor_group` column carries the grouping; the \
     human cohort has 15 groups and can be divided, the mouse cohort has one and cannot.",
    "Donor keys merge a named subject across studies, because one donor is deposited by several \
     studies -- CAP256 appears in three. Where two studies reuse a placeholder such as Donor-1 \
     for different people, this merges them too, which costs data rather than leaking it.",
    "Repertoires are not samples of what a B cell can make. They are shaped by the primers used, \
     the cells sorted, the immunisation, and the disease state of the donor, all of which are in \
     the columns and none of which is controlled for here.",
    "A cap of 5,000 sequences per donor keeps deeply sequenced donors from dominating. It also \
     means row counts per donor say nothing about repertoire size.",
];

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal OAS unit: the JSON metadata line, the AIRR header, and rows.
    fn unit(rows: &[&str]) -> String {
        let metadata = r#"{"Run": "ERR1", "Author": "Greiff et al., 2017", "Species": "mouse_C57BL/6", "Subject": "HP-HEL-Mouse-1", "Chain": "Heavy", "Isotype": "Bulk", "Disease": "None", "BType": "Plasma-B-Cells"}"#;
        let mut text = String::new();
        // The metadata line is a single CSV field, so its quotes are doubled.
        text.push('"');
        text.push_str(&metadata.replace('"', "\"\""));
        text.push_str("\"\n");
        text.push_str("sequence_alignment_aa,productive,stop_codon,vj_in_frame,locus,v_call,j_call,cdr3_aa,Redundancy\n");
        for row in rows {
            text.push_str(row);
            text.push('\n');
        }
        text
    }

    /// 117 residues, the length of the real first row in the fixture unit.
    fn good_sequence() -> String {
        "QVQLQQPGAELVRPGSSVKLSCKASGYTFTSYWMHWVKQRPIQGLEWIGNIDPSDSETHYNQKFKDKATLTVDKSSSTAYMQLSSLTSEDSAVYSVALLRWYPACYAYWGRGNLVTVS".to_string()
    }

    fn row(sequence: &str, productive: &str, stop: &str, frame: &str, locus: &str) -> String {
        format!("{sequence},{productive},{stop},{frame},{locus},IGHV1-52*01,IGHJ3*01,ALLRWYPACYAY,3")
    }

    #[test]
    fn a_unit_carries_its_own_organism_and_donor() {
        let text = unit(&[&row(&good_sequence(), "T", "F", "T", "H")]);
        let (metadata, sequences, stats) = read_unit(text.as_bytes()).unwrap();
        assert_eq!(metadata.species, "mouse_C57BL/6");
        assert_eq!(metadata.subject, "HP-HEL-Mouse-1");
        assert_eq!(metadata.chain, "Heavy");
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].redundancy, 3);
        assert_eq!(sequences[0].cdr3, "ALLRWYPACYAY");
        assert_eq!(stats.kept, 1);
        assert_eq!(stats.rows, 1);
    }

    /// Every rejection is counted, because a card that cannot account for the
    /// difference between rows read and rows kept is a card nobody can check.
    #[test]
    fn every_kind_of_rejection_is_counted() {
        let long = "A".repeat(MAX_RESIDUES + 1);
        let short = "A".repeat(MIN_RESIDUES - 1);
        let text = unit(&[
            &row(&good_sequence(), "T", "F", "T", "H"),
            &row(&good_sequence(), "F", "F", "T", "H"),
            &row(&good_sequence(), "T", "T", "T", "H"),
            &row(&good_sequence(), "T", "F", "F", "H"),
            &row(&good_sequence(), "T", "F", "T", "K"),
            &row("", "T", "F", "T", "H"),
            &row("QVQLXQPGAELVRPGSSVKLSCKASGYTFTSYWMHWVKQRPIQGLEWIGNIDPSDSETHYNQKFKDKATLTVDKSSSTAYMQLSSLTSEDSAVYSVALLRWYPACYAYWGRGNLVTVS", "T", "F", "T", "H"),
            &row(&long, "T", "F", "T", "H"),
            &row(&short, "T", "F", "T", "H"),
        ]);
        let (_, sequences, stats) = read_unit(text.as_bytes()).unwrap();
        assert_eq!(stats.rows, 9);
        assert_eq!(stats.kept, 1);
        assert_eq!(sequences.len(), 1);
        assert_eq!(stats.not_productive, 1);
        assert_eq!(stats.stop_codon, 1);
        assert_eq!(stats.out_of_frame, 1);
        assert_eq!(stats.wrong_locus, 1);
        assert_eq!(stats.empty_sequence, 1);
        assert_eq!(stats.nonstandard_residue, 1, "an X is not a residue to embed");
        assert_eq!(stats.out_of_length_band, 2);
        // Nothing may vanish between the two totals.
        let accounted = stats.not_productive
            + stats.stop_codon
            + stats.out_of_frame
            + stats.wrong_locus
            + stats.empty_sequence
            + stats.nonstandard_residue
            + stats.out_of_length_band
            + stats.duplicate_within_unit
            + stats.kept;
        assert_eq!(accounted, stats.rows);
    }

    #[test]
    fn an_identical_sequence_is_kept_once() {
        let text =
            unit(&[&row(&good_sequence(), "T", "F", "T", "H"), &row(&good_sequence(), "T", "F", "T", "H")]);
        let (_, sequences, stats) = read_unit(text.as_bytes()).unwrap();
        assert_eq!(sequences.len(), 1);
        assert_eq!(stats.duplicate_within_unit, 1);
    }

    /// The trap this dataset exists to avoid. A transgenic mouse's V genes are
    /// human, so it must not be labelled non-human -- and it must not be
    /// labelled human either, since the animal is a mouse. It gets no label.
    #[test]
    fn a_mouse_with_human_loci_is_never_given_a_species_label() {
        for species in ["HIS-Mouse", "Kymouse"] {
            let cohort = cohort_of(species).unwrap();
            assert_eq!(cohort, Cohort::HumanLociInMouse, "{species}");
            assert_eq!(cohort.is_human(), None, "{species} must carry no label");
            assert_eq!(cohort.table(), "diagnostic");
        }
    }

    #[test]
    fn ordinary_mice_and_humans_are_labelled_and_trained_on() {
        assert_eq!(cohort_of("human").unwrap().is_human(), Some(true));
        for strain in ["mouse", "mouse_C57BL/6", "mouse_BALB/c", "mouse_Swiss-Webster"] {
            assert_eq!(cohort_of(strain).unwrap().is_human(), Some(false), "{strain}");
            assert_eq!(cohort_of(strain).unwrap().table(), "labelled");
        }
    }

    /// Rat, rabbit and camel have one study each. They are kept, unlabelled,
    /// so a head can be reported against them without being trained on one
    /// laboratory's protocol.
    #[test]
    fn single_study_species_are_held_out_rather_than_trained_on() {
        for species in ["rat", "rabbit", "rhesus", "bactrian_camel"] {
            let cohort = cohort_of(species).unwrap();
            assert_eq!(cohort, Cohort::OutOfDistribution, "{species}");
            assert_eq!(cohort.is_human(), None);
            assert_eq!(cohort.table(), "out_of_distribution");
        }
    }

    /// A species nobody has looked at must stop the build. The alternative is
    /// a catch-all arm that quietly decides for us.
    #[test]
    fn an_unrecognised_species_is_an_error_not_a_guess() {
        let error = cohort_of("llama").unwrap_err().to_string();
        assert!(error.contains("llama"), "{error}");
        assert!(error.contains("must not be labelled non-human"), "{error}");
    }

    /// CAP256 is one HIV donor deposited by three studies. Keying by name
    /// keeps the three on the same side of every split.
    #[test]
    fn one_donor_in_several_studies_gets_one_key() {
        let metadata = |subject: &str| UnitMetadata {
            species: "human".into(),
            subject: subject.into(),
            chain: "Heavy".into(),
            author: "x".into(),
            isotype: String::new(),
            disease: String::new(),
            b_type: String::new(),
        };
        assert_eq!(
            donor_key(&metadata("CAP256"), "Bhiman_2015"),
            donor_key(&metadata("CAP256"), "Doria-Rose_2015")
        );
    }

    /// An absent subject is a missing field, not a person. Those units fall
    /// back to the study, so the boundary is still safe.
    #[test]
    fn a_missing_subject_falls_back_to_the_study() {
        let metadata = UnitMetadata {
            species: "rat".into(),
            subject: "no".into(),
            chain: "Heavy".into(),
            author: "x".into(),
            isotype: String::new(),
            disease: String::new(),
            b_type: String::new(),
        };
        assert_eq!(donor_key(&metadata, "VanDuijn_2017"), "study:VanDuijn_2017");
        // And two studies with no subject stay apart.
        assert_ne!(donor_key(&metadata, "A"), donor_key(&metadata, "B"));
    }

    #[test]
    fn v_genes_are_counted_by_gene_not_allele() {
        let make = |v: &str| AntibodySequence {
            sequence: good_sequence(),
            v_call: v.into(),
            j_call: "IGHJ3*01".into(),
            cdr3: "AAA".into(),
            redundancy: 1,
        };
        let sequences = [make("IGHV1-52*01"), make("IGHV1-52*02"), make("IGHV3-23*01")];
        let counts = v_gene_counts(sequences.iter().map(|s| s.v_call.as_str()));
        assert_eq!(counts["IGHV1-52"], 2);
        assert_eq!(counts["IGHV3-23"], 1);
        assert_eq!(counts.len(), 2);
    }

    /// The fixtures above state what an OAS unit looks like. This checks that
    /// against real ones, which is the only way to find out that the shape
    /// was written down wrong -- a fixture agrees with whatever built it.
    ///
    /// Skipped unless `OPENDNA_OAS_CACHE` points at a directory of downloaded
    /// `.csv.gz` units, because the units are not in the repository and a
    /// test that fails on every machine without them gets switched off.
    #[test]
    fn real_units_parse_and_every_species_is_recognised() {
        let Some(dir) = std::env::var_os("OPENDNA_OAS_CACHE") else {
            println!("OPENDNA_OAS_CACHE unset; skipping the real-unit check");
            return;
        };
        let mut checked = 0_usize;
        let mut cohorts: BTreeMap<String, usize> = BTreeMap::new();
        for entry in std::fs::read_dir(&dir).expect("OPENDNA_OAS_CACHE is not readable").flatten() {
            let path = entry.path();
            if path.extension().is_none_or(|e| e != "gz") {
                continue;
            }
            let (metadata, sequences, stats) = read_unit(open_unit(&path).expect("opening a cached unit"))
                .unwrap_or_else(|e| panic!("{} is not readable as an OAS unit: {e:#}", path.display()));
            let cohort = cohort_of(&metadata.species).unwrap_or_else(|e| panic!("{}: {e:#}", path.display()));
            *cohorts.entry(format!("{cohort:?}")).or_default() += 1;

            assert_eq!(metadata.chain, "Heavy", "{} is not a heavy chain", path.display());
            assert!(stats.rows > 0, "{} has no rows", path.display());
            assert_eq!(
                stats.kept as usize,
                sequences.len(),
                "{}: the kept count and the sequences disagree",
                path.display()
            );
            for s in &sequences {
                assert!(
                    (MIN_RESIDUES..=MAX_RESIDUES).contains(&s.sequence.len()),
                    "{}: a kept sequence is {} residues",
                    path.display(),
                    s.sequence.len()
                );
                assert!(is_standard_protein(&s.sequence), "{}: a kept sequence is not", path.display());
            }
            checked += 1;
        }
        assert!(checked > 0, "OPENDNA_OAS_CACHE held no .csv.gz units");
        println!("{checked} real units parsed; cohorts: {cohorts:?}");
    }

    /// Upstream types these fields two ways. Both are real: the strings come
    /// from Greiff 2017 and the list from Li 2017, in the same release.
    #[test]
    fn a_metadata_field_typed_two_ways_upstream_is_read_either_way() {
        let read = |isotype: &str| -> UnitMetadata {
            let json = format!(
                r#"{{"Species": "human", "Subject": "S1", "Chain": "Heavy", "Author": "a", "Isotype": {isotype}}}"#
            );
            serde_json::from_str(&json).unwrap_or_else(|e| panic!("{isotype}: {e}"))
        };
        assert_eq!(read(r#""Bulk""#).isotype, "Bulk");
        assert_eq!(read(r#"["Bulk"]"#).isotype, "Bulk");
        assert_eq!(read(r#"["IGHG","IGHM"]"#).isotype, "IGHG;IGHM");
        assert_eq!(read("null").isotype, "");
        assert_eq!(read("5").isotype, "5");
        // An absent field is fine too.
        let bare: UnitMetadata =
            serde_json::from_str(r#"{"Species": "human", "Subject": "S1", "Chain": "Heavy", "Author": "a"}"#)
                .unwrap();
        assert_eq!(bare.isotype, "");
    }

    /// The fields that decide something stay strict: a `Species` that is not
    /// a string is a change upstream that must stop the build.
    #[test]
    fn a_load_bearing_field_is_still_strict() {
        let json = r#"{"Species": ["human"], "Subject": "S1", "Chain": "Heavy", "Author": "a"}"#;
        assert!(serde_json::from_str::<UnitMetadata>(json).is_err());
    }

    /// The committed catalogue is data, and data in a repository rots. This
    /// checks the invariants the build relies on, so a hand edit that breaks
    /// one fails here rather than three hours into a download.
    #[test]
    fn the_committed_catalogue_is_consistent() {
        let rows = read_catalogue(CATALOGUE).expect("oas_units.tsv must parse");
        assert_eq!(rows.len(), 54, "the catalogue should hold the 54 chosen units");

        let mut units = BTreeSet::new();
        let mut labelled_studies: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
        for row in &rows {
            assert!(units.insert(row.unit.clone()), "{} is listed twice", row.unit);
            // The organism must agree with the cohort the purpose implies, or
            // a unit would land in the wrong cohort, labelled or not.
            let cohort = cohort_of(&row.organism).expect(&row.organism);
            let expected = match cohort {
                Cohort::Human | Cohort::Mouse => "labelled cohort",
                Cohort::HumanLociInMouse => "humanisation diagnostic",
                Cohort::OutOfDistribution => "out of distribution",
            };
            assert_eq!(row.purpose, expected, "{} is {}", row.unit, row.organism);
            if cohort.is_human().is_some() {
                labelled_studies.entry(row.group.as_str()).or_default().insert(row.study.as_str());
            }
            assert!(row.unit.ends_with(".csv.gz"), "{}", row.unit);
            assert!(!row.unit.starts_with('/') && !row.unit.contains(".."), "{}", row.unit);
            assert!(row.url().starts_with(OAS_UNIT_BASE), "{}", row.url());
            assert!(!row.file_name().contains('/'), "{}", row.file_name());
        }
        // Both labelled cohorts need several independent studies, so that a
        // number computed on either is not one laboratory's protocol. This
        // says nothing about splitting them: the mouse donors cannot be
        // split at all, which the module documentation explains.
        for cohort in ["human", "mouse"] {
            let studies = labelled_studies.get(cohort).map_or(0, BTreeSet::len);
            assert!(studies >= 6, "{cohort} has only {studies} studies in the catalogue");
        }
    }

    #[test]
    fn a_catalogue_row_knows_where_it_comes_from_and_where_it_lands() {
        let text = "group\tpurpose\torganism\tstudy\tdonor\tisotype\tunique_sequences\tunit\tsha256\n\
                    human\ttrained on\thuman\tGalson_2015\tperson:S1\tIGHM\t76027\tGalson_2015/csv/SRR3990841_Heavy_IGHM.csv.gz\tab\n";
        let rows = read_catalogue(text).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].file_name(), "Galson_2015__SRR3990841_Heavy_IGHM.csv.gz");
        assert_eq!(
            rows[0].url(),
            "https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Galson_2015/csv/SRR3990841_Heavy_IGHM.csv.gz"
        );
        assert_eq!(rows[0].sha256.as_deref(), Some("ab"));
    }

    /// An unpinned unit must read as unpinned, not as a digest of "-".
    #[test]
    fn a_dash_means_not_yet_pinned() {
        let text = "group\tpurpose\torganism\tstudy\tdonor\tisotype\tunique_sequences\tunit\tsha256\n\
                    human\ttrained on\thuman\tS\tperson:S1\tIGHM\t1\tS/csv/u.csv.gz\t-\n";
        assert_eq!(read_catalogue(text).unwrap()[0].sha256, None);
    }

    #[test]
    fn an_empty_catalogue_is_an_error() {
        let text = "group\tpurpose\torganism\tstudy\tdonor\tisotype\tunique_sequences\tunit\tsha256\n";
        assert!(read_catalogue(text).unwrap_err().to_string().contains("lists no units"));
    }

    fn record(donor: &str, sequence: &str) -> Record {
        Record {
            sequence: sequence.to_owned(),
            cohort: Cohort::Human,
            organism: "human".into(),
            study: "S".into(),
            donor: donor.to_owned(),
            isotype: "IGHM".into(),
            v_call: "IGHV3-23*01".into(),
            j_call: "IGHJ4*02".into(),
            cdr3: "ARDL".into(),
            redundancy: 1,
        }
    }

    /// The cap applies per donor, not across the dataset, so a shallow donor
    /// is never truncated because a deep one exists.
    #[test]
    fn the_cap_is_per_donor_and_does_not_depend_on_order() {
        let mut records = Vec::new();
        for i in 0..10 {
            records.push(record("person:deep", &format!("SEQUENCE{i}")));
        }
        records.push(record("person:shallow", "ONLYONE"));
        let (kept, dropped) = cap_per_donor(records.clone(), 4);
        assert_eq!(dropped, 6);
        assert_eq!(kept.iter().filter(|r| r.donor == "person:deep").count(), 4);
        assert_eq!(kept.iter().filter(|r| r.donor == "person:shallow").count(), 1);

        // Which four is decided by the sequence hash, so shuffling the input
        // cannot change the dataset.
        records.reverse();
        let (again, _) = cap_per_donor(records, 4);
        let one: Vec<&String> = kept.iter().map(|r| &r.sequence).collect();
        let two: Vec<&String> = again.iter().map(|r| &r.sequence).collect();
        assert_eq!(one, two, "the cap must not depend on input order");
    }

    #[test]
    fn under_the_cap_nothing_is_dropped() {
        let records = vec![record("person:a", "AAA"), record("person:b", "BBB")];
        let (kept, dropped) = cap_per_donor(records, 5);
        assert_eq!(kept.len(), 2);
        assert_eq!(dropped, 0);
    }

    #[test]
    fn a_stream_that_is_not_an_oas_unit_is_rejected() {
        assert!(read_unit(b"".as_slice()).unwrap_err().to_string().contains("empty OAS unit"));
        let error = read_unit(b"not json\nheader\n".as_slice()).unwrap_err().to_string();
        assert!(error.contains("expected JSON"), "{error}");
    }
}
