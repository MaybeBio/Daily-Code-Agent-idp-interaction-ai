//! The peptide humanness prevalence table.
//!
//! Not a split dataset. There is nothing to train and nothing to hold out:
//! the output is a reference table of how many human donors hold each
//! 9-residue peptide, and a score computed from it is an exact count, not a
//! prediction. So this writes its own three files rather than going through
//! [`crate::context::write_dataset`], which assumes Parquet splits and a
//! schema.
//!
//! The measure itself lives in [`opendna_analysis::humanness`]. This module
//! only decides which repertoires go into the table, and records honestly
//! what that choice does to the numbers.
//!
//! ## What the donor choice costs
//!
//! Prevalence is a fraction of *donors*, so the donor count sets how fine a
//! threshold the table can express: from `n` donors it takes `n + 1` values,
//! and a "1% of donors" threshold means nothing at all unless `n >= 100`.
//! That is why the catalogue is 150 donors wide and one unit deep rather than
//! a handful of deeply sequenced repertoires.
//!
//! The price is paid in depth. Each unit holds a few thousand sequences out
//! of a repertoire of millions, so a peptide the donor really does have, but
//! which did not turn up in the sample, is counted as absent. Every
//! prevalence in this table is therefore biased *down*, and by an amount that
//! depends on how deeply each donor happened to be sequenced. Scores from
//! this table are internally comparable and are **not** comparable with
//! published OASis values, which were computed from deep repertoires of 266
//! subjects. The card says so where a reader will see it before the numbers.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::BufWriter;

use anyhow::{bail, Context, Result};
use opendna_analysis::humanness::{PrevalenceTable, Stringency, FORMAT, PEPTIDE_LEN};
use serde::Serialize;

use crate::context::BuildContext;
use crate::datasets::antibody::{self, open_unit, read_unit, ReadStats, MAX_RESIDUES, MIN_RESIDUES};
use crate::sources::{self, Source};

pub const NAME: &str = "humanness";

/// The committed list of human donor repertoires the table is built from.
///
/// Compiled in for the same reason as the antibody catalogue: OAS publishes no
/// stable URL for its unit index, so this file *is* the catalogue and changing
/// which donors define "human" has to show up in a reviewed diff.
pub const DONORS: &str = include_str!("../../oas_human_donors.tsv");

/// Peptides held by fewer donors than this are left out of the written table.
///
/// A peptide seen in exactly one donor's shallow sample is not evidence about
/// human repertoires in general, and the single-donor tail is the bulk of the
/// table by count. Dropping it costs nothing at the thresholds anyone scores
/// at -- one donor in 150 is 0.67%, below even the loosest 1% -- and the card
/// reports how many peptides it was, so the omission is visible rather than
/// implicit.
pub const MIN_DONORS: u32 = 2;

/// Fewest sequences a donor must contribute to count as a donor.
///
/// Prevalence counts the donors who hold a peptide, so a donor who
/// contributed almost nothing votes "absent" on almost every peptide and
/// drags the whole table down by `1/donors` each time. Two donors in this
/// catalogue contribute zero after filtering and two more contribute two and
/// three sequences; counting them would divide every prevalence by donors who
/// could not possibly have held the peptide.
///
/// The value is not tuned. Measured on the committed catalogue, depth per
/// donor is bimodal: nine donors fall below 500 sequences -- all from two
/// studies, one depositing partial-length sequences and one tiny samples --
/// and the next donor up has 1,153. The cut sits in that empty gap, so
/// anywhere between 500 and 1,100 gives the same 141 donors, which is what
/// makes it worth stating as a number rather than a judgement.
pub const MIN_SEQUENCES_PER_DONOR: usize = 500;

/// One row of the donor catalogue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DonorRow {
    pub study: String,
    /// Donor key, as [`antibody::donor_key`] computes it from unit metadata.
    pub donor: String,
    pub isotype: String,
    pub unique_sequences: u64,
    pub unit: String,
    pub sha256: Option<String>,
}

impl DonorRow {
    /// The cache filename for this unit: its path, flattened.
    pub fn file_name(&self) -> String {
        self.unit.replace("/csv/", "__").replace('/', "__")
    }

    pub fn url(&self) -> String {
        format!("{}{}", antibody::OAS_UNIT_BASE, self.unit)
    }
}

fn column(header: &csv::StringRecord, name: &str) -> Result<usize> {
    header
        .iter()
        .position(|field| field.trim() == name)
        .with_context(|| format!("oas_human_donors.tsv has no column {name:?}"))
}

/// Parse the committed donor catalogue.
pub fn read_donors(text: &str) -> Result<Vec<DonorRow>> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').comment(Some(b'#')).from_reader(text.as_bytes());
    let header = rdr.headers()?.clone();
    let (c_study, c_donor, c_iso, c_seqs, c_unit, c_sha) = (
        column(&header, "study")?,
        column(&header, "donor")?,
        column(&header, "isotype")?,
        column(&header, "unique_sequences")?,
        column(&header, "unit")?,
        column(&header, "sha256")?,
    );
    let mut rows = Vec::new();
    for record in rdr.records() {
        let record = record?;
        let get = |i: usize| record.get(i).unwrap_or("").trim().to_owned();
        let unit = get(c_unit);
        // A unit path becomes a URL and a filename, so it must not be able to
        // escape either.
        if unit.starts_with('/') || unit.contains("..") || !unit.ends_with(".csv.gz") {
            bail!("oas_human_donors.tsv has an implausible unit path {unit:?}");
        }
        let sha = get(c_sha);
        rows.push(DonorRow {
            study: get(c_study),
            donor: get(c_donor),
            isotype: get(c_iso),
            unique_sequences: get(c_seqs).parse().unwrap_or(0),
            unit,
            sha256: (!sha.is_empty() && sha != "-").then_some(sha),
        });
    }
    if rows.is_empty() {
        bail!("oas_human_donors.tsv lists no units");
    }
    Ok(rows)
}

/// Donor keys the antibody dataset uses, from its committed catalogue.
///
/// Read from the `donor` column rather than from the units themselves, which
/// would mean downloading a gigabyte to answer a question the catalogue
/// already records. The antibody build checks that column against each unit's
/// own metadata every time it runs, so a stale value there fails that build
/// rather than quietly weakening this check.
fn evaluation_donors() -> Result<BTreeSet<String>> {
    Ok(antibody::read_catalogue(antibody::CATALOGUE)?.into_iter().map(|row| row.donor).collect())
}

/// Numbers the card reports about the table it describes.
struct TableStats {
    donors: u32,
    units: usize,
    studies: usize,
    sequences: u64,
    read: ReadStats,
    distinct_peptides: usize,
    kept_peptides: usize,
    /// How many peptides occur in at least each stringency's fraction of donors.
    at_stringency: Vec<(Stringency, usize)>,
    /// Donors too shallow to count, with what each contributed.
    dropped: Vec<(String, usize)>,
    /// Fewest and most sequences kept from any one donor.
    ///
    /// The spread is the size of the unevenness in the sampling bias: a
    /// peptide is more likely to be seen at all in the donor who contributed
    /// four times as many sequences, so prevalence is pulled towards the
    /// deeper donors' repertoires. Reported because "biased downwards" on its
    /// own does not say by how much it varies.
    depth: (u64, u64),
    unpinned: u64,
    bytes: u64,
}

#[derive(Serialize)]
struct Manifest<'a> {
    dataset: &'a str,
    pipeline: &'static str,
    pipeline_version: &'static str,
    /// The format version of `table.bin`, so a reader can refuse an old one.
    table_format: u32,
    peptide_length: usize,
    donors: u32,
    min_donors: u32,
    distinct_peptides: usize,
    peptides_written: usize,
    sequences: u64,
    catalogue_sha256: String,
    table_sha256: String,
    table_bytes: u64,
    units: BTreeMap<&'a str, String>,
    /// Donors listed in the catalogue but too shallow to be counted, with
    /// the number of sequences each contributed.
    donors_below_depth_floor: &'a [(String, usize)],
    min_sequences_per_donor: usize,
}

pub const SOURCES: &[Source] = &[sources::OAS];

pub fn build(ctx: &mut BuildContext) -> Result<()> {
    let catalogue = read_donors(DONORS)?;
    let catalogue_digest = crate::sha256_hex(DONORS.as_bytes());
    eprintln!("{} donor units in oas_human_donors.tsv (sha256 {catalogue_digest})", catalogue.len());

    // A table built from donors it is later scored against measures how well
    // it remembers its own input. Checked before anything is downloaded,
    // because the answer does not depend on the files.
    let evaluation = evaluation_donors()?;
    let overlap: Vec<&str> =
        catalogue.iter().map(|r| r.donor.as_str()).filter(|d| evaluation.contains(*d)).collect();
    if !overlap.is_empty() {
        bail!(
            "{} donor(s) appear in both oas_human_donors.tsv and oas_units.tsv: {}. A table built \
             from the donors it is evaluated against measures recall of its own input, not \
             humanness.",
            overlap.len(),
            overlap.join(", ")
        );
    }

    let mut repertoires: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut read = ReadStats::default();
    let mut digests: BTreeMap<&str, String> = BTreeMap::new();
    let mut unpinned = 0_u64;
    for row in &catalogue {
        let (path, sha) = ctx.cache.fetch_bulk("oas", &row.file_name(), &row.url())?;
        match &row.sha256 {
            Some(pinned) if pinned != &sha => bail!(
                "{} changed upstream: oas_human_donors.tsv pins sha256 {pinned}, the file is \
                 {sha}.\nDelete the cached file to re-download, or update the catalogue \
                 deliberately.",
                row.unit
            ),
            Some(_) => {}
            None => unpinned += 1,
        }
        digests.insert(row.unit.as_str(), sha);

        let (metadata, sequences, unit_stats) =
            read_unit(open_unit(&path)?).with_context(|| format!("reading {}", path.display()))?;
        antibody::accumulate(&mut read, &unit_stats);

        // This table is the definition of "human". A non-human unit in it
        // would not fail anything later -- it would just quietly become part
        // of what human means.
        if metadata.species != "human" {
            bail!(
                "{} is {} but oas_human_donors.tsv lists only human donors. This table defines \
                 what the measure calls human; a unit from another species belongs in the \
                 antibody dataset's catalogue, not here.",
                row.unit,
                metadata.species
            );
        }
        let donor = antibody::donor_key(&metadata, &row.study);
        if donor != row.donor {
            bail!(
                "{} is donor {donor} but oas_human_donors.tsv says {}. The donor column is what \
                 the disjointness check reads, so it cannot be stale.",
                row.unit,
                row.donor
            );
        }
        // Keyed by donor, not by unit: two units from one person are one
        // repertoire, or that person would count twice in every prevalence.
        repertoires.entry(donor).or_default().extend(sequences.into_iter().map(|s| s.sequence));
    }

    // Donors too shallow to be evidence come out of the denominator and get
    // named in the card. See MIN_SEQUENCES_PER_DONOR for why, and why 500.
    let (shallow, present): (Vec<_>, Vec<_>) =
        repertoires.into_iter().partition(|(_, sequences)| sequences.len() < MIN_SEQUENCES_PER_DONOR);
    let dropped: Vec<(String, usize)> =
        shallow.into_iter().map(|(donor, sequences)| (donor, sequences.len())).collect();
    let repertoires: BTreeMap<String, Vec<String>> = present.into_iter().collect();
    for (donor, kept) in &dropped {
        eprintln!(
            "{donor} contributed {kept} sequences, under the {MIN_SEQUENCES_PER_DONOR} floor; \
             not counted as a donor"
        );
    }

    let sequences: u64 = repertoires.values().map(|r| r.len() as u64).sum();
    if sequences == 0 {
        bail!("no sequences survived reading; nothing to build a table from");
    }
    // Below a hundred donors the loosest published threshold stops being
    // expressible at all, which is the whole reason this catalogue is as wide
    // as it is. Losing that many to filtering is a catalogue problem, not
    // something to report in a footnote.
    if repertoires.len() < 100 {
        bail!(
            "only {} of {} donors contributed any sequence, which is too few to express a 1% \
             threshold. Choose donor units whose sequences survive the {}-{} residue band rather \
             than widening the band for this table alone.",
            repertoires.len(),
            repertoires.len() + dropped.len(),
            MIN_RESIDUES,
            MAX_RESIDUES
        );
    }
    eprintln!("{} donors, {sequences} sequences; building {PEPTIDE_LEN}-mer prevalence", repertoires.len());

    // Built with no floor so the whole distribution can be reported, then
    // restricted. The denominator is the donor count either way.
    let mut table = PrevalenceTable::build(repertoires.values(), 1);
    let distinct = table.len();
    let at_stringency: Vec<(Stringency, usize)> = Stringency::ALL
        .iter()
        .map(|&s| {
            let floor = s.fraction() * f64::from(table.donors());
            (s, table.entries().filter(|(_, donors)| f64::from(*donors) >= floor).count())
        })
        .collect();
    table.restrict(MIN_DONORS);
    eprintln!("{distinct} distinct peptides, {} held by at least {MIN_DONORS} donors", table.len());

    let dir = ctx.dataset_dir(NAME);
    fs::create_dir_all(&dir)?;
    let table_path = dir.join("table.bin");
    {
        let file =
            fs::File::create(&table_path).with_context(|| format!("creating {}", table_path.display()))?;
        table.write(BufWriter::with_capacity(1 << 20, file))?;
    }
    let (table_digest, table_bytes) = crate::cache::hash_file(&table_path)?;

    let stats = TableStats {
        donors: table.donors(),
        units: catalogue.len(),
        // Studies behind the donors that are actually counted, not behind the
        // catalogue: the headline says "N donors from M studies", and a study
        // all of whose donors were dropped for depth is not one of them.
        studies: catalogue
            .iter()
            .filter(|r| repertoires.contains_key(&r.donor))
            .map(|r| r.study.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        sequences,
        read,
        distinct_peptides: distinct,
        kept_peptides: table.len(),
        at_stringency,
        dropped,
        depth: (
            repertoires.values().map(|r| r.len() as u64).min().unwrap_or(0),
            repertoires.values().map(|r| r.len() as u64).max().unwrap_or(0),
        ),
        unpinned,
        bytes: table_bytes,
    };

    fs::write(dir.join("README.md"), render_card(&stats, &catalogue_digest, &table_digest))?;
    let manifest = Manifest {
        dataset: NAME,
        pipeline: "opendna-data",
        pipeline_version: env!("CARGO_PKG_VERSION"),
        table_format: FORMAT,
        peptide_length: PEPTIDE_LEN,
        donors: stats.donors,
        min_donors: MIN_DONORS,
        distinct_peptides: stats.distinct_peptides,
        peptides_written: stats.kept_peptides,
        sequences,
        catalogue_sha256: catalogue_digest,
        table_sha256: table_digest,
        table_bytes,
        units: digests,
        donors_below_depth_floor: &stats.dropped,
        min_sequences_per_donor: MIN_SEQUENCES_PER_DONOR,
    };
    fs::write(dir.join("manifest.json"), serde_json::to_string_pretty(&manifest)? + "\n")?;

    println!("{NAME}: wrote {}", dir.display());
    println!("  {:<22} {:>12}", "donors", stats.donors);
    println!("  {:<22} {:>12}", "sequences", stats.sequences);
    println!("  {:<22} {:>12}", "distinct peptides", stats.distinct_peptides);
    println!("  {:<22} {:>12}", "peptides written", stats.kept_peptides);
    println!("  {:<22} {:>12}", "table bytes", stats.bytes);
    for (stringency, count) in &stats.at_stringency {
        println!("  {:<22} {count:>12}", format!("at {}", stringency.label()));
    }
    if stats.unpinned > 0 {
        println!("  {:<22} {:>12}", "units not pinned", stats.unpinned);
    }
    Ok(())
}

/// The dataset card.
///
/// Written here rather than through [`crate::card`] because that renderer
/// describes Parquet splits and a column schema, and this output has neither.
/// The sections a reader needs are different, not fewer.
fn render_card(stats: &TableStats, catalogue_digest: &str, table_digest: &str) -> String {
    use std::fmt::Write as _;
    let mut s = String::new();
    let smallest = 1.0 / f64::from(stats.donors.max(1)) * 100.0;

    s.push_str("---\nlicense: cc-by-4.0\n");
    s.push_str("pretty_name: OpenDNA peptide humanness prevalence table\n");
    s.push_str("tags:\n- protein\n- antibody\n- humanness\n- immunogenicity\n");
    s.push_str("---\n\n");

    s.push_str("# OpenDNA peptide humanness prevalence table\n\n");
    let _ = writeln!(
        s,
        "How many human donors hold each {PEPTIDE_LEN}-residue peptide, over {} donors from {} \
         studies of the Observed Antibody Space. Cut an antibody sequence into overlapping \
         {PEPTIDE_LEN}-mers, look each one up, and every window comes back with the fraction of \
         donors whose repertoire contains it. The per-window list is the useful output: it says \
         *where* a sequence is unlike anything human, which a single score cannot.\n",
        stats.donors, stats.studies
    );

    s.push_str("\n## Read this before using the numbers\n\n");
    let _ = writeln!(
        s,
        "Each donor contributes one shallow data unit: between {} and {} sequences survive \
         reading, out of a repertoire of millions. A peptide that donor really has, but that did \
         not appear in the sample, is counted as absent. **Every prevalence here is biased \
         downwards**, and unevenly: the deepest donor contributes {:.1} times as many sequences \
         as the shallowest, so its peptides are the more likely to be seen at all.\n\n\
         Scores from this table are comparable with each other. They are **not** comparable with \
         published OASis values: that work used deep repertoires from 266 subjects, so the same \
         sequence will score lower here. The choice was deliberate -- prevalence is a fraction of \
         donors, so donor breadth decides what thresholds mean at all -- but it is a trade, not a \
         free lunch.",
        stats.depth.0,
        stats.depth.1,
        if stats.depth.0 == 0 { 0.0 } else { stats.depth.1 as f64 / stats.depth.0 as f64 }
    );

    s.push_str("\n## Sources and licenses\n\n| Source | License | Index | SHA-256 of unit list |\n|---|---|---|---|\n");
    let src = sources::OAS;
    // Not `src.description`: that one names `oas_units.tsv`, which is the
    // antibody dataset's catalogue, and this table is built from a different
    // list of units entirely.
    let _ = writeln!(
        s,
        "| Observed Antibody Space, unpaired human heavy-chain data units (Olsen, Boyles and \
         Deane 2021); the units read are listed in `oas_human_donors.tsv` | [{}]({}) | \
         [OAS]({}) | `{catalogue_digest}` |",
        src.license.name(),
        src.license.url(),
        src.url
    );
    let _ = writeln!(
        s,
        "\nOAS publishes its data units under CC BY 4.0. This table is a derived work -- counts \
         of peptides, not sequences -- and is released under the same license, with attribution \
         to OAS and to the studies whose repertoires it was built from. Every unit's SHA-256 is \
         recorded in `manifest.json`; the exact list, with donor keys, is `oas_human_donors.tsv` \
         in the pipeline source.\n"
    );

    s.push_str("## Construction\n\n");
    let _ = writeln!(
        s,
        "1. One OAS unpaired heavy-chain unit per donor, {} of them, picked round-robin across \
         studies so no single study dominates, from a fixed size band and for breadth of donors \
         rather than depth per donor. {} of those donors are counted, across {} studies; the \
         others contributed too little to be evidence, and are named below.\n\
         2. Rows kept only where OAS marks them productive, stop-codon-free, in frame and \
         heavy-locus, with an amino-acid V(D)J alignment of {}-{} standard residues. \
         {} rows were read and {} kept; the drops are itemised below.\n\
         3. Sequences grouped by donor key, which is the subject identifier where OAS gives a \
         global one and the study otherwise. Two units from one person would be one repertoire, \
         not two donors.\n\
         4. Every overlapping {PEPTIDE_LEN}-mer of every sequence, counted **once per donor** \
         however many of that donor's sequences contain it. Counting sequences instead would let \
         one deeply sequenced donor make their own peptides look universal.\n\
         5. Peptides held by fewer than {MIN_DONORS} donors dropped: {} of {} distinct peptides. \
         One donor in {} is {:.2}%, below the loosest threshold anyone scores at, so this changes \
         no threshold answer -- but it is most of the table by count, and the reported prevalence \
         of such a peptide becomes zero rather than {:.2}%.\n",
        stats.units,
        stats.donors,
        stats.studies,
        antibody::MIN_RESIDUES,
        antibody::MAX_RESIDUES,
        stats.read.rows,
        stats.read.kept,
        stats.distinct_peptides - stats.kept_peptides,
        stats.distinct_peptides,
        stats.donors,
        smallest,
        smallest,
    );

    let _ = writeln!(
        s,
        "No donor in this table appears in the antibody dataset's catalogue (`oas_units.tsv`), \
         and the build fails if one does. The two donor sets have to be disjoint or an \
         evaluation measures how well the table remembers its own input."
    );

    if !stats.dropped.is_empty() {
        let _ = writeln!(
            s,
            "\n{} donor(s) listed in the catalogue contributed fewer than {} sequences and are \
             **not** among the {} counted: {}. A donor who contributed almost nothing votes \
             \"absent\" on almost every peptide, so counting them would divide every prevalence \
             by donors who could not have held it. The floor is not tuned: depth per donor here \
             is bimodal, nine donors fall below 500 and the next has 1,153, so any floor between \
             those two gives the same set.",
            stats.dropped.len(),
            MIN_SEQUENCES_PER_DONOR,
            stats.donors,
            stats
                .dropped
                .iter()
                .map(|(donor, kept)| format!("{donor} ({kept})"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    s.push_str("\n### Rows read and why they were dropped\n\n| Reason | Rows |\n|---|---|\n");
    let r = &stats.read;
    for (reason, n) in [
        ("read", r.rows),
        ("not productive", r.not_productive),
        ("stop codon", r.stop_codon),
        ("out of frame", r.out_of_frame),
        ("wrong locus", r.wrong_locus),
        ("empty sequence", r.empty_sequence),
        ("non-standard residue", r.nonstandard_residue),
        ("outside the length band", r.out_of_length_band),
        ("duplicate within the unit", r.duplicate_within_unit),
        ("**kept**", r.kept),
    ] {
        let _ = writeln!(s, "| {reason} | {n} |");
    }

    s.push_str("\n## Files\n\n");
    let _ = writeln!(
        s,
        "`table.bin` -- {} bytes, SHA-256 `{table_digest}`. Format {FORMAT}: a 16-byte magic, \
         then format version, peptide length, donor total and entry count as little-endian \
         integers, then one entry per peptide of {PEPTIDE_LEN} ASCII residues and a \
         little-endian `u32` donor count, sorted by peptide. Readers refuse a file whose entries \
         are unsorted or whose counts exceed the donor total, because either would produce wrong \
         scores silently rather than failing.\n\n\
         `manifest.json` -- the SHA-256 of every unit read, the catalogue digest, and the counts \
         above.",
        stats.bytes
    );

    s.push_str("\n\n## What the thresholds can mean\n\n");
    let _ = writeln!(
        s,
        "From {} donors a prevalence takes one of {} values, so the smallest non-zero prevalence \
         expressible is {:.2}%. A threshold finer than that cannot be evaluated at all, however \
         it is reported.\n",
        stats.donors,
        u64::from(stats.donors) + 1,
        smallest
    );
    s.push_str("| Stringency | Expressible here | Peptides at or above |\n|---|---|---|\n");
    for (stringency, count) in &stats.at_stringency {
        let _ = writeln!(
            s,
            "| {} | {} | {count} |",
            stringency.label(),
            if stringency.expressible_with(stats.donors) { "yes" } else { "**no**" }
        );
    }

    s.push_str("\n## Limitations\n\n");
    for limitation in LIMITATIONS {
        let _ = writeln!(s, "- {limitation}");
    }

    let _ = writeln!(s, "\n## Citations\n\n{}", CITATIONS.trim());
    s
}

/// What this table cannot tell you. Written to be read before the numbers are
/// used, not as a disclaimer after them.
pub const LIMITATIONS: &[&str] = &[
    "Prevalence is biased downwards by shallow sampling, and the bias is uneven across donors. \
     Absolute values are not comparable with OASis or with any deep-repertoire table.",
    "Heavy chains only. Nothing here says anything about a light chain, and a paired antibody \
     scored window-by-window on its heavy chain alone is half measured.",
    "The donors are whoever OAS has, which is not a sample of humanity: the studies are \
     overwhelmingly from North America and Europe, and several are patient cohorts rather than \
     healthy adults. A peptide common in an unrepresented population reads as rare here.",
    "A high prevalence means a peptide occurs in human repertoires, not that an antibody \
     carrying it is safe, non-immunogenic or developable. Immunogenicity depends on MHC \
     presentation, dose, route and the recipient, none of which this measures.",
    "A low prevalence means the peptide was not seen in these donors' sampled sequences. It is \
     not evidence that no human has it.",
    "The measure is a count over a fixed table, not a model. It will not generalise to a \
     peptide length other than nine or to a different table without being rebuilt, and nothing \
     about it is learned, so there is no accuracy figure to report.",
];

pub const CITATIONS: &str = r"The data are OAS; the measure follows OASis, published as
part of BioPhi. This table is built from different, shallower repertoires than
OASis and its values are not interchangeable with those published there.

```bibtex
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
}

@article{prihoda2022biophi,
  title   = {BioPhi: A platform for antibody design, humanization, and
             humanness evaluation based on natural antibody repertoires and
             deep learning},
  author  = {Prihoda, David and Maamary, Jad and Waight, Andrew and Juan,
             Veronica and Fayadat-Dilman, Laurence and Svozil, Daniel and
             Bitton, Danny A.},
  journal = {mAbs},
  volume  = {14},
  number  = {1},
  pages   = {2020203},
  year    = {2022},
  doi     = {10.1080/19420862.2021.2020203}
}
```";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_committed_catalogue_parses_and_is_wide_enough_to_mean_something() {
        let rows = read_donors(DONORS).expect("committed catalogue parses");
        // A table narrower than 100 donors cannot express a 1% threshold, which
        // is the whole reason this catalogue exists separately.
        assert!(rows.len() >= 100, "{} donors", rows.len());
        assert!(Stringency::Loose.expressible_with(u32::try_from(rows.len()).unwrap()));

        let donors: BTreeSet<&str> = rows.iter().map(|r| r.donor.as_str()).collect();
        assert_eq!(donors.len(), rows.len(), "a donor is listed twice");
        let studies: BTreeSet<&str> = rows.iter().map(|r| r.study.as_str()).collect();
        // Breadth of studies, not just of donors: one study's cohort shares
        // ancestry, exposure and library preparation.
        assert!(studies.len() >= 20, "{} studies", studies.len());

        for row in &rows {
            assert!(row.unit.ends_with(".csv.gz"), "{}", row.unit);
            assert!(!row.unit.starts_with('/') && !row.unit.contains(".."), "{}", row.unit);
            assert!(row.url().starts_with(antibody::OAS_UNIT_BASE), "{}", row.url());
            assert!(!row.file_name().contains('/'), "{}", row.file_name());
            assert!(row.donor.starts_with("person:") || row.donor.starts_with("study:"), "{}", row.donor);
        }
    }

    #[test]
    fn the_table_donors_and_the_evaluation_donors_are_disjoint() {
        let table: BTreeSet<String> = read_donors(DONORS).unwrap().into_iter().map(|r| r.donor).collect();
        let evaluation = evaluation_donors().unwrap();
        let shared: Vec<&String> = table.intersection(&evaluation).collect();
        assert!(shared.is_empty(), "donors in both catalogues: {shared:?}");
    }

    #[test]
    fn a_catalogue_with_a_traversing_unit_path_is_refused() {
        let text = "study\tdonor\tisotype\tunique_sequences\tunit\tsha256\n\
                    S\tperson:A\tIGHM\t100\t../../etc/passwd\t-\n";
        assert!(read_donors(text).is_err());
    }

    #[test]
    fn an_empty_catalogue_is_an_error_rather_than_an_empty_table() {
        let text = "study\tdonor\tisotype\tunique_sequences\tunit\tsha256\n";
        assert!(read_donors(text).is_err());
    }
}
