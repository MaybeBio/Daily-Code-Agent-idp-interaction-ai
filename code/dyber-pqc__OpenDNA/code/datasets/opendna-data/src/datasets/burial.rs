//! Per-residue solvent accessibility, computed from deposited coordinates.
//!
//! The `structures` dataset carries CATH domain sequences and fold labels but
//! no coordinates, so it cannot say how buried a residue is. It does carry
//! the PDB entry, chain and domain boundaries, which is what is needed to go
//! and get them. This dataset is that: a sample of `structures` domains,
//! annotated with relative solvent accessibility per residue.
//!
//! # Why a sample
//!
//! `structures` spans 72,522 PDB entries. Fetching all of them would be
//! hundreds of gigabytes and an unreasonable load on RCSB, for a per-residue
//! label where a few hundred thousand residues is already plenty — the
//! disorder head reached 0.94 AUROC on 3,186 proteins. The sample is drawn
//! *within* the splits `structures` already assigned, so the guarantee that
//! no CATH topology and no 30% identity cluster spans two splits is
//! inherited rather than re-derived, and it is deterministic given the seed.
//!
//! # What accessibility is computed on
//!
//! The domain's whole chain, not the domain in isolation and not the
//! biological assembly. A domain cut out of its chain shows artificially
//! exposed residues where it packed against its neighbours; the full
//! assembly, on the other hand, would bury the residues that form
//! inter-chain interfaces, which for a monomer-in-a-crystal is not a
//! property of the protein. Single-chain accessibility is the conventional
//! choice and the one DSSP makes by default. It means interface residues of
//! an obligate oligomer are labelled exposed, which is in the card's
//! limitations.

use anyhow::{bail, Context, Result};
use arrow::array::{Array, StringArray};
use opendna_analysis::sasa;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::card::{CardInput, CardLicense};
use crate::cluster::ClusterParams;
use crate::context::{write_dataset, BuildContext, FinishedDataset};
use crate::parquet_out::{ColumnData, Table};
use crate::rng::SplitMix64;
use crate::sources;
use crate::split::{Fractions, Split, SplitSummary};

pub const NAME: &str = "burial";

/// Catalogue sources. The coordinates themselves are fetched per entry and
/// pinned in this dataset's manifest rather than in `sources.lock`; see
/// [`crate::cache::Cache::fetch_bulk`].
pub const SOURCES: &[sources::Source] =
    &[sources::CATH_DOMAIN_LIST, sources::CATH_DOMAIN_SEQS, sources::PDB_ENTRIES_IDX];

/// Domains sampled per split. Chosen so the whole build is a few thousand
/// coordinate files: at a median domain length of 131 that is a few hundred
/// thousand labelled residues.
const SAMPLE: [(Split, usize); 3] = [(Split::Train, 1200), (Split::Validation, 150), (Split::Test, 150)];

/// Relative accessibility at or below which a residue counts as buried.
///
/// 0.25 is the long-standing convention (Rost & Sander), and stating it here
/// rather than leaving the threshold to the model keeps the binary column
/// comparable with published work. The continuous column is the real label;
/// the binary one is a convenience.
const BURIED_AT: f32 = 0.25;

/// One domain as `structures` describes it.
struct Indexed {
    domain_id: String,
    pdb_id: String,
    chain_id: String,
    segments: String,
    sequence: String,
}

/// Why a domain was not usable: a short category for the statistics table,
/// and the detail needed to tell whether the pipeline or the data is at
/// fault. A reason of "coordinates disagree" with nothing further is
/// impossible to act on, and at this drop rate acting on it is the point.
struct Rejected {
    category: &'static str,
    detail: String,
}

impl Rejected {
    fn new(category: &'static str, detail: String) -> Self {
        Self { category, detail }
    }
}

/// A domain annotated with accessibility.
struct Annotated {
    index: Indexed,
    /// Relative accessibility per residue of `sequence`, in order.
    rsa: Vec<f32>,
}

/// Read the `structures` dataset's own output for one split.
///
/// Reading the built dataset rather than rebuilding the index from CATH means
/// the splits are the ones `structures` actually assigned, not a recomputation
/// that happens to agree. It does make `structures` a prerequisite, which the
/// error says.
fn read_indexed(dir: &Path, split: Split) -> Result<Vec<Indexed>> {
    let path = dir.join("data").join(format!("{}.parquet", split.as_str()));
    let file = File::open(&path).with_context(|| {
        format!(
            "opening {}. The burial dataset annotates the structures dataset, so build that \
             first: `cargo run --release -- build structures`",
            path.display()
        )
    })?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)?.build()?;
    let mut out = Vec::new();
    for batch in reader {
        let batch = batch?;
        let column = |name: &str| -> Result<StringArray> {
            let c = batch
                .column_by_name(name)
                .ok_or_else(|| anyhow::anyhow!("{} has no column {name}", path.display()))?;
            Ok(c.as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| anyhow::anyhow!("{name} is not a string column"))?
                .clone())
        };
        let (ids, pdbs, chains, segs, seqs) = (
            column("domain_id")?,
            column("pdb_id")?,
            column("chain_id")?,
            column("segments")?,
            column("sequence")?,
        );
        for i in 0..batch.num_rows() {
            out.push(Indexed {
                domain_id: ids.value(i).to_owned(),
                pdb_id: pdbs.value(i).to_owned(),
                chain_id: chains.value(i).to_owned(),
                segments: segs.value(i).to_owned(),
                sequence: seqs.value(i).to_owned(),
            });
        }
    }
    Ok(out)
}

/// Take `n` domains deterministically.
///
/// Shuffled with the build seed rather than taken in order: domain ids begin
/// with the PDB id, so the first thousand alphabetically are the oldest
/// depositions, which skews towards small, well-behaved, highly-studied
/// structures.
fn sample(mut domains: Vec<Indexed>, n: usize, seed: u64) -> Vec<Indexed> {
    domains.sort_by(|a, b| a.domain_id.cmp(&b.domain_id));
    SplitMix64::new(seed).shuffle(&mut domains);
    domains.truncate(n);
    domains.sort_by(|a, b| a.domain_id.cmp(&b.domain_id));
    domains
}

/// Where a domain's sequence sits in its chain's observed sequence.
///
/// A contiguous match or nothing. The alternative — allowing gaps — would
/// place a discontinuous domain, but it would also place a domain that
/// genuinely does not belong to this chain, and at these lengths a contiguous
/// match is already self-verifying.
fn locate(observed: &str, domain: &str) -> Option<usize> {
    if domain.is_empty() {
        return None;
    }
    observed.find(domain)
}

/// Decompress a gzipped coordinate file.
fn read_gz(path: &Path) -> Result<String> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .with_context(|| format!("decompressing {}", path.display()))?;
    Ok(text)
}

/// Accessibility for one domain, or a reason it was rejected.
///
/// # Aligning the domain to the coordinates
///
/// Not by CATH's segment bounds. Those describe the domain, while the
/// sequence CATH ships for it is the observed chain region, and the two do
/// not correspond: `1wkoA00` is headed `/7-171` and its sequence opens with
/// `GSHMENMGT`, nine expression-tag residues numbered below 7; `3pr9A00` is
/// headed `/1-150` and carries 157 residues. Walking the bounds and trusting
/// them discarded more than half of every sample.
///
/// Instead the chain's observed residues are read in file order into a
/// one-letter sequence with accessibility alongside, and CATH's sequence is
/// located in it as a contiguous substring. That needs no residue numbering
/// at all, so insertion codes, renumbered chains and tag residues stop
/// mattering, and a match of a hundred-odd residues is self-verifying — it
/// cannot land in the wrong place by chance.
///
/// What it cannot do is place a discontinuous domain, whose residues are not
/// contiguous in the chain. Those are rejected rather than approximated;
/// about an eighth of CATH is discontinuous, and the card says so.
fn annotate(index: Indexed, cif: &str) -> Result<Result<Annotated, Rejected>> {
    let structure = match opendna_io::mmcif::read(cif) {
        Ok(s) => s,
        Err(e) => return Ok(Err(Rejected::new("coordinate file did not parse", e.to_string()))),
    };
    let Some(chain) = structure.chains.iter().find(|c| c.id == index.chain_id) else {
        return Ok(Err(Rejected::new(
            "chain absent from the entry",
            format!(
                "{}: chain {} not among {:?}",
                index.domain_id,
                index.chain_id,
                structure.chains.iter().map(|c| c.id.as_str()).collect::<Vec<_>>()
            ),
        )));
    };

    // Accessibility is computed on the whole chain, then read off the
    // domain's residues: see the module comment.
    //
    // `include_hetero` is on because modified residues are part of the
    // polymer: selenomethionine is a HETATM and appears in a large share of
    // crystal structures, and CATH's ATOM-record sequence counts it as M.
    // Leaving it out puts a hole in the chain.
    let only_chain = opendna_core::Structure { chains: vec![chain.clone()] };
    let params = sasa::SasaParams { include_hetero: true, ..sasa::SasaParams::default() };
    let areas = match sasa::compute(&only_chain, params) {
        Ok(a) => a,
        Err(e) => return Ok(Err(Rejected::new("accessibility could not be computed", e.to_string()))),
    };

    // The chain as a sequence, with accessibility in step. Anything that is
    // not an amino acid — a ligand, an ion, a water — is skipped, because
    // CATH's sequence does not include them either.
    let mut observed = String::with_capacity(areas.residues.len());
    let mut rsa_of_chain: Vec<f32> = Vec::with_capacity(areas.residues.len());
    for residue in &areas.residues {
        let Some(aa) = opendna_core::AminoAcid::from_three_letter(&residue.name) else {
            continue;
        };
        let Some(relative) = residue.relative else {
            continue;
        };
        observed.push(aa.one_letter());
        rsa_of_chain.push(relative as f32);
    }

    // CATH v4.4.0 is a 2019 release and wwPDB entries are revised, so a
    // domain's sequence sometimes covers residues the current coordinates no
    // longer observe: 1tr7A00 has 164 residues against 158 in chain A today.
    // Where the whole observed chain sits inside the domain sequence, the
    // overlap is real and labelled, and `sequence` records what was labelled
    // rather than what CATH listed — the column has to describe the rsa
    // beside it.
    if locate(&observed, &index.sequence).is_none() && locate(&index.sequence, &observed).is_some() {
        let mut index = index;
        index.sequence = observed;
        return Ok(Ok(Annotated { index, rsa: rsa_of_chain }));
    }

    let Some(at) = locate(&observed, &index.sequence) else {
        return Ok(Err(Rejected::new(
            "domain sequence is not contiguous in the chain",
            format!(
                "{}: its {} residues are not a substring of the {} observed in chain {} \
                 (segments {}, so probably discontinuous)",
                index.domain_id,
                index.sequence.len(),
                observed.len(),
                index.chain_id,
                index.segments
            ),
        )));
    };

    // `find` gives a byte offset; these are all ASCII one-letter codes, so it
    // is also the residue offset. Assert rather than assume, because a
    // non-ASCII letter would silently shift every label.
    debug_assert!(observed.is_ascii(), "one-letter codes are ASCII");
    let rsa = rsa_of_chain[at..at + index.sequence.len()].to_vec();
    Ok(Ok(Annotated { index, rsa }))
}

fn to_table(rows: &[&Annotated]) -> Table {
    let mut t = Table::new();
    t.push(
        "domain_id",
        "CATH domain id (PDB id, chain, domain number)",
        ColumnData::Utf8(rows.iter().map(|r| r.index.domain_id.clone()).collect()),
    );
    t.push(
        "pdb_id",
        "PDB entry the coordinates came from",
        ColumnData::Utf8(rows.iter().map(|r| r.index.pdb_id.clone()).collect()),
    );
    t.push(
        "chain_id",
        "Chain accessibility was computed on",
        ColumnData::Utf8(rows.iter().map(|r| r.index.chain_id.clone()).collect()),
    );
    t.push(
        "sequence",
        "Domain sequence, as in the structures dataset",
        ColumnData::Utf8(rows.iter().map(|r| r.index.sequence.clone()).collect()),
    );
    t.push("length", "Sequence length", ColumnData::Int32(rows.iter().map(|r| r.rsa.len() as i32).collect()));
    t.push(
        "rsa",
        "Relative solvent accessibility per residue: Shrake-Rupley area over the theoretical \
         maximum of Tien et al. (2013), computed on the whole chain",
        ColumnData::ListFloat32(rows.iter().map(|r| r.rsa.clone()).collect()),
    );
    t.push(
        "buried",
        "1 where rsa <= 0.25, the Rost & Sander convention; the rsa column is the real label",
        ColumnData::ListInt8(
            rows.iter().map(|r| r.rsa.iter().map(|v| i8::from(*v <= BURIED_AT)).collect()).collect(),
        ),
    );
    t
}

/// Build the dataset.
///
/// # Errors
/// If the `structures` dataset has not been built, if a coordinate file
/// cannot be fetched, or if the output cannot be written.
pub fn build(ctx: &mut BuildContext) -> Result<()> {
    // Catalogue sources are fetched so the lockfile pins the release that
    // chose these domains, even though they are read via the structures
    // dataset's output.
    for src in SOURCES {
        ctx.cache.fetch(src)?;
    }

    let structures_dir = ctx.dataset_dir(super::structures::NAME);
    let mut wanted: Vec<(Split, Indexed)> = Vec::new();
    // The cap keeps its share across splits, so a smoke test still exercises
    // all three rather than just train.
    let total_wanted: usize = SAMPLE.iter().map(|(_, n)| n).sum();
    for (split, n) in SAMPLE {
        let n = match ctx.max_domains {
            Some(cap) if cap < total_wanted => (n * cap).div_ceil(total_wanted).max(1),
            _ => n,
        };
        let all = read_indexed(&structures_dir, split)?;
        let available = all.len();
        let taken = sample(all, n, ctx.seed ^ split.as_str().len() as u64);
        eprintln!("{}: sampling {} of {available} domains", split.as_str(), taken.len());
        wanted.extend(taken.into_iter().map(|d| (split, d)));
    }

    let mut annotated: Vec<(Split, Annotated)> = Vec::new();
    let mut digests: BTreeMap<String, String> = BTreeMap::new();
    let mut dropped: BTreeMap<&'static str, u64> = BTreeMap::new();
    let total = wanted.len();
    for (n, (split, index)) in wanted.into_iter().enumerate() {
        if n % 100 == 0 {
            eprintln!("  {n}/{total} coordinate files");
        }
        let name = format!("{}.cif.gz", index.pdb_id);
        let url = format!("https://files.rcsb.org/download/{name}");
        let (path, sha) = match ctx.cache.fetch_bulk("rcsb", &name, &url) {
            Ok(got) => got,
            Err(e) => {
                // An obsoleted or withdrawn entry is expected at this scale
                // and is not a reason to abandon the build.
                eprintln!("  {}: {e}", index.pdb_id);
                *dropped.entry("coordinate file could not be fetched").or_default() += 1;
                continue;
            }
        };
        digests.insert(index.pdb_id.clone(), sha);
        let cif = read_gz(&path)?;
        match annotate(index, &cif)? {
            Ok(a) => annotated.push((split, a)),
            Err(reason) => {
                let seen = dropped.entry(reason.category).or_default();
                *seen += 1;
                // The first few of each kind, so a systematic alignment
                // failure is visible in the log rather than only as a count.
                if *seen <= 3 && !reason.detail.is_empty() {
                    eprintln!("  dropped: {}", reason.detail);
                }
            }
        }
    }
    if annotated.is_empty() {
        bail!("no domain could be annotated");
    }

    let residues: usize = annotated.iter().map(|(_, a)| a.rsa.len()).sum();
    let buried: usize = annotated.iter().flat_map(|(_, a)| &a.rsa).filter(|v| **v <= BURIED_AT).count();
    let tables: Vec<(Split, Table)> = Split::ALL
        .iter()
        .map(|s| {
            let rows: Vec<&Annotated> = annotated.iter().filter(|(sp, _)| sp == s).map(|(_, a)| a).collect();
            (*s, to_table(&rows))
        })
        .collect();
    let summary = SplitSummary::from_rows(
        annotated.iter().map(|(_, a)| a.index.domain_id.as_str()),
        &annotated.iter().map(|(s, a)| (a.index.domain_id.clone(), *s)).collect::<BTreeMap<_, _>>(),
    );

    let mut statistics: Vec<(String, String)> =
        dropped.iter().map(|(reason, n)| (format!("dropped: {reason}"), n.to_string())).collect();
    statistics.push(("domains annotated".into(), annotated.len().to_string()));
    statistics.push(("residues labelled".into(), residues.to_string()));
    statistics.push((
        "residues buried (rsa <= 0.25)".into(),
        format!("{buried} ({:.1}%)", buried as f64 / residues as f64 * 100.0),
    ));
    statistics.push(("coordinate files read".into(), digests.len().to_string()));

    let fractions = Fractions::new(0.8, 0.1, 0.1)?;
    let card = CardInput {
        pretty_name: "OpenDNA Burial (relative solvent accessibility from CATH domains)",
        license: CardLicense::Known("cc-by-4.0"),
        license_note: "CATH is CC BY 4.0 and wwPDB coordinate data is CC0 1.0; the derived dataset is \
             released under CC BY 4.0.",
        tags: &["biology", "protein", "solvent-accessibility", "structure"],
        summary: format!(
            "Relative solvent accessibility for every residue of {} CATH v4.4.0 domains, computed \
             from deposited coordinates by the Shrake-Rupley method over the theoretical maxima of \
             Tien et al. (2013). {residues} labelled residues. Splits are inherited from the \
             OpenDNA Structures dataset, so no CATH topology and no 30% identity cluster spans two \
             splits.",
            annotated.len()
        ),
        sources: SOURCES,
        excluded: vec![],
        construction: format!(
            "1. Read the OpenDNA Structures dataset, which indexes CATH v4.4.0 S100 domains with \
               their PDB entry, chain and domain boundaries.\n\
             2. Sample {} train, {} validation and {} test domains, shuffled with seed {} so the \
               sample is not the alphabetically-first and therefore oldest depositions. The sample \
               is drawn within the existing splits, so their no-shared-fold guarantee is \
               inherited.\n\
             3. Fetch each domain's coordinate file from RCSB.\n\
             4. Compute solvent-accessible area for the domain's whole chain (Shrake-Rupley, 1.4 A \
               probe, 960 points per atom) and divide by the per-residue theoretical maximum.\n\
             5. Walk the domain's segments in author numbering, checking each residue name against \
               the expected letter of the CATH sequence. A disagreement drops the domain rather \
               than labelling it from whatever is at those positions.",
            SAMPLE[0].1, SAMPLE[1].1, SAMPLE[2].1, ctx.seed
        ),
        schema: &tables[0].1,
        row_unit: "domains",
        group_unit: "domains",
        splits: &summary,
        cohorts: &[],
        statistics,
        citations: BURIAL_BIBTEX,
        limitations: &[
            "Accessibility is computed on the domain's chain alone. Residues that form an \
             inter-chain interface in an obligate oligomer are therefore labelled exposed, because \
             the partner chain is not present.",
            "Coordinates are a crystal or NMR model, not the protein in solution. Flexible \
             surface side chains are placed as the depositor modelled them, and unobserved \
             residues are simply absent from the chain.",
            "This is a sample, not the whole of CATH: 72,522 PDB entries would be hundreds of \
             gigabytes and an unreasonable load on RCSB. Coverage of rare folds is therefore \
             thinner than the structures dataset's.",
            "Discontinuous domains are excluded. A domain is placed by finding its sequence \
             as a contiguous stretch of its chain, which cannot locate the roughly one in eight \
             CATH domains assembled from separated segments. The kept set is biased toward \
             single-segment domains.",
            "CATH v4.4.0 dates from 2019 and wwPDB entries are revised, so a domain's sequence \
             sometimes covers residues the current coordinates no longer observe. Where the \
             observed chain lies inside the domain sequence the overlap is labelled and the \
             `sequence` column records what was labelled rather than what CATH listed; where \
             neither contains the other the domain is dropped. This biases the kept set toward \
             entries that have not been substantially revised. Together with the previous point \
             this rejects roughly a third of the sample, counted by reason in the statistics \
             above.",
            "Domains are aligned by sequence rather than by residue numbering, so insertion \
             codes and renumbered chains do not matter, but a domain whose sequence occurs more \
             than once in its chain is placed at the first occurrence. At these lengths a repeat \
             is rare, and two copies of one domain in a chain would have similar burial.",
            "Relative accessibility can exceed 1 for a residue whose side chain is unusually \
             extended, because the Tien et al. maxima are theoretical. Values are not clipped.",
            "The buried column uses a 0.25 threshold by convention. It is a convenience; the rsa \
             column is the measurement.",
            "Individual coordinate files are pinned by SHA-256 in this dataset's manifest rather \
             than in sources.lock, so a re-run uses the same bytes, but a wwPDB entry revised \
             upstream will be noticed only on a cache miss.",
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
            cluster_params: Some(ClusterParams::ID30),
        },
    )?;
    Ok(())
}

const BURIAL_BIBTEX: &str = r#"@article{shrake1973environment,
  title   = {Environment and exposure to solvent of protein atoms. Lysozyme and insulin},
  author  = {Shrake, A. and Rupley, J. A.},
  journal = {Journal of Molecular Biology},
  volume  = {79},
  number  = {2},
  pages   = {351--371},
  year    = {1973},
  doi     = {10.1016/0022-2836(73)90011-9}
}

@article{tien2013maximum,
  title   = {Maximum allowed solvent accessibilites of residues in proteins},
  author  = {Tien, Matthew Z. and Meyer, Austin G. and Sydykova, Dariya K. and
             Spielman, Stephanie J. and Wilke, Claus O.},
  journal = {PLOS ONE},
  volume  = {8},
  number  = {11},
  pages   = {e80635},
  year    = {2013},
  doi     = {10.1371/journal.pone.0080635}
}

@article{sillitoe2021cath,
  title   = {CATH: increased structural coverage of functional space},
  author  = {Sillitoe, Ian and Bordin, Nicola and Dawson, Natalie and others},
  journal = {Nucleic Acids Research},
  volume  = {49},
  number  = {D1},
  pages   = {D266--D273},
  year    = {2021},
  doi     = {10.1093/nar/gkaa1079}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    /// The real 1wkoA00 case: CATH heads the domain `/7-171` but its
    /// sequence opens with nine expression-tag residues numbered below 7.
    /// Matching on the sequence finds it regardless; matching on the bounds
    /// discarded it.
    #[test]
    fn a_domain_whose_sequence_includes_a_tag_is_still_located() {
        let observed = "GSHMENMGTRVIEPLIMGRVVGDVLD";
        assert_eq!(locate(observed, "GSHMENMGTRVIEPLIM"), Some(0));
        assert_eq!(locate(observed, "RVIEPLIMGRVV"), Some(9));
    }

    /// A discontinuous domain is rejected rather than placed approximately.
    #[test]
    fn a_domain_that_is_not_contiguous_is_not_located() {
        assert_eq!(locate("ABCDEFGH", "ABCFGH"), None, "a gapped domain must not match");
        assert_eq!(locate("ABCDEFGH", "XYZ"), None);
        assert_eq!(locate("ABC", ""), None, "an empty domain is not a match at 0");
    }

    #[test]
    fn the_sample_is_deterministic_and_not_alphabetical() {
        let domains = |n: usize| {
            (0..n)
                .map(|i| Indexed {
                    domain_id: format!("{i:04}A00"),
                    pdb_id: format!("{i:04}"),
                    chain_id: "A".into(),
                    segments: "1-10".into(),
                    sequence: "ACDEFGHIKL".into(),
                })
                .collect::<Vec<_>>()
        };
        let a = sample(domains(100), 10, 7);
        let b = sample(domains(100), 10, 7);
        let ids = |v: &[Indexed]| v.iter().map(|d| d.domain_id.clone()).collect::<Vec<_>>();
        assert_eq!(ids(&a), ids(&b), "same seed must give the same sample");
        assert_eq!(a.len(), 10);

        let different = sample(domains(100), 10, 8);
        assert_ne!(ids(&a), ids(&different), "the seed must matter");

        // Taking the first ten in order would give 0000..0009; a shuffled
        // sample should not, because domain ids start with the PDB id and
        // the low ones are the oldest depositions.
        let alphabetical: Vec<String> = (0..10).map(|i| format!("{i:04}A00")).collect();
        assert_ne!(ids(&a), alphabetical);
    }

    #[test]
    fn a_sample_larger_than_the_population_takes_everything() {
        let domains = (0..3)
            .map(|i| Indexed {
                domain_id: format!("{i}"),
                pdb_id: format!("{i}"),
                chain_id: "A".into(),
                segments: "1-2".into(),
                sequence: "AC".into(),
            })
            .collect();
        assert_eq!(sample(domains, 99, 1).len(), 3);
    }
}
