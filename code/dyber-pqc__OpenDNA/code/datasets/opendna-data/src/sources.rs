//! Catalogue of upstream sources, their licenses, and sources that were
//! evaluated and excluded.
//!
//! URLs point at versioned artefacts wherever the upstream offers one (a
//! Zenodo record, a named DisProt release, a CATH release directory, a git
//! commit). Where it does not (UniProt REST, the wwPDB index), the lockfile
//! checksum is what makes a build reproducible: a changed upstream file is
//! detected rather than silently used.

/// License of an upstream source, as verified on the source's own site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceLicense {
    Cc0,
    CcBy4,
    CcBySa21Jp,
    CcBy3,
}

impl SourceLicense {
    pub fn name(self) -> &'static str {
        match self {
            SourceLicense::Cc0 => "CC0 1.0",
            SourceLicense::CcBy4 => "CC BY 4.0",
            SourceLicense::CcBySa21Jp => "CC BY-SA 2.1 JP",
            SourceLicense::CcBy3 => "CC BY 3.0",
        }
    }

    pub fn url(self) -> &'static str {
        match self {
            SourceLicense::Cc0 => "https://creativecommons.org/publicdomain/zero/1.0/",
            SourceLicense::CcBy4 => "https://creativecommons.org/licenses/by/4.0/",
            SourceLicense::CcBySa21Jp => "https://creativecommons.org/licenses/by-sa/2.1/jp/deed.en",
            SourceLicense::CcBy3 => "https://creativecommons.org/licenses/by/3.0/",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Source {
    /// Stable id; also the cache subdirectory and the lockfile key.
    pub id: &'static str,
    pub url: &'static str,
    pub file_name: &'static str,
    pub license: SourceLicense,
    pub description: &'static str,
}

/// A source that was requested or considered and deliberately left out.
#[derive(Debug, Clone, Copy)]
pub struct ExcludedSource {
    pub name: &'static str,
    pub dataset: &'static str,
    pub reason: &'static str,
}

pub const TSUBOYAMA_2023: Source = Source {
    id: "tsuboyama2023",
    url: "https://zenodo.org/api/records/7992926/files/Processed_K50_dG_datasets.zip/content",
    file_name: "Processed_K50_dG_datasets.zip",
    license: SourceLicense::CcBy4,
    description: "Tsuboyama et al. 2023 mega-scale cDNA display proteolysis, processed dG/ddG tables (Zenodo 10.5281/zenodo.7992926)",
};

pub const DISPROT_2025_12: Source = Source {
    id: "disprot-2025_12",
    url: "https://disprot.org/api/search?release=2025_12&show_ambiguous=false&show_obsolete=false&format=json&page_size=100000",
    file_name: "disprot_2025_12.json",
    license: SourceLicense::CcBy4,
    description: "DisProt release 2025_12, all entries, ambiguous and obsolete annotations excluded",
};

/// CAID3 Disorder-PDB reference, from the official CAID assessment repository,
/// pinned to a commit.
pub const CAID3_DISORDER_PDB: Source = Source {
    id: "caid3-disorder-pdb",
    url: "https://raw.githubusercontent.com/BioComputingUP/CAID/7bab7e8880d7f949e480fdb377c5a5c8546ae336/demo-data/references/disorder_pdb.fasta",
    file_name: "disorder_pdb.fasta",
    license: SourceLicense::CcBy3,
    description: "CAID3 Disorder-PDB reference set (BioComputingUP/CAID @ 7bab7e8); only the protein ids and sequences are used, to hold them out",
};

pub const ESOL_NBDC: Source = Source {
    id: "esol-nbdc",
    url: "https://dbarchive.biosciencedbc.jp/data/esol/LATEST/esol.zip",
    file_name: "esol.zip",
    license: SourceLicense::CcBySa21Jp,
    description:
        "eSOL (Niwa et al. 2009) as archived by the NBDC LSDB Archive, doi:10.18908/lsdba.nbdc00440-000",
};

pub const UNIPROT_ECOLI_K12: Source = Source {
    id: "uniprot-up000000625",
    url: "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000625&fields=accession,reviewed,gene_oln,sequence&format=tsv",
    file_name: "up000000625.tsv",
    license: SourceLicense::CcBy4,
    description: "UniProtKB entries of the E. coli K-12 MG1655 reference proteome (UP000000625): accession, review status, ordered locus names, sequence",
};

pub const CATH_DOMAIN_LIST: Source = Source {
    id: "cath-4_4_0-domain-list-s100",
    url: "https://download.cathdb.info/cath/releases/all-releases/v4_4_0/cath-classification-data/cath-domain-list-S100-v4_4_0.txt",
    file_name: "cath-domain-list-S100-v4_4_0.txt",
    license: SourceLicense::CcBy4,
    description: "CATH v4.4.0 domain classification, S100 non-redundant representatives",
};

pub const CATH_DOMAIN_SEQS: Source = Source {
    id: "cath-4_4_0-domain-seqs-s100",
    url: "https://download.cathdb.info/cath/releases/all-releases/v4_4_0/sequence-data/cath-domain-seqs-S100-v4_4_0.fa",
    file_name: "cath-domain-seqs-S100-v4_4_0.fa",
    license: SourceLicense::CcBy4,
    description: "CATH v4.4.0 domain sequences (from ATOM records), S100 representatives",
};

pub const PDB_ENTRIES_IDX: Source = Source {
    id: "wwpdb-entries-idx",
    url: "https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx",
    file_name: "entries.idx",
    license: SourceLicense::Cc0,
    description: "wwPDB entries index: deposition date, title, resolution and experimental method per entry",
};

/// Observed Antibody Space.
///
/// Not fetched through the lockfile, and the card says "(not locked)" for it
/// accordingly. OAS publishes no stable URL for its unit index -- the index
/// comes out of a POST to a search form -- so there is nothing upstream to pin
/// the way the CATH release pins the burial dataset's file list. The catalogue
/// is `oas_units.tsv`, committed and compiled into the builder, and it carries
/// a SHA-256 per data unit. The dataset's own manifest records the digest of
/// every unit actually read, and the card records the digest of the catalogue.
pub const OAS: Source = Source {
    id: "oas-unpaired",
    url: "https://opig.stats.ox.ac.uk/webapps/oas/",
    file_name: "oas_units.tsv",
    license: SourceLicense::CcBy4,
    description: "Observed Antibody Space, unpaired heavy-chain data units (Olsen, Boyles and Deane 2021); the units read are listed in oas_units.tsv",
};

pub const ALL_SOURCES: &[Source] = &[
    TSUBOYAMA_2023,
    DISPROT_2025_12,
    CAID3_DISORDER_PDB,
    ESOL_NBDC,
    UNIPROT_ECOLI_K12,
    CATH_DOMAIN_LIST,
    CATH_DOMAIN_SEQS,
    PDB_ENTRIES_IDX,
    OAS,
];

/// SAbDab, considered for the antibody dataset and not used.
pub const SABDAB_EXCLUDED: ExcludedSource = ExcludedSource {
    name: "SAbDab / SAbDab-nano",
    dataset: "antibody",
    reason: "Not used, and so not relied on. SAbDab holds antibody structures and this is a sequence dataset, so it would add nothing a heavy-chain V region needs. Its licence was therefore not reviewed: third-party pages and its own papers read as CC BY 4.0, but the OPIG site states no terms on a page that renders without JavaScript, and this project does not record a licence it has not seen the maintainers state. If a structural antibody dataset is built later, that review has to happen first.",
};

pub const EXCLUDED_SOURCES: &[ExcludedSource] = &[
    SABDAB_EXCLUDED,
    ExcludedSource {
        name: "FireProtDB 2.0",
        dataset: "stability",
        reason: "The FireProtDB terms of use (https://loschmidt.chemi.muni.cz/fireprotdb/terms-of-use/, checked 2026-11-14) grant a non-exclusive, non-transferable license and prohibit sublicensing or transferring rights to the service output and creating derivative works. That does not permit redistribution on Hugging Face. Third-party mirrors label it CC BY 4.0, but that is not the maintainers' statement. FireProtDB 2.0 also already incorporates the Tsuboyama 2023 measurements (its MEGASCALE columns), so including both would double-count them.",
    },
    ExcludedSource {
        name: "NESG solubility scores (Price et al. 2011)",
        dataset: "solubility",
        reason: "The article (doi:10.1186/2042-5783-1-6) is CC BY, but its only supplementary file is a methods document with no per-protein data. The per-target expression and solubility scores used by later work were obtained from the authors on request and have no public release with a stated license.",
    },
    ExcludedSource {
        name: "PSI:Biology / TargetTrack",
        dataset: "solubility",
        reason: "The PSI Structural Biology Knowledgebase and TargetTrack are retired. No maintained public download with a stated redistribution license was found.",
    },
];

pub fn by_id(id: &str) -> Option<&'static Source> {
    ALL_SOURCES.iter().find(|s| s.id == id)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn ids_and_cache_paths_are_unique() {
        let ids: BTreeSet<_> = ALL_SOURCES.iter().map(|s| s.id).collect();
        assert_eq!(ids.len(), ALL_SOURCES.len());
        for s in ALL_SOURCES {
            assert!(s.url.starts_with("https://"), "{}", s.id);
            assert!(!s.file_name.contains('/'), "{}", s.id);
            assert!(by_id(s.id).is_some());
        }
    }
}
