//! File import: detect format, parse, store one protein per polymer chain.

use crate::commands::fail;
use opendna_app_protocol::CommandError;
use opendna_core::Protein;
use opendna_store::{Database, StructureFile, StructureFormat};
use std::path::Path;

/// File extensions offered in the import dialog.
pub const EXTENSIONS: &[&str] = &["pdb", "ent", "cif", "mmcif", "fasta", "fa", "faa"];

/// Largest file accepted, to keep a stray multi-gigabyte file from freezing the app.
const MAX_BYTES: u64 = 256 * 1024 * 1024;

/// Import a file and return the ids of the proteins stored.
pub fn import_path(db: &Database, path: &Path) -> Result<Vec<String>, CommandError> {
    let size = std::fs::metadata(path).map_err(fail)?.len();
    if size > MAX_BYTES {
        return Err(fail(format!("{} is {} MB; the limit is 256 MB", path.display(), size / 1_048_576)));
    }
    let text = std::fs::read_to_string(path).map_err(fail)?;
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("structure").to_string();
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();

    if matches!(ext.as_str(), "fasta" | "fa" | "faa") {
        return import_fasta(db, &text);
    }
    import_structure(db, &stem, text)
}

fn import_fasta(db: &Database, text: &str) -> Result<Vec<String>, CommandError> {
    let proteins = opendna_io::fasta::read(text).map_err(fail)?;
    proteins.iter().map(|p| db.save(p).map_err(fail)).collect()
}

fn import_structure(db: &Database, stem: &str, text: String) -> Result<Vec<String>, CommandError> {
    let structure = opendna_io::read_structure(&text).map_err(fail)?;
    let format =
        if text.trim_start().starts_with("data_") { StructureFormat::MmCif } else { StructureFormat::Pdb };

    let mut ids = Vec::new();
    for chain in &structure.chains {
        let sequence = chain.sequence();
        if sequence.is_empty() {
            continue; // ligand- or water-only chain
        }
        let label = if chain.id.is_empty() { "_" } else { chain.id.as_str() };
        let id = db.save(&Protein::new(format!("{stem}:{label}"), sequence)).map_err(fail)?;
        // Each chain keeps a copy of the source file so the viewer can show
        // its context. Deduplicating shared files is planned for the project store.
        db.save_structure(&id, &StructureFile { chain: chain.id.clone(), format, content: text.clone() })
            .map_err(fail)?;
        ids.push(id);
    }
    if ids.is_empty() {
        return Err(fail("the file contains no amino-acid chains"));
    }
    Ok(ids)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../crates/opendna-analysis/tests/data").join(name)
    }

    #[test]
    fn pdb_import_stores_chain_and_structure() {
        let db = Database::open_in_memory().unwrap();
        let ids = import_path(&db, &fixture("1ubq.pdb")).unwrap();
        assert_eq!(ids.len(), 1);
        let stored = db.get(&ids[0]).unwrap().unwrap();
        assert_eq!(stored.protein.name, "1ubq:A");
        assert_eq!(stored.protein.sequence.len(), 76);
        let file = db.structure(&ids[0]).unwrap().unwrap();
        assert_eq!((file.chain.as_str(), file.format), ("A", StructureFormat::Pdb));
    }

    #[test]
    fn mmcif_is_detected_by_content() {
        let db = Database::open_in_memory().unwrap();
        let ids = import_path(&db, &fixture("1ubq.cif")).unwrap();
        assert_eq!(db.structure(&ids[0]).unwrap().unwrap().format, StructureFormat::MmCif);
    }

    #[test]
    fn unreadable_file_is_an_error_not_a_panic() {
        let db = Database::open_in_memory().unwrap();
        assert!(import_path(&db, Path::new("does/not/exist.pdb")).is_err());
    }
}
