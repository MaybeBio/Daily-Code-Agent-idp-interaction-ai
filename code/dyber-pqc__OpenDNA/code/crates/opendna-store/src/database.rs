use chrono::{DateTime, Utc};
use opendna_core::{CoreError, Protein, Sequence};
use rusqlite::{params, Connection, OptionalExtension};
use std::path::Path;
use thiserror::Error;

/// Errors from [`Database`].
#[derive(Debug, Error)]
pub enum StoreError {
    /// Underlying `SQLite` failure.
    #[error("database error: {0}")]
    Sqlite(#[from] rusqlite::Error),
    /// A stored sequence no longer parses (corrupted or hand-edited file).
    #[error("stored sequence for {id} is invalid: {source}")]
    CorruptSequence {
        /// Row id.
        id: String,
        /// Parse error.
        source: CoreError,
    },
    /// A stored structure format is not recognised.
    #[error("stored structure format for {0} is invalid")]
    CorruptFormat(String),
    /// A stored timestamp is not RFC 3339.
    #[error("stored timestamp for {0} is invalid")]
    CorruptTimestamp(String),
}

/// A protein row with bookkeeping fields.
#[derive(Debug, Clone, PartialEq)]
pub struct StoredProtein {
    /// Content id (see [`Protein::id`]).
    pub id: String,
    /// The protein. Structures are not stored in this table.
    pub protein: Protein,
    /// First time this sequence was saved.
    pub created_at: DateTime<Utc>,
}

/// Coordinate file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructureFormat {
    /// wwPDB fixed-column format.
    Pdb,
    /// `PDBx`/mmCIF.
    MmCif,
}

impl StructureFormat {
    /// Lower-case name used in storage and by viewers (`pdb`, `mmcif`).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pdb => "pdb",
            Self::MmCif => "mmcif",
        }
    }

    /// Inverse of [`StructureFormat::as_str`].
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "pdb" => Some(Self::Pdb),
            "mmcif" => Some(Self::MmCif),
            _ => None,
        }
    }
}

/// An original coordinate file kept alongside a protein.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureFile {
    /// Chain the protein's sequence was taken from.
    pub chain: String,
    /// File format.
    pub format: StructureFormat,
    /// File text, unmodified.
    pub content: String,
}

/// SQLite-backed protein store.
pub struct Database {
    conn: Connection,
}

const SCHEMA_VERSION: i32 = 2;

impl Database {
    /// Open (creating if needed) a database file.
    ///
    /// # Errors
    /// [`StoreError::Sqlite`] if the file cannot be opened or migrated.
    pub fn open(path: &Path) -> Result<Self, StoreError> {
        Self::init(Connection::open(path)?)
    }

    /// An in-memory database, for tests and scratch work.
    ///
    /// # Errors
    /// [`StoreError::Sqlite`] on failure.
    pub fn open_in_memory() -> Result<Self, StoreError> {
        Self::init(Connection::open_in_memory()?)
    }

    fn init(conn: Connection) -> Result<Self, StoreError> {
        conn.execute_batch(
            "PRAGMA foreign_keys = ON;
             PRAGMA journal_mode = WAL;
             CREATE TABLE IF NOT EXISTS proteins (
                 id         TEXT PRIMARY KEY,
                 name       TEXT NOT NULL,
                 sequence   TEXT NOT NULL,
                 created_at TEXT NOT NULL
             );
             CREATE TABLE IF NOT EXISTS structures (
                 protein_id TEXT PRIMARY KEY REFERENCES proteins(id) ON DELETE CASCADE,
                 chain      TEXT NOT NULL,
                 format     TEXT NOT NULL,
                 content    TEXT NOT NULL
             );",
        )?;
        conn.pragma_update(None, "user_version", SCHEMA_VERSION)?;
        Ok(Self { conn })
    }

    /// Insert a protein, or rename it if the same sequence is already stored.
    /// Returns its id.
    ///
    /// # Errors
    /// [`StoreError::Sqlite`] on failure.
    pub fn save(&self, protein: &Protein) -> Result<String, StoreError> {
        let id = protein.id();
        self.conn.execute(
            "INSERT INTO proteins (id, name, sequence, created_at) VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT(id) DO UPDATE SET name = excluded.name",
            params![id, protein.name, protein.sequence.to_string(), Utc::now().to_rfc3339()],
        )?;
        Ok(id)
    }

    /// Fetch by id.
    ///
    /// # Errors
    /// [`StoreError`] on `SQLite` failure or a corrupted row.
    pub fn get(&self, id: &str) -> Result<Option<StoredProtein>, StoreError> {
        self.conn
            .query_row(
                "SELECT id, name, sequence, created_at FROM proteins WHERE id = ?1",
                params![id],
                |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?)),
            )
            .optional()?
            .map(Self::decode)
            .transpose()
    }

    /// All proteins, newest first.
    ///
    /// # Errors
    /// [`StoreError`] on `SQLite` failure or a corrupted row.
    pub fn list(&self) -> Result<Vec<StoredProtein>, StoreError> {
        let mut stmt = self
            .conn
            .prepare("SELECT id, name, sequence, created_at FROM proteins ORDER BY created_at DESC, id")?;
        let rows = stmt.query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?)))?;
        rows.map(|r| Self::decode(r?)).collect()
    }

    /// Attach the source coordinate file to a stored protein, replacing any
    /// previous one. `chain` names the chain the protein's sequence came from.
    ///
    /// # Errors
    /// [`StoreError::Sqlite`] on failure, including when `protein_id` is unknown.
    pub fn save_structure(&self, protein_id: &str, file: &StructureFile) -> Result<(), StoreError> {
        self.conn.execute(
            "INSERT INTO structures (protein_id, chain, format, content) VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT(protein_id) DO UPDATE SET
                 chain = excluded.chain, format = excluded.format, content = excluded.content",
            params![protein_id, file.chain, file.format.as_str(), file.content],
        )?;
        Ok(())
    }

    /// The coordinate file attached to a protein, if any.
    ///
    /// # Errors
    /// [`StoreError`] on `SQLite` failure or an unknown stored format.
    pub fn structure(&self, protein_id: &str) -> Result<Option<StructureFile>, StoreError> {
        let row: Option<(String, String, String)> = self
            .conn
            .query_row(
                "SELECT chain, format, content FROM structures WHERE protein_id = ?1",
                params![protein_id],
                |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
            )
            .optional()?;
        row.map(|(chain, format, content)| {
            let format = StructureFormat::parse(&format)
                .ok_or_else(|| StoreError::CorruptFormat(protein_id.to_string()))?;
            Ok(StructureFile { chain, format, content })
        })
        .transpose()
    }

    /// Delete by id. Returns whether a row was removed.
    ///
    /// # Errors
    /// [`StoreError::Sqlite`] on failure.
    pub fn delete(&self, id: &str) -> Result<bool, StoreError> {
        Ok(self.conn.execute("DELETE FROM proteins WHERE id = ?1", params![id])? > 0)
    }

    fn decode(
        (id, name, sequence, created): (String, String, String, String),
    ) -> Result<StoredProtein, StoreError> {
        let sequence = Sequence::parse(&sequence)
            .map_err(|source| StoreError::CorruptSequence { id: id.clone(), source })?;
        let created_at = DateTime::parse_from_rfc3339(&created)
            .map_err(|_| StoreError::CorruptTimestamp(id.clone()))?
            .with_timezone(&Utc);
        Ok(StoredProtein { id, protein: Protein::new(name, sequence), created_at })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn protein(name: &str, seq: &str) -> Protein {
        Protein::new(name, Sequence::parse(seq).unwrap())
    }

    #[test]
    fn save_get_list_delete() {
        let db = Database::open_in_memory().unwrap();
        let id = db.save(&protein("one", "MKTV")).unwrap();
        db.save(&protein("two", "ACDE")).unwrap();
        assert_eq!(db.get(&id).unwrap().unwrap().protein.name, "one");
        assert_eq!(db.list().unwrap().len(), 2);
        assert!(db.delete(&id).unwrap());
        assert!(db.get(&id).unwrap().is_none());
    }

    #[test]
    fn structures_attach_and_cascade_on_delete() {
        let db = Database::open_in_memory().unwrap();
        let id = db.save(&protein("p", "MKTV")).unwrap();
        let file =
            StructureFile { chain: "A".into(), format: StructureFormat::Pdb, content: "ATOM ...".into() };
        db.save_structure(&id, &file).unwrap();
        assert_eq!(db.structure(&id).unwrap(), Some(file));
        db.delete(&id).unwrap();
        assert_eq!(db.structure(&id).unwrap(), None);
        assert!(db
            .save_structure(
                "missing",
                &StructureFile { chain: "A".into(), format: StructureFormat::Pdb, content: String::new() }
            )
            .is_err());
    }

    #[test]
    fn saving_same_sequence_renames_instead_of_duplicating() {
        let db = Database::open_in_memory().unwrap();
        db.save(&protein("old", "MKTV")).unwrap();
        let id = db.save(&protein("new", "MKTV")).unwrap();
        assert_eq!(db.list().unwrap().len(), 1);
        assert_eq!(db.get(&id).unwrap().unwrap().protein.name, "new");
    }
}
