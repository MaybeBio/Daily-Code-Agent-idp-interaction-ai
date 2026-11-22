use opendna_store::{Database, StoreError};
use std::path::Path;
use std::sync::{Mutex, MutexGuard, PoisonError};

/// Process-wide state shared by all commands.
pub struct AppState {
    db: Mutex<Database>,
}

impl AppState {
    pub fn open(path: &Path) -> Result<Self, StoreError> {
        Ok(Self { db: Mutex::new(Database::open(path)?) })
    }

    /// Lock the database. A panic in another command can't leave the database in a
    /// half-written state (`SQLite` statements are atomic), so a poisoned lock is safe to reuse.
    pub fn db(&self) -> MutexGuard<'_, Database> {
        self.db.lock().unwrap_or_else(PoisonError::into_inner)
    }
}
