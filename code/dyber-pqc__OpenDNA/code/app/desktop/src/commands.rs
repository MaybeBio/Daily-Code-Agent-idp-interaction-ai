//! Tauri command handlers. Each is a thin adapter: validate input, call into
//! the library crates, convert the result to a protocol type.

use crate::{import, predict::Predictor, reports, state::AppState};
use opendna_app_protocol::{
    BurialPrediction, CommandError, DisorderPrediction, ModelStatus, ProteinDetail, ProteinSummary,
    SequenceReport, SolubilityPrediction, StabilityScan, StructurePayload, StructureReport,
};
use std::fmt::Display;
use std::path::PathBuf;
use tauri::{AppHandle, State};
use tauri_plugin_dialog::DialogExt;

type CommandResult<T> = Result<T, CommandError>;

pub(crate) fn fail(e: impl Display) -> CommandError {
    CommandError { message: e.to_string() }
}

fn summary(
    db: &opendna_store::Database,
    stored: &opendna_store::StoredProtein,
) -> CommandResult<ProteinSummary> {
    Ok(ProteinSummary {
        id: stored.id.clone(),
        name: stored.protein.name.clone(),
        length: stored.protein.sequence.len(),
        has_structure: db.structure(&stored.id).map_err(fail)?.is_some(),
    })
}

#[tauri::command]
pub fn list_proteins(state: State<'_, AppState>) -> CommandResult<Vec<ProteinSummary>> {
    let db = state.db();
    db.list().map_err(fail)?.iter().map(|p| summary(&db, p)).collect()
}

/// Import a file. With no `path`, shows a native open dialog; cancelling
/// returns an empty list. Runs off the main thread so the dialog can block.
#[tauri::command(async)]
pub fn import_file(
    app: AppHandle,
    state: State<'_, AppState>,
    path: Option<String>,
) -> CommandResult<Vec<ProteinSummary>> {
    let path = if let Some(p) = path {
        PathBuf::from(p)
    } else {
        let picked = app
            .dialog()
            .file()
            .set_title("Import structure or sequence")
            .add_filter("Structures and sequences", import::EXTENSIONS)
            .blocking_pick_file();
        match picked {
            Some(file) => file.into_path().map_err(fail)?,
            None => return Ok(Vec::new()),
        }
    };
    let db = state.db();
    let ids = import::import_path(&db, &path)?;
    ids.iter()
        .map(|id| {
            let stored = db.get(id).map_err(fail)?.ok_or_else(|| fail("imported protein vanished"))?;
            summary(&db, &stored)
        })
        .collect()
}

#[tauri::command]
pub fn get_protein(state: State<'_, AppState>, id: String) -> CommandResult<ProteinDetail> {
    let db = state.db();
    let stored = db.get(&id).map_err(fail)?.ok_or_else(|| fail(format!("no protein with id {id}")))?;
    let structure = db.structure(&id).map_err(fail)?.map(|f| StructurePayload {
        format: f.format.as_str().to_string(),
        chain: f.chain,
        content: f.content,
    });
    Ok(ProteinDetail {
        summary: summary(&db, &stored)?,
        sequence: stored.protein.sequence.to_string(),
        structure,
    })
}

#[tauri::command]
pub fn delete_protein(state: State<'_, AppState>, id: String) -> CommandResult<()> {
    state.db().delete(&id).map_err(fail).map(|_| ())
}

#[tauri::command]
pub fn analyze_sequence(state: State<'_, AppState>, id: String) -> CommandResult<SequenceReport> {
    let stored =
        state.db().get(&id).map_err(fail)?.ok_or_else(|| fail(format!("no protein with id {id}")))?;
    reports::sequence(&stored.protein.sequence)
}

#[tauri::command(async)]
pub fn analyze_structure(state: State<'_, AppState>, id: String) -> CommandResult<StructureReport> {
    let file =
        state.db().structure(&id).map_err(fail)?.ok_or_else(|| fail("this protein has no structure"))?;
    reports::structure(&file)
}

/// What the prediction models can currently do. Loads nothing and downloads
/// nothing, so the UI can call it at startup.
///
/// Infallible by design: a head that cannot be used is a value in the result,
/// not an error, so the UI can report *which* head is unusable and why rather
/// than losing that detail in a single failure.
#[tauri::command(async)]
pub fn model_status(predictor: State<'_, Predictor>) -> ModelStatus {
    predictor.status()
}

/// Score all 19 substitutions at one position.
///
/// The first prediction of a session loads the embedding model, which may
/// download it; later calls reuse it.
#[tauri::command(async)]
pub fn predict_stability(
    state: State<'_, AppState>,
    predictor: State<'_, Predictor>,
    id: String,
    position: usize,
) -> CommandResult<StabilityScan> {
    predictor.stability_scan(&sequence_of(&state, &id)?, position)
}

#[tauri::command(async)]
pub fn predict_solubility(
    state: State<'_, AppState>,
    predictor: State<'_, Predictor>,
    id: String,
) -> CommandResult<SolubilityPrediction> {
    predictor.solubility(&sequence_of(&state, &id)?)
}

#[tauri::command(async)]
pub fn predict_disorder(
    state: State<'_, AppState>,
    predictor: State<'_, Predictor>,
    id: String,
) -> CommandResult<DisorderPrediction> {
    predictor.disorder(&sequence_of(&state, &id)?)
}

#[tauri::command(async)]
pub fn predict_burial(
    state: State<'_, AppState>,
    predictor: State<'_, Predictor>,
    id: String,
) -> CommandResult<BurialPrediction> {
    predictor.burial(&sequence_of(&state, &id)?)
}

/// Save text the UI has rendered, to a file the user picks.
///
/// The UI renders the table because it already holds the prediction and the
/// limitations that came with it; this only writes bytes. Keeping the dialog
/// in Rust means the front end needs no filesystem permission of its own —
/// the window's capability grants `core:default` and nothing more, so a bug
/// in the UI cannot reach the disk except through a path the user chose in a
/// native dialog.
///
/// Returns the path written, or `None` if the user cancelled, which is not an
/// error.
#[tauri::command(async)]
pub fn save_text_file(app: AppHandle, name: String, contents: String) -> CommandResult<Option<String>> {
    let picked = app
        .dialog()
        .file()
        .set_title("Save predictions")
        .set_file_name(&name)
        .add_filter("Comma-separated values", &["csv"])
        .blocking_save_file();
    let Some(file) = picked else {
        return Ok(None);
    };
    let path = file.into_path().map_err(fail)?;
    std::fs::write(&path, contents).map_err(fail)?;
    Ok(Some(path.display().to_string()))
}

/// Look up a stored protein's sequence, releasing the database lock before
/// the long-running model work starts.
fn sequence_of(state: &State<'_, AppState>, id: &str) -> CommandResult<String> {
    let stored = state.db().get(id).map_err(fail)?.ok_or_else(|| fail(format!("no protein with id {id}")))?;
    Ok(stored.protein.sequence.to_string())
}
