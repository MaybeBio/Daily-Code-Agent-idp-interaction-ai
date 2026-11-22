//! Typed wrappers over Tauri's `invoke`. Every call uses a command name from
//! `opendna_app_protocol::command` and protocol types for arguments and results.

use opendna_app_protocol::{
    command, BurialPrediction, CommandError, DisorderPrediction, ModelStatus, ProteinDetail, ProteinSummary,
    SequenceReport, SolubilityPrediction, StabilityScan, StructureReport,
};
use serde::{de::DeserializeOwned, Serialize};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = ["window", "__TAURI__", "core"], js_name = invoke, catch)]
    async fn tauri_invoke(cmd: &str, args: JsValue) -> Result<JsValue, JsValue>;
}

/// `true` when running inside the Tauri shell (`window.__TAURI__.core` exists).
fn tauri_available() -> bool {
    let get = |target: &JsValue, key: &str| js_sys::Reflect::get(target, &JsValue::from_str(key)).ok();
    get(&js_sys::global(), "__TAURI__")
        .filter(|t| !t.is_undefined())
        .and_then(|t| get(&t, "core"))
        .is_some_and(|c| !c.is_undefined())
}

async fn call<A: Serialize, R: DeserializeOwned>(cmd: &str, args: &A) -> Result<R, String> {
    if !tauri_available() {
        return Err("Not connected to the OpenDNA backend. Start the app with `cargo tauri dev`.".into());
    }
    let args = serde_wasm_bindgen::to_value(args).map_err(|e| e.to_string())?;
    match tauri_invoke(cmd, args).await {
        Ok(value) => serde_wasm_bindgen::from_value(value).map_err(|e| format!("{cmd}: bad response: {e}")),
        Err(err) => Err(serde_wasm_bindgen::from_value::<CommandError>(err.clone())
            .map_or_else(|_| format!("{cmd} failed: {err:?}"), |e| e.message)),
    }
}

#[derive(Serialize)]
struct Id<'a> {
    id: &'a str,
}

#[derive(Serialize)]
struct Import {
    path: Option<String>,
}

pub async fn list_proteins() -> Result<Vec<ProteinSummary>, String> {
    call(command::LIST_PROTEINS, &()).await
}

/// Show the native file dialog and import the chosen file.
pub async fn import_with_dialog() -> Result<Vec<ProteinSummary>, String> {
    call(command::IMPORT_FILE, &Import { path: None }).await
}

pub async fn get_protein(id: &str) -> Result<ProteinDetail, String> {
    call(command::GET_PROTEIN, &Id { id }).await
}

pub async fn delete_protein(id: &str) -> Result<(), String> {
    call(command::DELETE_PROTEIN, &Id { id }).await
}

pub async fn analyze_sequence(id: &str) -> Result<SequenceReport, String> {
    call(command::ANALYZE_SEQUENCE, &Id { id }).await
}

pub async fn analyze_structure(id: &str) -> Result<StructureReport, String> {
    call(command::ANALYZE_STRUCTURE, &Id { id }).await
}

/// Which heads are installed and whether the embedding model is cached.
/// Downloads nothing, so this is safe to call on startup.
pub async fn model_status() -> Result<ModelStatus, String> {
    call(command::MODEL_STATUS, &()).await
}

#[derive(Serialize)]
struct Scan<'a> {
    id: &'a str,
    position: usize,
}

pub async fn predict_stability(id: &str, position: usize) -> Result<StabilityScan, String> {
    call(command::PREDICT_STABILITY, &Scan { id, position }).await
}

pub async fn predict_solubility(id: &str) -> Result<SolubilityPrediction, String> {
    call(command::PREDICT_SOLUBILITY, &Id { id }).await
}

pub async fn predict_disorder(id: &str) -> Result<DisorderPrediction, String> {
    call(command::PREDICT_DISORDER, &Id { id }).await
}

pub async fn predict_burial(id: &str) -> Result<BurialPrediction, String> {
    call(command::PREDICT_BURIAL, &Id { id }).await
}

// Both fields are single words, so the key is the same whether Tauri looks
// for the camelCase or the snake_case spelling of the parameter name. The
// other payloads here are single words for the same reason.
#[derive(Serialize)]
struct SaveText<'a> {
    name: &'a str,
    contents: &'a str,
}

/// Save rendered text to a file the user picks. `Ok(None)` means they
/// cancelled, which is not a failure.
pub async fn save_text_file(name: &str, contents: &str) -> Result<Option<String>, String> {
    call(command::SAVE_TEXT_FILE, &SaveText { name, contents }).await
}
