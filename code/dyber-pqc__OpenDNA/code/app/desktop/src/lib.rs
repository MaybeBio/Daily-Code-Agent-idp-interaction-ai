//! `OpenDNA` desktop backend.
//!
//! The UI talks to this process only through the Tauri commands in
//! [`commands`]; there is no HTTP server and no open port.

#![cfg_attr(test, allow(clippy::unwrap_used))]

mod commands;
mod import;
mod predict;
mod reports;
mod state;

use tauri::Manager;

/// Start the application.
///
/// # Panics
/// If the Tauri runtime fails to start or the data directory can't be created;
/// both are unrecoverable at launch.
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .setup(|app| {
            let dir = app.path().app_data_dir()?;
            std::fs::create_dir_all(&dir)?;
            app.manage(state::AppState::open(&dir.join("opendna.sqlite3"))?);
            // Models are found now but not loaded: the embedding model is a
            // few hundred megabytes, so the first prediction pays for it
            // rather than every launch.
            app.manage(predict::Predictor::new(predict::models_dir(&dir)));
            Ok(())
        })
        .invoke_handler(tauri::generate_handler![
            commands::list_proteins,
            commands::import_file,
            commands::get_protein,
            commands::delete_protein,
            commands::analyze_sequence,
            commands::analyze_structure,
            commands::model_status,
            commands::predict_stability,
            commands::predict_solubility,
            commands::predict_disorder,
            commands::predict_burial,
            commands::save_text_file,
        ])
        .run(tauri::generate_context!())
        .expect("failed to start OpenDNA");
}
