//! Root component and shared workspace state.

use crate::components::{AnalysisPanel, Explorer, PredictionPanel, SequenceRuler, ViewerPane};
use crate::ipc;
use leptos::prelude::*;
use leptos::task::spawn_local;
use opendna_app_protocol::{ProteinDetail, ProteinSummary, SequenceReport, StructureReport};

/// Message shown in the status bar.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Status {
    pub text: String,
    pub is_error: bool,
}

/// State shared by every component, provided through context.
#[derive(Clone, Copy)]
pub struct Workspace {
    pub proteins: RwSignal<Vec<ProteinSummary>>,
    pub selected: RwSignal<Option<String>>,
    pub detail: RwSignal<Option<ProteinDetail>>,
    pub sequence_report: RwSignal<Option<Result<SequenceReport, String>>>,
    pub structure_report: RwSignal<Option<Result<StructureReport, String>>>,
    pub status: RwSignal<Status>,
    pub busy: RwSignal<bool>,
}

impl Workspace {
    fn new() -> Self {
        Self {
            proteins: RwSignal::new(Vec::new()),
            selected: RwSignal::new(None),
            detail: RwSignal::new(None),
            sequence_report: RwSignal::new(None),
            structure_report: RwSignal::new(None),
            status: RwSignal::new(Status::default()),
            busy: RwSignal::new(false),
        }
    }

    pub fn info(self, text: impl Into<String>) {
        self.status.set(Status { text: text.into(), is_error: false });
    }

    pub fn error(self, text: impl Into<String>) {
        self.status.set(Status { text: text.into(), is_error: true });
    }

    pub fn refresh(self) {
        spawn_local(async move {
            match ipc::list_proteins().await {
                Ok(list) => self.proteins.set(list),
                Err(e) => self.error(e),
            }
        });
    }

    /// Select a protein and load its detail and analyses. Responses for a
    /// protein that is no longer selected are discarded.
    pub fn select(self, id: String) {
        self.selected.set(Some(id.clone()));
        self.detail.set(None);
        self.sequence_report.set(None);
        self.structure_report.set(None);
        let still_selected = move |id: &str| self.selected.get_untracked().as_deref() == Some(id);

        spawn_local(async move {
            let detail = match ipc::get_protein(&id).await {
                Ok(d) => d,
                Err(e) => return self.error(e),
            };
            if !still_selected(&id) {
                return;
            }
            let has_structure = detail.structure.is_some();
            self.detail.set(Some(detail));

            let seq = ipc::analyze_sequence(&id).await;
            if still_selected(&id) {
                self.sequence_report.set(Some(seq));
            }
            if has_structure {
                let report = ipc::analyze_structure(&id).await;
                if still_selected(&id) {
                    self.structure_report.set(Some(report));
                }
            }
        });
    }

    pub fn import(self) {
        self.busy.set(true);
        spawn_local(async move {
            match ipc::import_with_dialog().await {
                Ok(added) if added.is_empty() => {}
                Ok(added) => {
                    self.info(format!(
                        "Imported {} {}",
                        added.len(),
                        if added.len() == 1 { "chain" } else { "chains" }
                    ));
                    if let Some(first) = added.first() {
                        self.select(first.id.clone());
                    }
                    self.refresh();
                }
                Err(e) => self.error(e),
            }
            self.busy.set(false);
        });
    }

    pub fn delete(self, id: String) {
        spawn_local(async move {
            if let Err(e) = ipc::delete_protein(&id).await {
                return self.error(e);
            }
            if self.selected.get_untracked().as_deref() == Some(id.as_str()) {
                self.selected.set(None);
                self.detail.set(None);
                self.sequence_report.set(None);
                self.structure_report.set(None);
            }
            self.refresh();
        });
    }
}

#[component]
pub fn App() -> impl IntoView {
    let ws = Workspace::new();
    provide_context(ws);
    ws.refresh();

    view! {
        <div class="app">
            <header class="topbar">
                <span class="brand">"OpenDNA"</span>
                <span class="spacer"></span>
                <button class="primary" disabled=move || ws.busy.get() on:click=move |_| ws.import()>
                    "Import…"
                </button>
            </header>
            <aside class="sidebar">
                <Explorer />
            </aside>
            <main class="main">
                <ViewerPane />
                <SequenceRuler />
            </main>
            <aside class="panel">
                <AnalysisPanel />
                <PredictionPanel />
            </aside>
            <footer class="statusbar">
                <span>
                    {move || {
                        let n = ws.proteins.get().len();
                        format!("{n} {}", if n == 1 { "protein" } else { "proteins" })
                    }}
                </span>
                <span class:error=move || ws.status.get().is_error>{move || ws.status.get().text}</span>
            </footer>
        </div>
    }
}
