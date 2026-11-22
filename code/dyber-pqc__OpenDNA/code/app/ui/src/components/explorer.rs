use crate::app::Workspace;
use leptos::prelude::*;
use opendna_app_protocol::ProteinSummary;

#[component]
pub fn Explorer() -> impl IntoView {
    let ws = expect_context::<Workspace>();
    view! {
        <div class="explorer-header">
            <h2 class="section-title">"Proteins"</h2>
        </div>
        <Show
            when=move || !ws.proteins.get().is_empty()
            fallback=|| view! { <p class="empty">"Nothing imported yet. Use Import… to open a PDB, mmCIF or FASTA file."</p> }
        >
            <ul class="explorer-list">
                <For
                    each=move || ws.proteins.get()
                    key=|p| p.id.clone()
                    children=move |p: ProteinSummary| view! { <ExplorerItem protein=p /> }
                />
            </ul>
        </Show>
    }
}

#[component]
fn ExplorerItem(protein: ProteinSummary) -> impl IntoView {
    let ws = expect_context::<Workspace>();
    // Deleting takes two clicks: the first arms the button, the second deletes.
    let armed = RwSignal::new(false);
    let id = protein.id.clone();
    let is_selected = {
        let id = id.clone();
        move || ws.selected.get().as_deref() == Some(id.as_str())
    };
    let select_id = id.clone();

    view! {
        <li
            class="explorer-item"
            class:selected=is_selected
            title=protein.name.clone()
            on:click=move |_| ws.select(select_id.clone())
            on:mouseleave=move |_| armed.set(false)
        >
            <span class="name">{protein.name.clone()}</span>
            {protein.has_structure.then(|| view! { <span class="badge" title="Has 3D coordinates">"3D"</span> })}
            <span class="meta">{format!("{} aa", protein.length)}</span>
            <button
                class="danger"
                class:armed=move || armed.get()
                title="Delete"
                on:click=move |ev| {
                    ev.stop_propagation();
                    if armed.get() {
                        ws.delete(id.clone());
                    } else {
                        armed.set(true);
                    }
                }
            >
                {move || if armed.get() { "Delete" } else { "×" }}
            </button>
        </li>
    }
}
