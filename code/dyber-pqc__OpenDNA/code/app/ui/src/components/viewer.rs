use crate::{app::Workspace, bridge};
use leptos::{html, prelude::*, task::spawn_local};

#[component]
pub fn ViewerPane() -> impl IntoView {
    let ws = expect_context::<Workspace>();
    let host: NodeRef<html::Div> = NodeRef::new();

    // Reload the viewer whenever the selected protein's structure changes.
    Effect::new(move |_| {
        let structure = ws.detail.with(|d| d.as_ref().and_then(|d| d.structure.clone()));
        let Some(el) = host.get() else { return };
        match structure {
            Some(s) => spawn_local(async move {
                if let Err(e) = bridge::load(&el, &s.content, &s.format).await {
                    ws.error(e);
                }
            }),
            None => bridge::clear(),
        }
    });

    let placeholder = move || {
        let detail = ws.detail.get();
        match detail {
            None if ws.selected.get().is_some() => Some("Loading…"),
            None => Some("Select a protein to view it."),
            Some(d) if d.structure.is_none() => {
                Some("No 3D coordinates for this protein. Import a PDB or mmCIF file to view a structure.")
            }
            Some(_) => None,
        }
    };

    view! {
        <div class="viewer-frame">
            <div class="viewer-host" node_ref=host></div>
            {move || placeholder().map(|text| view! { <div class="viewer-placeholder">{text}</div> })}
        </div>
    }
}
