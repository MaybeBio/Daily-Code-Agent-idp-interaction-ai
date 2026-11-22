use crate::{app::Workspace, bridge};
use leptos::prelude::*;
use std::fmt::Write;

const ROW: usize = 50;
const BLOCK: usize = 10;

#[derive(Clone)]
struct Cell {
    letter: char,
    dssp: char,
    /// Residue number in the structure, when the ruler is built from one.
    number: Option<i32>,
    tooltip: String,
}

fn dssp_name(code: char) -> &'static str {
    match code {
        'H' => "α-helix",
        'G' => "3₁₀ helix",
        'I' => "π helix",
        'E' => "strand",
        'B' => "bridge",
        'T' => "turn",
        'S' => "bend",
        _ => "loop",
    }
}

/// Residue cells: from the structure report when available (so numbering and
/// secondary structure match the coordinates), otherwise from the sequence.
fn cells(ws: Workspace) -> Vec<Cell> {
    if let Some(Ok(report)) = ws.structure_report.get() {
        return report
            .residues
            .iter()
            .map(|r| {
                let mut tooltip = format!("{}{} · {}", r.amino_acid, r.label, dssp_name(r.dssp));
                if let Some(s) = r.relative_sasa {
                    let _ = write!(tooltip, " · {:.0}% exposed", s * 100.0);
                }
                if let (Some(phi), Some(psi)) = (r.phi, r.psi) {
                    let _ = write!(tooltip, " · φ {phi:.0}° ψ {psi:.0}°");
                }
                Cell { letter: r.amino_acid, dssp: r.dssp, number: Some(r.number), tooltip }
            })
            .collect();
    }
    ws.detail
        .get()
        .map(|d| {
            d.sequence
                .chars()
                .enumerate()
                .map(|(i, c)| Cell { letter: c, dssp: '-', number: None, tooltip: format!("{c}{}", i + 1) })
                .collect()
        })
        .unwrap_or_default()
}

#[component]
pub fn SequenceRuler() -> impl IntoView {
    let ws = expect_context::<Workspace>();
    let chain =
        move || ws.detail.with(|d| d.as_ref().and_then(|d| d.structure.as_ref().map(|s| s.chain.clone())));

    let rows = move || {
        let cells = cells(ws);
        let chain = chain();
        cells
            .chunks(ROW)
            .enumerate()
            .map(|(row, row_cells)| {
                let blocks = row_cells
                    .chunks(BLOCK)
                    .map(|block| {
                        block
                            .iter()
                            .map(|c| {
                                let (enter, leave, click) = (chain.clone(), chain.clone(), chain.clone());
                                let n = c.number;
                                view! {
                                    <span
                                        class=format!("res ss-{}", c.dssp)
                                        title=c.tooltip.clone()
                                        on:mouseenter=move |_| if let (Some(ch), Some(n)) = (&enter, n) { bridge::highlight(ch, Some(n)) }
                                        on:mouseleave=move |_| if let Some(ch) = &leave { bridge::highlight(ch, None) }
                                        on:click=move |_| if let (Some(ch), Some(n)) = (&click, n) { bridge::focus(ch, n) }
                                    >
                                        {c.letter.to_string()}
                                    </span>
                                }
                            })
                            .collect_view()
                    })
                    .map(|residues| view! { <span class="ruler-block">{residues}</span> })
                    .collect_view();
                view! {
                    <div class="ruler-row">
                        <span class="ruler-index">{row * ROW + 1}</span>
                        <span>{blocks}</span>
                    </div>
                }
            })
            .collect_view()
    };

    view! {
        <Show when=move || ws.detail.with(Option::is_some)>
            <section class="ruler">
                <Show when=move || matches!(ws.structure_report.get(), Some(Ok(_)))>
                    <div class="legend">
                        <span style="--swatch: var(--ss-helix)">"Helix"</span>
                        <span style="--swatch: var(--ss-strand)">"Strand"</span>
                        <span style="--swatch: var(--ss-turn)">"Turn / bend"</span>
                        <span>"Hover to highlight, click to focus"</span>
                    </div>
                </Show>
                {rows}
            </section>
        </Show>
    }
}
