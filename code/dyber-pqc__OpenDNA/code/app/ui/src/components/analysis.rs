use crate::{app::Workspace, format};
use leptos::prelude::*;
use opendna_app_protocol::{Contact, SequenceReport, StructureReport};

#[component]
pub fn AnalysisPanel() -> impl IntoView {
    let ws = expect_context::<Workspace>();
    view! {
        {move || match ws.sequence_report.get() {
            None if ws.selected.get().is_none() => view! { <p class="empty">"No protein selected."</p> }.into_any(),
            None => view! { <p class="empty">"Analysing…"</p> }.into_any(),
            Some(Err(e)) => view! {
                <section class="section">
                    <h2 class="section-title">"Sequence"</h2>
                    <p class="note">{e}</p>
                </section>
            }.into_any(),
            Some(Ok(r)) => view! { <SequenceSection report=r /> }.into_any(),
        }}
        {move || ws.structure_report.get().map(|r| match r {
            Ok(r) => view! { <StructureSection report=r /> }.into_any(),
            Err(e) => view! {
                <section class="section">
                    <h2 class="section-title">"Structure"</h2>
                    <p class="note">{e}</p>
                </section>
            }.into_any(),
        })}
    }
}

fn row(label: &'static str, value: String, unit: &'static str) -> impl IntoView {
    view! {
        <tr>
            <th>{label}</th>
            <td>{value}{(!unit.is_empty()).then(|| view! { <span class="unit">{unit}</span> })}</td>
        </tr>
    }
}

#[component]
fn SequenceSection(report: SequenceReport) -> impl IntoView {
    view! {
        <section class="section">
            <h2 class="section-title">"Sequence properties"</h2>
            <table class="kv">
                <tbody>
                    {row("Length", report.length.to_string(), "aa")}
                    {row("Molecular weight", format::number(report.molecular_weight, 1), "Da")}
                    {row("Isoelectric point", format!("{:.2}", report.isoelectric_point), "")}
                    {row("Net charge, pH 7", format!("{:+.2}", report.charge_at_ph7), "")}
                    {row("GRAVY", format!("{:.3}", report.gravy), "")}
                    {row("Aromaticity", format::percent(report.aromaticity), "%")}
                    {row("Aliphatic index", format!("{:.1}", report.aliphatic_index), "")}
                    {row("ε280, Cys reduced", format::number(report.extinction_reduced, 0), "M⁻¹cm⁻¹")}
                    {row("ε280, cystines", format::number(report.extinction_cystines, 0), "M⁻¹cm⁻¹")}
                </tbody>
            </table>
            <p class="note">"Computed with ExPASy ProtParam conventions."</p>
        </section>
        <section class="section">
            <h2 class="section-title">"Hydropathy (Kyte–Doolittle, window 9)"</h2>
            <Sparkline values=report.hydropathy />
        </section>
    }
}

#[component]
fn Sparkline(values: Vec<f64>) -> impl IntoView {
    // Fixed y-range of the Kyte–Doolittle scale so plots are comparable between proteins.
    const LIMIT: f64 = 4.5;
    let (w, h) = (300.0, 56.0);
    let y = |v: f64| h / 2.0 - (v / LIMIT) * (h / 2.0);
    let step = if values.len() > 1 { w / (values.len() - 1) as f64 } else { 0.0 };
    let points: String = values
        .iter()
        .enumerate()
        .map(|(i, v)| format!("{:.1},{:.1}", i as f64 * step, y(*v)))
        .collect::<Vec<_>>()
        .join(" ");
    view! {
        <svg class="sparkline" viewBox=format!("0 0 {w} {h}") preserveAspectRatio="none" role="img" aria-label="Hydropathy profile">
            <line class="zero" x1="0" x2=w y1=h / 2.0 y2=h / 2.0 />
            <polyline class="line" points=points />
        </svg>
    }
}

fn contact_list(contacts: Vec<Contact>) -> impl IntoView {
    view! {
        <ul class="contact-list">
            {contacts
                .into_iter()
                .map(|c| view! {
                    <li>
                        <span>{format!("{} – {}", c.a, c.b)}</span>
                        <span class="dist">{format!("{:.2} Å", c.distance)}</span>
                    </li>
                })
                .collect_view()}
        </ul>
    }
}

#[component]
fn StructureSection(report: StructureReport) -> impl IntoView {
    let n = report.residues.len().max(1) as f64;
    let share =
        |codes: &[char]| report.residues.iter().filter(|r| codes.contains(&r.dssp)).count() as f64 / n;
    let (helix, strand) = (share(&['H', 'G', 'I']), share(&['E', 'B']));
    let disulfide_count = report.disulfides.len();
    let salt_count = report.salt_bridges.len();

    view! {
        <section class="section">
            <h2 class="section-title">{format!("Structure, chain {}", report.chain)}</h2>
            <table class="kv">
                <tbody>
                    {row("Observed residues", report.residues.len().to_string(), "")}
                    {row("Helix", format::percent(helix), "%")}
                    {row("Strand", format::percent(strand), "%")}
                    {row("Radius of gyration (Cα)", format!("{:.2}", report.radius_of_gyration), "Å")}
                    {row("Solvent-accessible area", format::number(report.total_sasa, 0), "Å²")}
                </tbody>
            </table>
            <p class="note">"Secondary structure by DSSP; area by Shrake–Rupley, 1.4 Å probe, whole entry without waters."</p>
        </section>
        <section class="section">
            <h2 class="section-title">{format!("Disulfide bonds ({disulfide_count})")}</h2>
            {if disulfide_count == 0 { view! { <p class="note">"None."</p> }.into_any() } else { contact_list(report.disulfides).into_any() }}
        </section>
        <section class="section">
            <h2 class="section-title">{format!("Salt bridges ({salt_count})")}</h2>
            {if salt_count == 0 { view! { <p class="note">"None."</p> }.into_any() } else { contact_list(report.salt_bridges).into_any() }}
        </section>
    }
}
