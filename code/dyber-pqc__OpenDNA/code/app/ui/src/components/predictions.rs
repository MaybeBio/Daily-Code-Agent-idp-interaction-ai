//! Predictions from the trained heads.
//!
//! Every number shown here comes from a trained model. When a head is not
//! installed the panel says so and names what is missing; it never falls back
//! to a heuristic, because a plausible-looking number with no model behind it
//! is worse than no number at all.

use crate::{app::Workspace, ipc};
use leptos::prelude::*;
use leptos::task::spawn_local;
use opendna_app_protocol::{
    BurialPrediction, DisorderPrediction, HeadStatus, ModelStatus, SolubilityPrediction, StabilityScan,
};
use std::fmt::Write as _;

/// A request that has been sent but not yet answered, or its result.
type Pending<T> = RwSignal<Option<Result<T, String>>>;

#[component]
pub fn PredictionPanel() -> impl IntoView {
    let ws = expect_context::<Workspace>();
    let status: RwSignal<Option<Result<ModelStatus, String>>> = RwSignal::new(None);
    let solubility: Pending<SolubilityPrediction> = RwSignal::new(None);
    let disorder: Pending<DisorderPrediction> = RwSignal::new(None);
    let burial: Pending<BurialPrediction> = RwSignal::new(None);
    let stability: Pending<StabilityScan> = RwSignal::new(None);
    let running = RwSignal::new(false);
    let position = RwSignal::new(1usize);

    spawn_local(async move {
        status.set(Some(ipc::model_status().await));
    });

    // A prediction belongs to the protein it was requested for; switching
    // selection must not leave another protein's numbers on screen.
    Effect::new(move |_| {
        ws.selected.track();
        solubility.set(None);
        disorder.set(None);
        burial.set(None);
        stability.set(None);
        position.set(1);
    });

    view! {
        <section class="section">
            <h2 class="section-title">"Predictions"</h2>
            {move || match status.get() {
                None => view! { <p class="empty">"Checking for models…"</p> }.into_any(),
                Some(Err(e)) => view! { <p class="note error">{e}</p> }.into_any(),
                Some(Ok(s)) => view! { <Models status=s ws=ws solubility=solubility
                    disorder=disorder burial=burial stability=stability running=running
                    position=position /> }
                    .into_any(),
            }}
        </section>
    }
}

#[component]
fn Models(
    status: ModelStatus,
    ws: Workspace,
    solubility: Pending<SolubilityPrediction>,
    disorder: Pending<DisorderPrediction>,
    burial: Pending<BurialPrediction>,
    stability: Pending<StabilityScan>,
    running: RwSignal<bool>,
    position: RwSignal<usize>,
) -> impl IntoView {
    let any_ready = status.stability.is_ready()
        || status.solubility.is_ready()
        || status.disorder.is_ready()
        || status.burial.is_ready();
    let models_dir = status.models_dir.clone();
    let warn_download = !status.embedding_cached;
    let embedding_model = status.embedding_model.clone();
    let device = status.device.clone();

    let selected = move || ws.selected.get();
    let sequence_length = move || ws.detail.get().map_or(0, |d| d.sequence.len());

    // Exported tables are rendered from whatever is on screen, so a file can
    // never hold a different number than the panel showed. Each carries the
    // head version and the card's limitations.
    let (sol_head, dis_head, stab_head, bur_head) =
        (status.solubility.clone(), status.disorder.clone(), status.stability.clone(), status.burial.clone());
    let name_of = move || ws.detail.get().map_or_else(|| "protein".to_string(), |d| d.summary.name.clone());

    let for_csv = status.clone();
    let head = for_csv.solubility.clone();
    let sol_csv = Memo::new(move |_| {
        let Some(Ok(p)) = solubility.get() else { return String::new() };
        let mut out = provenance("solubility", &for_csv, &head);
        out.push_str("name,solubility_percent\n");
        let _ = writeln!(out, "{},{:.1}", field(&name_of()), p.percent);
        out
    });

    let for_csv = status.clone();
    let head = for_csv.disorder.clone();
    let dis_csv = Memo::new(move |_| {
        let Some(Ok(p)) = disorder.get() else { return String::new() };
        let sequence = ws.detail.get().map(|d| d.sequence.clone()).unwrap_or_default();
        let residues: Vec<char> = sequence.chars().collect();
        let mut out = provenance("disorder", &for_csv, &head);
        let _ = writeln!(out, "# threshold: {}", p.threshold);
        if p.windowed {
            out.push_str("# covered by overlapping windows: longer than the model's limit\n");
        }
        out.push_str("name,position,residue,probability,disordered\n");
        let who = field(&name_of());
        for (i, probability) in p.probabilities.iter().enumerate() {
            let _ = writeln!(
                out,
                "{who},{},{},{probability:.4},{}",
                i + 1,
                residues.get(i).copied().unwrap_or('?'),
                u8::from(*probability >= p.threshold)
            );
        }
        out
    });

    let for_csv = status.clone();
    let head = for_csv.burial.clone();
    let bur_csv = Memo::new(move |_| {
        let Some(Ok(p)) = burial.get() else { return String::new() };
        let sequence = ws.detail.get().map(|d| d.sequence.clone()).unwrap_or_default();
        let residues: Vec<char> = sequence.chars().collect();
        let mut out = provenance("burial", &for_csv, &head);
        let _ = writeln!(out, "# buried at or below: {}", p.buried_at);
        if p.windowed {
            out.push_str("# covered by overlapping windows: longer than the model's limit\n");
        }
        out.push_str("name,position,residue,accessibility,buried\n");
        let who = field(&name_of());
        for (i, value) in p.accessibility.iter().enumerate() {
            let _ = writeln!(
                out,
                "{who},{},{},{value:.4},{}",
                i + 1,
                residues.get(i).copied().unwrap_or('?'),
                u8::from(*value <= p.buried_at)
            );
        }
        out
    });

    let for_csv = status.clone();
    let head = for_csv.stability.clone();
    let stab_csv = Memo::new(move |_| {
        let Some(Ok(s)) = stability.get() else { return String::new() };
        let mut out = provenance("stability", &for_csv, &head);
        out.push_str("name,position,wild_type,mutant,ddg_kcal_per_mol,log_odds\n");
        let who = field(&name_of());
        for score in &s.scores {
            let _ = writeln!(
                out,
                "{who},{},{},{},{:+.3},{:+.3}",
                s.position, s.wild_type, score.mutant, score.ddg, score.log_odds
            );
        }
        out
    });

    // Each request runs one at a time: they share one embedding model, and a
    // second concurrent call would simply wait on the same lock in the
    // backend while looking responsive here.
    let run = move |kind: Kind| {
        let Some(id) = selected() else { return };
        running.set(true);
        spawn_local(async move {
            match kind {
                Kind::Solubility => solubility.set(Some(ipc::predict_solubility(&id).await)),
                Kind::Disorder => disorder.set(Some(ipc::predict_disorder(&id).await)),
                Kind::Burial => burial.set(Some(ipc::predict_burial(&id).await)),
                Kind::Stability => {
                    let at = position.get_untracked();
                    stability.set(Some(ipc::predict_stability(&id, at).await));
                }
            }
            running.set(false);
        });
    };

    view! {
        {(!any_ready).then(|| view! {
            <p class="note">
                "No trained heads are installed, so predictions are unavailable. "
                "OpenDNA does not estimate these properties without a model. Put the "
                "trained files in "<code>{models_dir.clone()}</code>", or set "
                <code>"OPENDNA_MODELS"</code>" to where they are."
            </p>
        })}
        {any_ready.then(move || view! {
            {warn_download.then(|| view! {
                <p class="note">
                    "The first prediction downloads the embedding model ("
                    {embedding_model.clone()}"), several hundred megabytes. "
                    "After that it is cached and works offline."
                </p>
            })}

            <StatusRow label="Stability (ΔΔG)" status=status.stability.clone() />
            <StatusRow label="Solubility" status=status.solubility.clone() />
            <StatusRow label="Disorder" status=status.disorder.clone() />
            <StatusRow label="Accessibility" status=status.burial.clone() />

            {status.solubility.is_ready().then(move || view! {
                <div class="predict-block">
                    <button class="secondary"
                        disabled=move || running.get() || selected().is_none()
                        on:click=move |_| run(Kind::Solubility)>
                        "Predict solubility"
                    </button>
                    {move || solubility.get().map(|r| match r {
                        Err(e) => view! { <p class="note error">{e}</p> }.into_any(),
                        Ok(p) => view! {
                            <table class="kv">
                                <tbody>
                                    <tr>
                                        <th>"Predicted solubility"</th>
                                        <td>{format!("{:.1}", p.percent)}<span class="unit">"%"</span></td>
                                    </tr>
                                </tbody>
                            </table>
                            <ExportButton name="solubility.csv".to_string()
                                contents=sol_csv />
                            <Limitations head=sol_head.clone() />
                        }.into_any(),
                    })}
                </div>
            })}

            {status.disorder.is_ready().then(move || view! {
                <div class="predict-block">
                    <button class="secondary"
                        disabled=move || running.get() || selected().is_none()
                        on:click=move |_| run(Kind::Disorder)>
                        "Predict disorder"
                    </button>
                    {move || disorder.get().map(|r| match r {
                        Err(e) => view! { <p class="note error">{e}</p> }.into_any(),
                        Ok(p) => view! {
                            <DisorderTrack prediction=p />
                            <ExportButton name="disorder.csv".to_string() contents=dis_csv />
                            <Limitations head=dis_head.clone() />
                        }.into_any(),
                    })}
                </div>
            })}

            {status.burial.is_ready().then(move || view! {
                <div class="predict-block">
                    <button class="secondary"
                        disabled=move || running.get() || selected().is_none()
                        on:click=move |_| run(Kind::Burial)>
                        "Predict solvent accessibility"
                    </button>
                    {move || burial.get().map(|r| match r {
                        Err(e) => view! { <p class="note error">{e}</p> }.into_any(),
                        Ok(p) => view! {
                            <BurialTrack prediction=p />
                            <ExportButton name="burial.csv".to_string() contents=bur_csv />
                            <Limitations head=bur_head.clone() />
                        }.into_any(),
                    })}
                </div>
            })}

            {status.stability.is_ready().then(move || view! {
                <div class="predict-block">
                    <label class="field">
                        "Scan position"
                        <input type="number" min="1"
                            max=move || sequence_length().max(1).to_string()
                            prop:value=move || position.get().to_string()
                            on:input=move |ev| {
                                if let Ok(v) = event_target_value(&ev).parse::<usize>() {
                                    position.set(v.max(1));
                                }
                            } />
                    </label>
                    <button class="secondary"
                        disabled=move || running.get() || selected().is_none()
                        on:click=move |_| run(Kind::Stability)>
                        "Score substitutions"
                    </button>
                    {move || stability.get().map(|r| match r {
                        Err(e) => view! { <p class="note error">{e}</p> }.into_any(),
                        Ok(s) => view! {
                            <StabilityTable scan=s />
                            <ExportButton name="stability.csv".to_string() contents=stab_csv />
                            <Limitations head=stab_head.clone() />
                        }.into_any(),
                    })}
                </div>
            })}

            <p class="note">"Predictions run on "{device.clone()}"."</p>
        })}
    }
}

#[derive(Clone, Copy)]
enum Kind {
    Solubility,
    Disorder,
    Stability,
    Burial,
}

#[component]
fn StatusRow(label: &'static str, status: HeadStatus) -> impl IntoView {
    let (text, class) = match &status {
        HeadStatus::Ready { version, .. } => (format!("v{version}"), "ok"),
        HeadStatus::Missing => ("not installed".to_string(), "muted"),
        HeadStatus::Unusable { reason } => (reason.clone(), "error"),
    };
    view! {
        <div class="model-row">
            <span class="model-name">{label}</span>
            <span class=format!("model-state {class}")>{text}</span>
        </div>
    }
}

/// Per-residue disorder as a bar per residue, with the calling threshold drawn
/// across it so a reader can see how marginal each call is.
#[component]
fn DisorderTrack(prediction: DisorderPrediction) -> impl IntoView {
    let (w, h) = (300.0_f64, 48.0_f64);
    let n = prediction.probabilities.len();
    let step = if n > 0 { w / n as f64 } else { 0.0 };
    let threshold_y = h - prediction.threshold * h;
    let disordered = prediction.probabilities.iter().filter(|p| **p >= prediction.threshold).count();
    let fraction = if n > 0 { disordered as f64 / n as f64 * 100.0 } else { 0.0 };

    let bars: Vec<_> = prediction
        .probabilities
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let height = p * h;
            view! {
                <rect x=format!("{:.2}", i as f64 * step) y=format!("{:.2}", h - height)
                    width=format!("{:.2}", step.max(0.5)) height=format!("{:.2}", height)
                    class={if *p >= prediction.threshold { "bar disordered" } else { "bar" }} />
            }
        })
        .collect();

    view! {
        <svg class="track" viewBox=format!("0 0 {w} {h}") preserveAspectRatio="none" role="img"
            aria-label="Per-residue disorder probability">
            {bars}
            <line x1="0" x2=format!("{w}") y1=format!("{threshold_y:.2}") y2=format!("{threshold_y:.2}")
                class="threshold" />
        </svg>
        <table class="kv">
            <tbody>
                <tr>
                    <th>"Residues called disordered"</th>
                    <td>{disordered}" / "{n}<span class="unit">{format!("({fraction:.0}%)")}</span></td>
                </tr>
            </tbody>
        </table>
        <p class="note">
            "Threshold "{format!("{:.2}", prediction.threshold)}". "
            {prediction.windowed.then_some("This sequence is longer than the embedding model's limit and was covered by overlapping windows. ")}
            "Trained on DisProt, which records disorder that has been demonstrated: unannotated \
             residues count as ordered, so the negative class is not exhaustive."
        </p>
    }
}

/// Per-residue accessibility as a bar per residue, with the buried threshold
/// drawn across it.
///
/// Drawn the opposite way round from the disorder track: a *short* bar is a
/// buried residue, because the quantity is exposure. The threshold line is
/// therefore near the bottom, and what sits under it is the core.
#[component]
fn BurialTrack(prediction: BurialPrediction) -> impl IntoView {
    let (w, h) = (300.0_f64, 48.0_f64);
    let n = prediction.accessibility.len();
    let step = if n > 0 { w / n as f64 } else { 0.0 };
    // Accessibility can exceed 1, so the track is scaled to whatever the
    // tallest bar is rather than clipping it out of sight.
    let ceiling = prediction.accessibility.iter().copied().fold(1.0_f64, f64::max);
    let threshold_y = h - prediction.buried_at / ceiling * h;
    let buried = prediction.accessibility.iter().filter(|v| **v <= prediction.buried_at).count();
    let fraction = if n > 0 { buried as f64 / n as f64 * 100.0 } else { 0.0 };

    let bars: Vec<_> = prediction
        .accessibility
        .iter()
        .enumerate()
        .map(|(i, value)| {
            let height = value / ceiling * h;
            view! {
                <rect x=format!("{:.2}", i as f64 * step) y=format!("{:.2}", h - height)
                    width=format!("{:.2}", step.max(0.5)) height=format!("{:.2}", height)
                    class={if *value <= prediction.buried_at { "bar buried" } else { "bar" }} />
            }
        })
        .collect();

    view! {
        <svg class="track" viewBox=format!("0 0 {w} {h}") preserveAspectRatio="none" role="img"
            aria-label="Per-residue relative solvent accessibility">
            {bars}
            <line x1="0" x2=format!("{w}") y1=format!("{threshold_y:.2}") y2=format!("{threshold_y:.2}")
                class="threshold" />
        </svg>
        <table class="kv">
            <tbody>
                <tr>
                    <th>"Residues buried"</th>
                    <td>{buried}" / "{n}<span class="unit">{format!("({fraction:.0}%)")}</span></td>
                </tr>
            </tbody>
        </table>
        <p class="note">
            "Relative accessibility: 0 is fully buried, 1 fully exposed. Buried at or below "
            {format!("{:.2}", prediction.buried_at)}", the usual convention. "
            {prediction.windowed.then_some("This sequence is longer than the embedding model's limit and was covered by overlapping windows. ")}
            "Accessibility describes a folded structure, so read it beside the disorder track: \
             for a region with no stable fold the number is still produced and still means \
             nothing."
        </p>
    }
}

/// The 19 substitutions at one position, most stabilising first.
#[component]
fn StabilityTable(scan: StabilityScan) -> impl IntoView {
    let rows: Vec<_> = scan
        .scores
        .iter()
        .map(|s| {
            view! {
                <tr>
                    <td class="mono">{s.mutant.to_string()}</td>
                    <td class="num">{format!("{:+.2}", s.ddg)}</td>
                    <td class="num muted">{format!("{:+.2}", s.log_odds)}</td>
                </tr>
            }
        })
        .collect();

    view! {
        <p class="scan-title">
            "Position "{scan.position}", wild type "
            <span class="mono">{scan.wild_type.to_string()}</span>
        </p>
        <table class="scan">
            <thead>
                <tr><th>"Mut"</th><th class="num">"ΔΔG"</th><th class="num">"log-odds"</th></tr>
            </thead>
            <tbody>{rows}</tbody>
        </table>
        <p class="note">
            "ΔΔG in kcal/mol; negative is destabilising. The model was trained with a squared-error \
             loss, so predictions are compressed toward the training mean and span a narrower range \
             than real measurements. Read the ranking, not the absolute value. The log-odds column \
             is the language model's own zero-shot signal, shown because disagreement between the \
             two is worth noticing."
        </p>
    }
}

/// The `#` provenance block that precedes an exported table.
///
/// The same convention the CLI uses, for the same reason: a prediction that
/// travels without its model version and its limitations becomes a number in
/// a spreadsheet that nobody can trace. The limitations come from the head's
/// model card by way of `model_status`, not from text written here, so the
/// panel, the file and the published card all say the same thing.
fn provenance(task: &str, status: &ModelStatus, head: &HeadStatus) -> String {
    let mut out = format!("# OpenDNA {task} prediction\n");
    if let Some(version) = head.version() {
        let _ = writeln!(out, "# head version: {version}");
    }
    let _ = writeln!(out, "# embedding model: {}", status.embedding_model);
    let _ = writeln!(out, "# device: {}", status.device);
    if head.limitations().is_empty() {
        out.push_str("# limitations: see the model card next to the weights\n");
    } else {
        out.push_str("# limitations:\n");
        for limitation in head.limitations() {
            let _ = writeln!(out, "#   - {limitation}");
        }
    }
    out
}

/// A CSV field, quoted if it could otherwise shift the columns.
fn field(value: &str) -> String {
    if value.contains(',') || value.contains('"') || value.contains('\n') {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

/// Offer a rendered table as a file, and report where it went.
#[component]
fn ExportButton(name: String, contents: Memo<String>) -> impl IntoView {
    let saved: RwSignal<Option<Result<Option<String>, String>>> = RwSignal::new(None);
    view! {
        <div class="export">
            <button class="secondary" on:click=move |_| {
                let (name, contents) = (name.clone(), contents.get_untracked());
                spawn_local(async move {
                    saved.set(Some(ipc::save_text_file(&name, &contents).await));
                });
            }>"Export CSV"</button>
            {move || saved.get().map(|r| match r {
                Err(e) => view! { <span class="note error">{e}</span> }.into_any(),
                // Cancelling the dialog is not a failure and gets no message.
                Ok(None) => ().into_any(),
                Ok(Some(path)) => view! { <span class="note ok">"Saved to "{path}</span> }.into_any(),
            })}
        </div>
    }
}

/// The limitations a head's card states, shown next to its numbers because
/// they determine how the numbers may be read.
#[component]
fn Limitations(head: HeadStatus) -> impl IntoView {
    let items: Vec<_> = head.limitations().iter().map(|l| view! { <li>{l.clone()}</li> }).collect();
    if items.is_empty() {
        return view! {
            <p class="note">"See this head's model card for its limitations."</p>
        }
        .into_any();
    }
    view! { <ul class="note limitations">{items}</ul> }.into_any()
}
