"use client";

import { Compass, ExternalLink } from "lucide-react";
import Link from "next/link";

// The five saturating axes of the Structural Impact Score (SIS v2).
// Kept in sync with backend/scripts/backfill_discovery_index.py.
const axes = [
  {
    axis: "Domain disruption",
    weight: 30,
    formula: "1 − 0.4^(lost + 0.5·gained)",
    desc: "Loss of a conserved Pfam domain is the strongest signal of functional disruption. The term saturates (1 lost → 0.60, 2 → 0.84), so an isoform that sheds many domains can't dominate the ranking.",
  },
  {
    axis: "Length change",
    weight: 20,
    formula: "min(|Δlen| / 150, 1)",
    desc: "A proxy for gross structural change, and the only structural signal available for novel isoforms with no annotated Pfam domain. Saturates at a 150-aa change relative to the canonical.",
  },
  {
    axis: "Neoantigen priority",
    weight: 20,
    formula: "P1→1 · P2→0.75 · P3→0.5 · P4→0.25",
    desc: "The best neoepitope priority tier per isoform (presentation consensus × immunogenicity) — a P1 hit contributes full weight, a P4 a quarter.",
  },
  {
    axis: "NMD fate",
    weight: 15,
    formula: "target→1, else→0",
    desc: "Transcripts predicted to be degraded by nonsense-mediated decay likely never yield a stable protein — a strong regulatory consequence.",
  },
  {
    axis: "PTM accessibility",
    weight: 15,
    formula: "lost→1 · exposed/buried→0.6",
    desc: "Deleting a post-translational-modification site outweighs merely repositioning it, where the residue is retained but its solvent accessibility changes.",
  },
];

const impactProfiles = [
  { label: "All", desc: "No filter — every isoform, ranked by impact score." },
  { label: "High Confidence", desc: "Only isoforms with mean pLDDT ≥ 70, where the structural prediction is reliable." },
  { label: "NMD Targets", desc: "Transcripts flagged as NMD targets — likely never translated." },
  { label: "Domain Loss", desc: "At least one Pfam domain lost relative to the canonical protein." },
  { label: "Neoantigens", desc: "Carrying at least one ranked, presented neoepitope (priority tier P1–P4)." },
  { label: "PTM Disrupted", desc: "Deleting or altering the accessibility of a PTM site." },
];

const interpretation = [
  { band: "< 20", label: "Low", desc: "canonical or near-canonical isoforms cluster here." },
  { band: "20–40", label: "Moderate", desc: "a single clear perturbation (one lost domain, or NMD targeting)." },
  { band: "40–60", label: "High", desc: "several converging lines of evidence — strong candidates for follow-up." },
  { band: "≥ 60", label: "Severe", desc: "multiple domains lost plus NMD / neoantigen / PTM disruption on a confident structure." },
];

export default function AnalyticalDiscoveryDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4 px-3">
        <div className="flex items-center gap-3 text-indigo-500">
          <Compass className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">Discovery &amp; Browsing</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          One score to rank the whole proteome — and three lenses to explore it.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-12 text-foreground">
        <section className="space-y-4">
          <p className="leading-relaxed">
            With well over a hundred thousand isoforms in the knowledgebase, the hard part isn't
            computing consequences — it's finding the ones that matter. The{" "}
            <strong>Discovery Index</strong> pre-computes a single{" "}
            <strong>Structural Impact Score (SIS)</strong> for every isoform, so the entire proteome
            can be ranked and filtered in real time to surface the highest-impact splice events.
          </p>
        </section>

        {/* ── The Structural Impact Score ────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">The Structural Impact Score (SIS)</h2>
          <p className="leading-relaxed">
            The SIS is a <strong>normalised 0–100 composite</strong> of five biological axes. Each
            axis is mapped to a sub-score between 0 and 1 with a <em>saturating</em> function, then
            combined with fixed weights that sum to 100. Because every axis is bounded, no single
            event type can dominate the ranking, and the number is directly interpretable — an
            isoform scoring 62 carries roughly 62&nbsp;% of the maximum evidence for structural
            consequence tracked in the knowledgebase.
          </p>

          {/* Formula */}
          <div className="rounded-xl bg-indigo-500/5 border border-indigo-500/15 p-5 space-y-4 not-prose">
            <div>
              <div className="text-[10px] font-bold uppercase tracking-widest text-indigo-400 mb-1.5">
                Step 1 — weighted composite (0–100)
              </div>
              <div className="font-mono text-sm text-foreground font-semibold">
                base = 30·s<sub>domain</sub> + 20·s<sub>length</sub> + 20·s<sub>neo</sub> + 15·s<sub>nmd</sub> + 15·s<sub>ptm</sub>
              </div>
              <p className="text-xs text-muted-foreground mt-1.5">
                Each s<sub>i</sub> ∈ [0, 1] (see the axes below); the weights sum to 100, so{" "}
                <span className="font-mono">base</span> is bounded 0–100.
              </p>
            </div>
            <div className="border-t border-indigo-500/15 pt-4">
              <div className="text-[10px] font-bold uppercase tracking-widest text-indigo-400 mb-1.5">
                Step 2 — confidence temper
              </div>
              <div className="font-mono text-sm text-foreground font-semibold">
                SIS = base × (0.7 + 0.3 · pLDDT / 100)
              </div>
              <p className="text-xs text-muted-foreground mt-1.5">
                The model's mean pLDDT scales the score by a factor of 0.70–1.00: low-confidence
                models are discounted but never zeroed, so a real domain loss on a lower-confidence
                novel isoform stays visible while genuine noise is down-weighted. Reference isoforms
                without an AlphaFold model are given a neutral pLDDT of 70 for this step.
              </p>
            </div>
          </div>

          {/* Axes table */}
          <div className="overflow-x-auto rounded-xl border border-border/40">
            <table className="w-full text-sm">
              <thead className="text-left">
                <tr className="border-b border-border bg-muted/40">
                  <th className="px-4 py-3 whitespace-nowrap">Axis</th>
                  <th className="px-4 py-3 whitespace-nowrap text-right">Weight</th>
                  <th className="px-4 py-3 whitespace-nowrap">Sub-score (0–1)</th>
                  <th className="px-4 py-3">What it captures</th>
                </tr>
              </thead>
              <tbody>
                {axes.map((a) => (
                  <tr key={a.axis} className="border-b border-border/40 align-top">
                    <td className="px-4 py-3 font-semibold text-foreground whitespace-nowrap">{a.axis}</td>
                    <td className="px-4 py-3 text-right font-mono text-muted-foreground">{a.weight}</td>
                    <td className="px-4 py-3 font-mono text-xs text-muted-foreground whitespace-nowrap">{a.formula}</td>
                    <td className="px-4 py-3 text-muted-foreground leading-relaxed">{a.desc}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {/* Interpretation */}
          <div className="rounded-xl bg-muted/20 border border-border/40 p-4 space-y-2 not-prose">
            <div className="text-sm font-bold text-foreground">Reading the score</div>
            {interpretation.map((b) => (
              <div key={b.band} className="flex gap-3 text-sm text-muted-foreground">
                <span className="font-mono text-foreground w-16 shrink-0">{b.band}</span>
                <span className="font-semibold text-foreground/80 w-20 shrink-0">{b.label}</span>
                <span className="leading-relaxed">{b.desc}</span>
              </div>
            ))}
          </div>
        </section>

        {/* ── Exploring the knowledgebase (three lenses) ─────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Exploring the Knowledgebase</h2>
          <p className="leading-relaxed">
            Three complementary lenses open the same data from different angles, each reachable from
            the site header. They share the underlying models but rank and filter for a different
            question.
          </p>

          <div className="grid gap-4 md:grid-cols-3 not-prose">
            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 flex flex-col gap-2">
              <div className="font-bold text-foreground">Impact Explorer</div>
              <p className="text-xs text-muted-foreground leading-relaxed flex-1">
                Ranks every isoform proteome-wide by its Structural Impact Score. Filter by impact
                profile and structural confidence, then open any row for the full isoform page.
              </p>
              <Link href="/explore/impact" className="text-xs font-semibold text-primary hover:underline">
                Open the explorer →
              </Link>
            </div>
            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 flex flex-col gap-2">
              <div className="font-bold text-foreground">Neoepitope Atlas</div>
              <p className="text-xs text-muted-foreground leading-relaxed flex-1">
                A proteome-wide table of splice-junction neoantigens, ranked into P1–P4 priority
                tiers with filters for HLA allele, binding, and junction type.
              </p>
              <div className="text-xs font-semibold flex gap-3">
                <Link href="/explore/neoepitopes" className="text-primary hover:underline">Open the atlas →</Link>
                <Link href="/documentation/neoantigens" className="text-muted-foreground hover:text-primary hover:underline">Guide</Link>
              </div>
            </div>
            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 flex flex-col gap-2">
              <div className="font-bold text-foreground">Drug–Isoform Atlas</div>
              <p className="text-xs text-muted-foreground leading-relaxed flex-1">
                Ranks drugged genes by how much their splicing subtracts from each drug's binding
                pocket — preserved, disrupted, or lost.
              </p>
              <div className="text-xs font-semibold flex gap-3">
                <Link href="/explore/drugs" className="text-primary hover:underline">Open the atlas →</Link>
                <Link href="/documentation/drug-interactions" className="text-muted-foreground hover:text-primary hover:underline">Guide</Link>
              </div>
            </div>
          </div>

          <h3 className="text-lg font-bold text-foreground pt-2">Impact profiles</h3>
          <p className="leading-relaxed">
            The Impact Explorer ships with preset filters that answer common questions in one click,
            on top of the free sliders for structural confidence and isoforms-per-gene:
          </p>
          <div className="rounded-xl border border-border/40 divide-y divide-border/40 not-prose">
            {impactProfiles.map((p) => (
              <div key={p.label} className="flex flex-col sm:flex-row sm:gap-4 px-4 py-2.5">
                <span className="font-semibold text-foreground sm:w-40 shrink-0">{p.label}</span>
                <span className="text-sm text-muted-foreground leading-relaxed">{p.desc}</span>
              </div>
            ))}
          </div>
        </section>

        {/* Footer nav */}
        <div className="pt-10 border-t border-border/20 flex items-center justify-between gap-4">
          <Link href="/documentation/drug-interactions" className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2">
            <span className="opacity-60">&larr;</span> Prev: Drug Interactions
          </Link>
          <div className="flex flex-col items-end gap-1">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest text-right">Back to Start</span>
            <Link href="/documentation" className="text-lg font-bold text-primary hover:underline flex items-center gap-2">
              Overview <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
