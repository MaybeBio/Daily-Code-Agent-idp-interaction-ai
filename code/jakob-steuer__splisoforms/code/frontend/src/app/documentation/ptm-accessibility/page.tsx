"use client";

import {
  Activity,
  ExternalLink,
  HelpCircle,
  Maximize,
  Target,
  Thermometer,
} from "lucide-react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// Status definitions kept in sync with
// backend/app/analysis/ptm_accessibility.py (after the e4b8c1a7d9f2 redesign).
const classifications = [
  {
    name: "Intact",
    color: "text-emerald-400",
    bg: "bg-emerald-400/10",
    rule: "|Δ SASA| < 10 Å² OR |Δ rel| < 25%",
    desc: "Site retained with no meaningful accessibility change.",
  },
  {
    name: "Exposed",
    color: "text-rose-400",
    bg: "bg-rose-400/10",
    rule: "Δ SASA > 10 Å² AND Δ rel > +25%",
    desc: "Site moved into a more solvent-accessible state — both the absolute and relative thresholds must clear.",
  },
  {
    name: "Buried",
    color: "text-blue-400",
    bg: "bg-blue-400/10",
    rule: "Δ SASA < -10 Å² AND Δ rel < -25%",
    desc: "Site moved into a more solvent-shielded state.  A conformational reshuffle or domain insertion.",
  },
  {
    name: "Lost",
    color: "text-slate-400",
    bg: "bg-slate-400/10",
    rule: "Canonical residue deleted",
    desc: "Splice-aware alignment maps the canonical position to a gap in the novel — the residue is not in the protein at all.",
  },
  {
    name: "Undetermined",
    color: "text-amber-400",
    bg: "bg-amber-400/10",
    rule: "Missing SASA OR min pLDDT < 70",
    desc: "SASA could not be computed for one or both structures, OR the residue lies in a low-confidence region of the AF prediction. Side-chain placement is unreliable here — we deliberately refuse to call exposed/buried. Every site on a canonical isoform compared against itself is also undetermined.",
  },
];

// Row counts in the October 2026 release (splisoforms_v2): one row per site × isoform.
const ptmCounts = [
  { label: "site × isoform rows", value: "2,087,743", cls: "text-foreground" },
  { label: "isoforms with sites", value: "107,603", cls: "text-foreground" },
  { label: "lost (residue spliced out)", value: "920,368", cls: "text-slate-400" },
  { label: "undetermined", value: "670,835", cls: "text-amber-400" },
  { label: "intact", value: "435,689", cls: "text-emerald-400" },
  { label: "exposed + buried", value: "60,851", cls: "text-rose-400" },
];

export default function RegulatoryImpactsDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4 px-3">
        <div className="flex items-center gap-3 text-rose-500">
          <Target className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">PTM Accessibility</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          Quantifying how splicing rewires PTM accessibility — with confidence masking and a
          relative-change threshold so calls reflect genuine conformational change, not noise.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-12 text-foreground">

        <section className="space-y-4">
          <p className="leading-relaxed">
            Post-translational modifications (PTMs) are chemical marks that switch a protein's
            activity, stability, or localisation on and off. A modification can only occur if its
            target residue is physically reachable — so when splicing reshapes the fold around a PTM
            site, it can quietly enable or disable that regulatory switch. For every known PTM site on
            the canonical isoform, SPLISOFORMS projects the residue onto the novel isoform through a{" "}
            <strong>splice-aware alignment</strong> (the same k-mer-anchored algorithm used by the MSA
            panel) and recomputes its solvent-accessible surface area (SASA) with{" "}
            <a href="https://freesasa.github.io/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">FreeSASA</a>{" "}
            on both AF3 models.
          </p>
          <p className="leading-relaxed">
            Sites come from two sources keyed on the gene's reviewed{" "}
            <a href="https://www.uniprot.org/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">UniProt</a>{" "}
            accession: the curated <em>Modified residue</em> features of the UniProt entry, and the
            large-scale mass-spectrometry sites served by the{" "}
            <a href="https://www.ebi.ac.uk/proteins/api/doc/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">EBI Proteins API</a>{" "}
            (proteomics PTM set). The proteomics set contributes about four sites in five, and all
            ubiquitination sites. We track the three best-annotated regulatory modifications:{" "}
            <strong>phosphorylation</strong> (on Ser/Thr/Tyr), <strong>ubiquitination</strong>, and{" "}
            <strong>acetylation</strong> (both on Lys) — these dominate the annotated regulatory
            landscape and have the coverage needed for reliable comparison across isoforms.
          </p>
          <div className="rounded-xl border border-rose-500/20 bg-rose-500/5 p-4 text-sm text-muted-foreground leading-relaxed">
            <strong className="text-foreground">Coordinate-frame guard.</strong> A site is only used
            when the residue letter at its annotated position matches our canonical sequence. Sites
            annotated on another sequence version or isoform fail this check and are dropped — about
            11 % of the raw sites. Before October 2026 those sites were measured on the wrong residue.
          </div>
        </section>

        <section className="space-y-8">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Thermometer className="w-5 h-5 text-rose-500" />
            Five-way classification
          </h2>
          <p>
            A site is reclassified <em>only</em> when both the absolute Δ
            (&gt; 10 Å²) and the relative Δ (&gt; 25% of the canonical SASA)
            clear their thresholds, and when the minimum pLDDT across the two
            structures is &ge; 70.  Otherwise the call stays at "intact" or
            "undetermined" — we'd rather be honest about uncertainty than
            inflate the count of regulatory rewirings.
          </p>

          <div className="overflow-x-auto rounded-xl border border-border/40">
            <table className="w-full text-sm">
              <thead className="text-left">
                <tr className="border-b border-border bg-muted/40">
                  <th className="px-4 py-3 whitespace-nowrap">Class</th>
                  <th className="px-4 py-3 whitespace-nowrap">Condition</th>
                  <th className="px-4 py-3">What it means</th>
                </tr>
              </thead>
              <tbody>
                {classifications.map((c) => (
                  <tr key={c.name} className="border-b border-border/40 align-top">
                    <td className="px-4 py-3">
                      <Badge className={`${c.bg} ${c.color} hover:${c.bg} text-[10px] uppercase font-bold tracking-widest px-2 py-0.5 ring-1 ring-border/50 whitespace-nowrap`}>
                        {c.name}
                      </Badge>
                    </td>
                    <td className="px-4 py-3 font-mono text-xs text-muted-foreground whitespace-nowrap">{c.rule}</td>
                    <td className="px-4 py-3 text-muted-foreground leading-relaxed">{c.desc}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          <div className="rounded-xl border border-border/40 bg-muted/10 p-4 not-prose">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              PTM accessibility · October 2026 release
            </p>
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
              {ptmCounts.map((c) => (
                <div key={c.label}>
                  <div className={`text-xl font-bold font-mono ${c.cls}`}>{c.value}</div>
                  <div className="text-xs text-muted-foreground">{c.label}</div>
                </div>
              ))}
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              Read <em>lost</em> and <em>exposed / buried</em> separately: most PTM "alterations" are
              the modified residue being spliced out, not a change in the accessibility of a retained
              residue. Counts fell from the July 2026 release because of the coordinate-frame guard.
            </p>
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <HelpCircle className="w-5 h-5 text-amber-500" />
            Why a relative threshold matters
          </h2>
          <p className="leading-relaxed">
            Baseline SASA varies an order of magnitude across residues — a
            fully-buried residue can have ~5 Å², a fully-exposed surface
            residue ~150 Å².  An absolute ±10 Å² rule labels the same physical
            shift as "huge" in one case and "noise" in the other.  Requiring a
            simultaneous ±25% relative change calibrates the call to the
            residue's own baseline.  The combined rule was chosen because:
          </p>
          <ul className="list-disc list-outside pl-6 space-y-1 text-foreground/85 text-sm">
            <li>
              Either threshold alone over-calls at the extremes — a relative-only rule inflates
              changes on tiny buried residues, an absolute-only rule inflates them on large exposed
              ones. Requiring both keeps only shifts that are real in absolute <em>and</em> relative
              terms.
            </li>
            <li>
              <em>Undetermined</em> is kept as its own class, separate from <em>intact</em>, so that
              "we couldn't measure this reliably" is never silently reported as "no change".
            </li>
          </ul>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Activity className="w-5 h-5 text-rose-500" />
            Reading the lollipop map
          </h2>
          <p className="leading-relaxed">
            On each isoform page, PTM sites are drawn as a <strong>lollipop map</strong> along the
            protein sequence (aa&nbsp;1 on the left to the C-terminus on the right), in two stacked
            tracks:
          </p>
          <ul className="list-disc list-outside pl-6 space-y-2 text-foreground/85">
            <li>
              <strong>Canonical reference</strong> — every known site, coloured by modification type
              (phosphorylation, ubiquitination, acetylation).
            </li>
            <li>
              <strong>Isoform accessibility</strong> — the same sites at their canonical positions,
              coloured by their accessibility call (exposed, buried, or unchanged); sites spliced out
              of the isoform appear as an <strong>✕</strong>. A coverage ribbon above the track marks
              which stretches of the canonical protein are present, deleted, or outside the isoform.
            </li>
          </ul>
          <div className="rounded-xl border border-amber-500/20 bg-amber-500/5 p-4 text-sm text-muted-foreground leading-relaxed">
            <strong className="text-foreground">Stem height carries no meaning.</strong> The varying
            lollipop heights only stagger sites that sit close together so they don't overlap — height
            is not a magnitude. Read a site's <em>colour</em> for its call, and hover it for the exact
            ΔSASA.
          </div>
          <div className="grid gap-3 sm:grid-cols-3">
            <div className="p-4 rounded-lg border border-emerald-500/20 bg-emerald-500/5">
              <h4 className="font-bold text-emerald-400 text-sm mb-1">Newly exposed</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                A phosphosite buried in the canonical fold ends up on the surface of the isoform — a
                kinase that couldn't reach it now can.
              </p>
            </div>
            <div className="p-4 rounded-lg border border-red-500/20 bg-red-500/5">
              <h4 className="font-bold text-red-400 text-sm mb-1">Newly buried</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                A site exposed in the canonical is shielded in the isoform — the modification is
                effectively blocked.
              </p>
            </div>
            <div className="p-4 rounded-lg border border-orange-500/20 bg-orange-500/5">
              <h4 className="font-bold text-orange-400 text-sm mb-1">Lost (✕)</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                The modified residue is spliced out of the isoform entirely — the regulatory site is
                gone.
              </p>
            </div>
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Maximize className="w-5 h-5 text-indigo-500" />
            Structural context (TED)
          </h2>
          <p>
            We also tag each site by its position relative to the{" "}
            <a href="https://ted.cathdb.info/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">TED</a>{" "}
            structural-domain boundaries. This is contextual only — it does <em>not</em> change the
            accessibility classification, and it is drawn purely from the structure-based TED domains,
            not from Pfam sequence domains.
          </p>
          <div className="grid gap-6 sm:grid-cols-3">
            <div className="p-4 rounded-lg bg-indigo-500/5 border border-indigo-500/10">
              <h4 className="font-bold text-indigo-400 text-sm mb-1 uppercase tracking-tight">
                Intra-domain
              </h4>
              <p className="text-xs text-slate-400">
                The site sits within a stable, folded domain in the novel isoform.
              </p>
            </div>
            <div className="p-4 rounded-lg bg-rose-500/5 border border-rose-500/10">
              <h4 className="font-bold text-rose-400 text-sm mb-1 uppercase tracking-tight">
                Disrupted
              </h4>
              <p className="text-xs text-slate-400">
                The site was inside a canonical domain that splicing fractured —
                no longer in a domain in the novel.
              </p>
            </div>
            <div className="p-4 rounded-lg bg-slate-500/5 border border-slate-500/10">
              <h4 className="font-bold text-slate-400 text-sm mb-1 uppercase tracking-tight">
                Inter-domain
              </h4>
              <p className="text-xs text-slate-400">
                The site is in a flexible linker / loop region in both isoforms,
                OR no TED domains were annotated.
              </p>
            </div>
          </div>
        </section>



        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "UniProt", note: "curated Modified-residue features", href: "https://www.uniprot.org/" },
              { label: "EBI Proteins API", note: "large-scale proteomics PTM sites", href: "https://www.ebi.ac.uk/proteins/api/doc/" },
              { label: "FreeSASA", note: "solvent-accessible surface area", href: "https://freesasa.github.io/" },
              { label: "TED", note: "structural-domain context", href: "https://ted.cathdb.info/" },
            ].map((s) => (
              <li key={s.label}>
                <a href={s.href} target="_blank" rel="noopener noreferrer"
                   className="flex items-center gap-1.5 text-sm text-muted-foreground hover:text-primary transition-colors">
                  <ExternalLink className="w-3.5 h-3.5 shrink-0" />
                  <span className="font-semibold text-foreground">{s.label}</span>
                  <span className="text-xs">— {s.note}</span>
                </a>
              </li>
            ))}
          </ul>
        </section>

        <div className="pt-10 border-t border-border/20 flex items-center justify-between gap-4">
          <Link
            href="/documentation/neoantigens"
            className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2"
          >
            <span className="opacity-60">&larr;</span> Prev: Neoantigens
          </Link>
          <div className="flex flex-col items-end gap-1">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest text-right">
              Next Section
            </span>
            <Link
              href="/documentation/interactions-expression"
              className="text-lg font-bold text-primary hover:underline flex items-center gap-2"
            >
              Interactions &amp; Expression <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
