"use client";

import { Puzzle, ExternalLink, Hammer, ShieldQuestion, Tags } from "lucide-react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// ── Alteration geometry (replaces the legacy four-status mapping) ───────────
const alterations = [
  {
    name: "Intact",
    color: "bg-emerald-500/10 text-emerald-500",
    desc: "Full-coverage HMM hit (>= 80% of the model length). The domain looks complete in this isoform.",
  },
  {
    name: "N-terminal truncation",
    color: "bg-amber-500/10 text-amber-500",
    desc: "The hit starts well into the domain model, so the N-terminal part of the canonical domain is missing. Depending on the domain this can strip an auto-inhibitory region or remove catalytic machinery (e.g. the kinase N-lobe) — the effect depends on exactly which residues are gone.",
  },
  {
    name: "C-terminal truncation",
    color: "bg-amber-500/10 text-amber-500",
    desc: "The hit ends well before the domain model does, so the C-terminal part is missing. This can drop regulatory tails, docking sites or substrate-recognition loops while sparing the catalytic core — or remove the core itself.",
  },
  {
    name: "Both termini truncated",
    color: "bg-amber-500/10 text-amber-500",
    desc: "An internal slice of the HMM model is retained while both ends are missing. Often the result of an exon-skipping event that lands inside the domain.",
  },
  {
    name: "Internal disruption",
    color: "bg-orange-500/10 text-orange-500",
    desc: "Both ends are inside the HMM model and HMM coverage is below 50% — consistent with a small fragment surviving inside an otherwise lost domain, or an internal in-frame deletion.",
  },
  {
    name: "Lost",
    color: "bg-rose-500/10 text-rose-500",
    desc: "The domain is present in the canonical isoform but entirely absent from the novel one. It is recorded with the canonical's coordinates and shown only on the canonical reference track — a lost row is a domain the isoform does not have, not a domain instance on the isoform.",
  },
  {
    name: "Gained",
    color: "bg-blue-500/10 text-blue-500",
    desc: "The domain matches in the novel isoform but not the canonical — usually a novel exon or a frameshift that exposes a different reading frame.",
  },
];

// ── Curated functional classes (display side; source of truth is
// backend/app/analysis/functional_classes.py) ─────────────────────────────────
const classes = [
  { name: "Kinase",                     family: "catalytic", examples: "Pkinase (PF00069), Pkinase_Tyr (PF07714)" },
  { name: "Phosphatase",                family: "catalytic", examples: "Y_phosphatase (PF00102), DSPc (PF00782), PP2C (PF00481)" },
  { name: "Protease",                   family: "catalytic", examples: "Trypsin (PF00089), Peptidase_C1 (PF00112), Caspase (PF00656)" },
  { name: "Catalytic (other)",          family: "catalytic", examples: "DEAD helicase (PF00270), Rossmann-fold dehydrogenases, P-loop NTPases" },
  { name: "Ubiquitin ligase / E2",      family: "catalytic", examples: "RING (PF00097), HECT (PF00632), UQ_con / E2 (PF00179)" },
  { name: "DNA-binding",                family: "binding",   examples: "Homeobox (PF00046), HLH (PF00010), bZIP (PF00170), C2H2 zinc finger (PF00096)" },
  { name: "RNA-binding",                family: "binding",   examples: "RRM (PF00076), KH (PF00013), dsrm (PF00035)" },
  { name: "Signaling adapter",          family: "adapter",   examples: "SH2 (PF00017), SH3 (PF00018), PDZ (PF00595), PH (PF00169), WW (PF00397)" },
  { name: "Oligomerization / scaffold", family: "scaffold",  examples: "Ankyrin (PF00023), WD40 (PF00400), TPR clan (CL0020)" },
  { name: "Other / unclassified",       family: "other",     examples: "Anything outside the curated list — explicitly NOT a low-priority bucket; many TFs and RNA-processing factors fall through here when their Pfam family lacks a clan assignment." },
];

const FAMILY_STYLES: Record<string, { badge: string; row: string }> = {
  catalytic: { badge: "bg-violet-500/15 text-violet-600 dark:text-violet-300 border-violet-500/30", row: "bg-violet-500/5" },
  binding:   { badge: "bg-sky-500/15    text-sky-600    dark:text-sky-300    border-sky-500/30",    row: "bg-sky-500/5" },
  adapter:   { badge: "bg-teal-500/15   text-teal-600   dark:text-teal-300   border-teal-500/30",   row: "bg-teal-500/5" },
  scaffold:  { badge: "bg-indigo-500/15 text-indigo-600 dark:text-indigo-300 border-indigo-500/30", row: "bg-indigo-500/5" },
  other:     { badge: "bg-muted-foreground/10 text-muted-foreground border-border", row: "" },
};

export default function FunctionalDomainsDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4">
        <div className="flex items-center gap-3 text-emerald-500">
          <Puzzle className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">Functional Domains</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          Two complementary domain definitions — sequence-based (Pfam) and structure-based (TED) —
          mapped onto splicing events without prejudging the outcome.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-10 text-foreground">

        <section className="space-y-4">
          <p className="leading-relaxed">
            SPLISOFORMS scans every isoform's protein sequence against the full{" "}
            <a href="https://www.ebi.ac.uk/interpro/entry/pfam/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">Pfam-A</a>{" "}
            HMM library with{" "}
            <a href="http://hmmer.org/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">HMMER</a>{" "}
            (via <a href="https://github.com/althonos/pyhmmer" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">pyhmmer</a>),
            recording every significant domain hit with residue-level coordinates. Each hit is
            then annotated along two independent axes:
          </p>
          <ul className="list-disc list-inside space-y-2 text-foreground/85">
            <li>
              <strong>Functional class</strong> — what the domain <em>does</em>,
              based on a curated Pfam clan + accession mapping.
            </li>
            <li>
              <strong>Alteration type</strong> — <em>where</em> the splicing event
              has perturbed the domain (N-terminal, C-terminal, both, internal,
              or full loss / gain).
            </li>
          </ul>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Hammer className="w-5 h-5 text-emerald-500" />
            The Mapping Pipeline
          </h2>
          <p className="leading-relaxed">
            The key idea is that an isoform and its canonical reference go through the{" "}
            <em>same</em> scan — so any difference between them reflects the splice event, not a
            change in method. Mapping runs in three steps:
          </p>
          <ol className="space-y-3 leading-relaxed list-decimal list-outside pl-6 marker:font-bold marker:text-emerald-500">
            <li>
              <strong>Scan.</strong> Each isoform's protein sequence is searched against Pfam-A
              with pyhmmer using Pfam's own gathering (GA) score thresholds. Families in the same
              Pfam <em>clan</em> are related models that often match the same stretch of sequence;
              where same-clan hits overlap, only the best-scoring one (lowest E-value) is kept — the
              <code className="font-mono text-xs mx-1">pfam_scan</code> convention — so a zinc-finger
              array is counted once per finger, not once per related model. Every surviving hit is
              recorded with 1-based residue coordinates and how much of the domain model it covers,
              so a partial hit can be pinned to a specific part of the domain.
            </li>
            <li>
              <strong>Compare.</strong> The same scan is run on the gene's designated canonical
              isoform (the MANE Select transcript where one exists). Matching each hit to the
              canonical hit for the same Pfam family — and comparing how much of the model each one
              covers — is what determines the <em>alteration type</em>.
            </li>
            <li>
              <strong>Classify.</strong> Every hit is also tagged with a <em>functional class</em>
              {" "}derived from its Pfam clan and accession, describing what the domain does.
            </li>
          </ol>

          <div className="rounded-xl border border-border/40 bg-muted/10 p-4 not-prose">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              Pfam · October 2026 release
            </p>
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
              <div>
                <div className="text-xl font-bold font-mono text-emerald-400">216,364</div>
                <div className="text-xs text-muted-foreground">domain hits on isoforms</div>
              </div>
              <div>
                <div className="text-xl font-bold font-mono text-foreground">104,135</div>
                <div className="text-xs text-muted-foreground">isoforms with ≥ 1 hit</div>
              </div>
              <div>
                <div className="text-xl font-bold font-mono text-rose-400">107,855</div>
                <div className="text-xs text-muted-foreground">canonical domains recorded as lost</div>
              </div>
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              The table holds 324,219 rows in total; the lost rows are placeholders on the canonical,
              so they are not counted as domain instances. Clan de-duplication removed roughly 70,000
              redundant overlapping hits relative to the July 2026 release.
            </p>
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Tags className="w-5 h-5 text-violet-500" />
            Functional Classes
          </h2>
          <p className="leading-relaxed">
            The classes below carve up the Pfam universe into broad
            <em> activity buckets</em> that surface together in the UI and the
            API.  The taxonomy is intentionally coarse and conservative —
            curated assignments only.  Domains whose Pfam family has no clan
            assignment, or whose clan we haven't reviewed yet, land in
            <strong> Other / unclassified</strong>.  That bucket is not "low
            importance" — it's <em>unsorted</em>, and many transcription
            co-factors, RNA-processing factors, and structural adaptors live
            there.  Treat the class as a filter / facet, not a verdict.
          </p>

          <div className="overflow-x-auto rounded-xl border border-border/40">
            <table className="w-full text-sm">
              <thead className="text-left">
                <tr className="border-b border-border bg-muted/40">
                  <th className="px-4 py-3">Class</th>
                  <th className="px-4 py-3">Family</th>
                  <th className="px-4 py-3">Example Pfam IDs</th>
                </tr>
              </thead>
              <tbody>
                {classes.map((c) => {
                  const fs = FAMILY_STYLES[c.family];
                  return (
                    <tr key={c.name} className={`border-b border-border/40 ${fs.row}`}>
                      <td className="px-4 py-3 font-medium text-foreground">
                        <Badge className={`${fs.badge} text-[10px] uppercase tracking-widest`}>
                          {c.name}
                        </Badge>
                      </td>
                      <td className="px-4 py-3 text-muted-foreground capitalize text-xs">
                        {c.family}
                      </td>
                      <td className="px-4 py-3 text-foreground/80 font-mono text-xs">{c.examples}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Alteration Geometry</h2>
          <p className="leading-relaxed">
            Each domain hit carries an <strong>alteration_type</strong> describing
            where the splicing event has changed the domain relative to the
            canonical isoform's hit on the same Pfam family.  This matters
            because the <em>same</em> family can be activated, deactivated, or
            converted into a dominant-negative version depending on which
            piece of the model is removed.
          </p>

          <div className="grid gap-3">
            {alterations.map((a) => (
              <div
                key={a.name}
                className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start"
              >
                <Badge
                  className={`${a.color} hover:${a.color} text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap`}
                >
                  {a.name}
                </Badge>
                <p className="text-sm leading-relaxed text-muted-foreground font-medium">
                  {a.desc}
                </p>
              </div>
            ))}
          </div>

          <div className="p-4 rounded-xl bg-muted/20 border border-border/40 text-xs text-muted-foreground space-y-1">
            <p className="font-bold text-foreground text-sm mb-2">How the cut-offs work</p>
            <p className="leading-relaxed">
              A domain counts as <strong>intact</strong> when the isoform covers at least{" "}
              <strong>80%</strong> of the domain model (or at least <strong>85%</strong> of what
              the canonical isoform covers for the same family). Whether a truncation is called
              N- or C-terminal depends on which end is missing: a hit that only starts after the
              first <strong>10%</strong> of the model has lost its N-terminus, and one that ends
              before the last <strong>10%</strong> has lost its C-terminus.
            </p>
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">TED — Structural Domains</h2>
          <p className="leading-relaxed">
            While Pfam answers "which sequence motifs are present",{" "}
            <a href="https://ted.cathdb.info/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">TED (The Encyclopedia of Domains)</a>{" "}
            answers the complementary, structure-based question: how does the folded protein partition into
            independent, compact domain units?
          </p>
          <p className="leading-relaxed">
            Two kinds of structural domain are stored, and the isoform page labels them differently:
          </p>
          <div className="grid gap-4 sm:grid-cols-2">
            <div className="p-4 rounded-lg bg-muted/20 border border-border/40 space-y-2">
              <h4 className="font-bold text-foreground text-sm">TED consensus domains (canonical)</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                TED's published human set (259,068 domains, a three-method consensus of Chainsaw,
                Merizo and UniDoc on AlphaFold Database models, with CATH labels) is mapped onto each
                gene's designated canonical isoform through the gene's UniProt accession. Only high-
                and medium-consensus domains are kept, and a domain that would extend past the end of
                our canonical protein is skipped rather than clipped.
              </p>
              <p className="text-xs text-foreground/80 font-mono">42,599 domains · 21,154 isoforms</p>
            </div>
            <div className="p-4 rounded-lg bg-muted/20 border border-border/40 space-y-2">
              <h4 className="font-bold text-foreground text-sm">Chainsaw de novo segmentation (all isoforms)</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Because every isoform has its own AlphaFold 3 structure, canonical and novel isoforms
                are all segmented <i>de novo</i> with{" "}
                <a href="https://github.com/JudeWells/chainsaw" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">Chainsaw</a>.
                Boundaries are drawn purely from the 3D fold, so a novel isoform is evaluated on the
                same footing as the reference. Chainsaw reports no fold classification, only
                boundaries; these are stored as 1-based residue positions like every other domain.
              </p>
              <p className="text-xs text-foreground/80 font-mono">193,290 domains · 101,830 isoforms</p>
            </div>
          </div>
          <div className="rounded-xl bg-muted/20 border border-border/40 p-4 text-xs text-muted-foreground mt-4">
            <p className="font-bold text-foreground text-sm mb-1">Why two systems?</p>
            <p className="leading-relaxed">
              A splice event can leave a Pfam motif intact while collapsing the structural domain that hosts
              it (or vice-versa). Reading the two side by side separates "the sequence signature survives" from
              "the fold survives".
            </p>
          </div>
        </section>



        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "Pfam-A", note: "sequence domain families (InterPro)", href: "https://www.ebi.ac.uk/interpro/entry/pfam/" },
              { label: "HMMER", note: "profile-HMM domain search", href: "http://hmmer.org/" },
              { label: "pyhmmer", note: "HMMER Python bindings", href: "https://github.com/althonos/pyhmmer" },
              { label: "TED", note: "The Encyclopedia of Domains", href: "https://ted.cathdb.info/" },
              { label: "Chainsaw", note: "de novo structural domain segmentation", href: "https://github.com/JudeWells/chainsaw" },
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
            href="/documentation/protein-structure"
            className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2"
          >
            <span className="opacity-60">&larr;</span> Prev: Protein Structure
          </Link>
          <div className="flex flex-col items-end gap-1">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest text-right">
              Next Section
            </span>
            <Link
              href="/documentation/nmd-transcript-fate"
              className="text-lg font-bold text-primary hover:underline flex items-center gap-2"
            >
              NMD &amp; Transcript Fate <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
