"use client";

import { AlertTriangle, ExternalLink, Shield } from "lucide-react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// Ranking v2 — priority tier is the conjunction of presentation consensus × immunogenicity
const priorityTiers = [
  { name: "P1 · High-confidence",     color: "text-emerald-400",      bg: "bg-emerald-500/15", desc: "Strong multi-tool presentation (≥ 2 of 3 tools call the peptide a strong binder) AND predicted immunogenic (BigMHC-IM ≥ 0.5). The candidates most likely to be presented and elicit a T-cell response." },
  { name: "P2 · Strong presentation", color: "text-sky-400",          bg: "bg-sky-500/10",     desc: "Strong multi-tool presentation, but not predicted immunogenic. Robustly presented; immunogenicity support is absent." },
  { name: "P3 · Supported",           color: "text-amber-400",        bg: "bg-amber-500/10",   desc: "Moderate presentation (a single strong call, or ≥ 2 weak calls) AND predicted immunogenic." },
  { name: "P4 · Candidate",           color: "text-muted-foreground", bg: "bg-muted/30",       desc: "Presentation evidence from a single tool or without immunogenicity support. Retained as a candidate but the weakest immunological evidence." },
];

// Presentation tier (the ranking backbone) — multi-tool consensus, %rank_EL cutoffs
const presentationTiers = [
  { name: "Strong",   color: "text-emerald-400",      bg: "bg-emerald-500/15", desc: "≥ 2 of 3 tools call the peptide a strong binder: MHCflurry %rank ≤ 0.5, NetMHCpan-4.2 EL %rank ≤ 0.5, or BigMHC EL ≥ 0.5." },
  { name: "Moderate", color: "text-amber-400",        bg: "bg-amber-500/15",   desc: "Exactly one strong call, or ≥ 2 tools presenting at the weak stringency (%rank ≤ 2.0 / BigMHC EL ≥ 0.25)." },
  { name: "Weak",     color: "text-orange-300",       bg: "bg-orange-500/10",  desc: "A single tool presents the peptide at the weak stringency." },
  { name: "None",     color: "text-muted-foreground", bg: "bg-muted/30",       desc: "No tool presents the peptide — excluded from the candidate set." },
];

// Decorating flags — orthogonal annotations that describe candidates, they do NOT gate the tier
const decoratingFlags = [
  { name: "Immunogenic",   color: "text-violet-400", bg: "bg-violet-500/10", desc: "BigMHC-IM immunogenicity score ≥ 0.5 — predicted T-cell recognition. This is the immunogenicity axis of the P-tier and is also surfaced as a standalone flag." },
  { name: "Proteome-unique", color: "text-rose-400", bg: "bg-rose-500/10", desc: "The peptide occurs in no GENCODE v47 reference protein at all — not in another isoform of the same gene, not in any other gene (see Proteome uniqueness below)." },
  { name: "Splice-enabled", color: "text-orange-400", bg: "bg-orange-500/10", desc: "The epitope exists only because of the splice event (a canonical inhibitory PTM site is lost at/near the junction). The resource's unique angle — reported, not used to rank." },
  { name: "Expressed",     color: "text-emerald-400", bg: "bg-emerald-500/10", desc: "The isoform has at least one full-length read in the ccRCC per-sample expression table. Only the ccRCC set carries such counts, so this is a dataset flag, not an abundance filter." },
  { name: "PTM-disrupting", color: "text-cyan-400",  bg: "bg-cyan-500/10",   desc: "The novel sequence ablates a canonical PTM site inside the epitope (ptm_canonical_residue_lost = TRUE)." },
  { name: "High-fidelity fold", color: "text-sky-400", bg: "bg-sky-500/10",  desc: "The epitope window is confidently modelled (avg pLDDT > 70, min pLDDT > 50) — the structural context is trustworthy." },
];

// Proteome-uniqueness classes — set by backfill_neoantigen_uniqueness.py after every load.
const uniquenessClasses = [
  { chip: "UNQ",  name: "unique",     color: "text-rose-400 bg-rose-500/10 border-rose-500/25",            p1: "60,656", desc: "The peptide occurs in no GENCODE v47 reference protein. The only class that is novel to the normal proteome." },
  { chip: "SELF", name: "self_gene",  color: "text-muted-foreground bg-muted/30 border-border/40",          p1: "31,511", desc: "Also present in another annotated isoform of the same gene — a reference self-peptide the immune system is likely tolerant to." },
  { chip: "REF",  name: "other_gene", color: "text-amber-400 bg-amber-500/10 border-amber-500/25",          p1: "10,673", desc: "Present in a reference protein of a different gene." },
];

// Table-level counts in the October 2026 release (splisoforms_v2).
const neoCounts = [
  { label: "peptide × allele × isoform rows", value: "7,930,448" },
  { label: "distinct peptides", value: "1,941,536" },
  { label: "isoforms with predictions", value: "66,747" },
  { label: "P1 peptides", value: "102,547" },
  { label: "P1 isoforms", value: "42,280" },
  { label: "rows missing a predictor", value: "0" },
];

// HLA allele panel — mirrors DEFAULT_HLA_ALLELES in backend/app/analysis/neoantigen.py
const hlaAlleles = {
  A: [
    { allele: "A*01:01", freq: "~16% European" },
    { allele: "A*02:01", freq: "~28% European, ~10% Asian" },
    { allele: "A*03:01", freq: "~14% European" },
    { allele: "A*11:01", freq: "~22% Asian, ~5% European" },
    { allele: "A*23:01", freq: "~9% African" },
    { allele: "A*24:02", freq: "~18% Asian, ~5% European" },
    { allele: "A*26:01", freq: "~5% European" },
    { allele: "A*29:02", freq: "~6% African, ~4% European" },
    { allele: "A*30:01", freq: "~8% African" },
    { allele: "A*31:01", freq: "~7% Asian, ~5% European" },
    { allele: "A*32:01", freq: "~6% European" },
    { allele: "A*33:01", freq: "~8% Asian" },
    { allele: "A*68:01", freq: "~7% African, ~3% European" },
  ],
  B: [
    { allele: "B*07:02", freq: "~13% European" },
    { allele: "B*08:01", freq: "~10% European" },
    { allele: "B*13:02", freq: "~7% Asian" },
    { allele: "B*15:01", freq: "~9% European" },
    { allele: "B*15:02", freq: "~8% Asian (SJS-associated)" },
    { allele: "B*15:03", freq: "~10% sub-Saharan African" },
    { allele: "B*18:01", freq: "~8% European" },
    { allele: "B*27:05", freq: "~8% European (AS-associated)" },
    { allele: "B*35:01", freq: "~8% European / Latino" },
    { allele: "B*38:01", freq: "~4% European" },
    { allele: "B*40:01", freq: "~9% Asian" },
    { allele: "B*44:02", freq: "~10% European" },
    { allele: "B*44:03", freq: "~7% African" },
    { allele: "B*51:01", freq: "~7% Asian / Mediterranean" },
    { allele: "B*53:01", freq: "~10% sub-Saharan African" },
    { allele: "B*58:01", freq: "~8% Asian (allopurinol-HSR)" },
  ],
  C: [
    { allele: "C*03:04", freq: "~10% European / Asian" },
    { allele: "C*04:01", freq: "~12% African / European" },
    { allele: "C*06:02", freq: "~9% European" },
    { allele: "C*07:01", freq: "~18% European" },
    { allele: "C*07:02", freq: "~14% European" },
    { allele: "C*12:03", freq: "~5% European" },
  ],
};

function AlleleTable({ locus, rows, note }: { locus: string; rows: { allele: string; freq: string }[]; note?: string }) {
  const colors: Record<string, string> = {
    A: "text-sky-400 border-sky-500/30 bg-sky-500/5",
    B: "text-violet-400 border-violet-500/30 bg-violet-500/5",
    C: "text-amber-400 border-amber-500/30 bg-amber-500/5",
  };
  const headerColor: Record<string, string> = {
    A: "bg-sky-500/10 text-sky-300",
    B: "bg-violet-500/10 text-violet-300",
    C: "bg-amber-500/10 text-amber-300",
  };

  return (
    <div className={`rounded-xl border ${colors[locus]} overflow-hidden`}>
      <div className={`px-4 py-2 flex items-center justify-between ${headerColor[locus]}`}>
        <span className="font-bold text-sm tracking-widest uppercase">HLA-{locus}</span>
        <span className="text-xs font-mono opacity-70">{rows.length} alleles</span>
      </div>
      <table className="w-full text-xs">
        <thead>
          <tr className="border-b border-border/30 text-muted-foreground">
            <th className="text-left py-1.5 px-3 font-semibold">Allele</th>
            <th className="text-left py-1.5 px-3 font-semibold">Approx. population frequency</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((r, i) => (
            <tr key={r.allele} className={`border-b border-border/20 ${i % 2 === 0 ? "bg-muted/10" : ""}`}>
              <td className="py-1.5 px-3 font-mono font-medium text-foreground/90">HLA-{r.allele}</td>
              <td className="py-1.5 px-3 text-muted-foreground">{r.freq}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {note && (
        <p className="px-3 py-2 text-[10px] text-muted-foreground/70 italic border-t border-border/20">{note}</p>
      )}
    </div>
  );
}

export default function EpitopeDiscoveryDoc() {
  const total = hlaAlleles.A.length + hlaAlleles.B.length + hlaAlleles.C.length;

  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4 px-3">
        <div className="flex items-center gap-3 text-purple-500">
          <Shield className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">Neoantigens</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          Identifying novel immunological targets created by tumor-specific splicing events.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-12 text-foreground">
        <section className="space-y-4">
          <p className="leading-relaxed">
            Alternative splicing can produce protein sequences the immune system has never seen.
            Where a novel junction falls inside a short peptide that MHC-I can display, that peptide
            becomes a candidate <strong>neoantigen</strong> — a tumour-specific target for T-cell
            recognition. SPLISOFORMS extracts these junction-spanning peptides genome-wide, scores
            each for MHC-I <strong>presentation</strong> and <strong>immunogenicity</strong>, and
            ranks them into <strong>P1–P4</strong> priority tiers. A distinctive layer additionally
            flags epitopes that exist <em>only</em> because splicing removed a canonical regulatory
            (PTM) site — the resource's unique splice-enabled angle.
          </p>
        </section>

        {/* Discovery engine */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">The Discovery Engine</h2>
          <p className="leading-relaxed">
            From every novel splice junction to a ranked list of candidate epitopes, the pipeline
            runs in four stages:
          </p>
          <ol className="space-y-4 leading-relaxed list-decimal list-outside pl-6 marker:font-bold marker:text-purple-400">
            <li>
              <strong>Junction extraction.</strong> Overlapping 8–11-mer peptides are tiled across
              each novel junction — novel exon, alternative splice site, intron retention, or
              frameshift — covering every reading frame that spans the junction.
            </li>
            <li>
              <strong>Flanking context.</strong> Each peptide is extracted with 10 residues of
              upstream and downstream sequence, giving the antigen-processing model the context it
              needs to score proteasomal cleavage right at the junction.
            </li>
            <li>
              <strong>MHC-I prediction.</strong> Every peptide × allele pair is scored by a
              three-tool presentation ensemble —{" "}
              <a href="https://github.com/openvax/mhcflurry" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">MHCflurry&nbsp;2.2</a>,{" "}
              <a href="https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">NetMHCpan-4.2</a>, and{" "}
              <a href="https://github.com/KarchinLab/bigmhc" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">BigMHC</a>{" "}
              — plus BigMHC-IM for immunogenicity. MHCflurry additionally returns a 0–1 presentation
              score (proteasomal processing × MHC binding). Since the October 2026 re-run every stored
              row carries all three predictors' scores, so the vote described below is the same
              three-of-three vote everywhere in the table. How these combine into the P1–P4 ranking is
              detailed below.
            </li>
            <li>
              <strong>Candidate set.</strong> MHCflurry acts as the entry filter: a peptide × allele
              pair is stored when MHCflurry calls it a binder (IC50 ≤ 500 nM, categorised as strong
              &lt; 50 nM or weak 50–500 nM) or presented (presentation score ≥ 0.5). The other two
              tools then score every stored pair. Pairs that MHCflurry rejects outright are not in the
              table, so the candidate set is MHCflurry-anchored rather than a union of the three tools.
            </li>
          </ol>
        </section>

        {/* HLA allele panel */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">HLA Allele Panel</h2>
          <p className="text-sm text-slate-400 leading-relaxed">
            Predictions are run across a curated panel of <strong className="text-foreground/80">{total} HLA class-I alleles</strong> selected
            to maximise global population coverage across HLA-A, HLA-B, and HLA-C loci. The panel is
            derived from the{" "}
            <a href="https://www.iedb.org/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">IEDB</a>{" "}
            class-I reference set (22 of its 27 alleles are included) and adds alleles with high
            frequency in African and Asian populations plus six HLA-C alleles to reduce Eurocentric
            bias. It is a fixed panel: predictions are not patient-specific.
          </p>

          <div className="grid gap-6 lg:grid-cols-3">
            <AlleleTable locus="A" rows={hlaAlleles.A} />
            <AlleleTable locus="B" rows={hlaAlleles.B} />
            <AlleleTable
              locus="C"
              rows={hlaAlleles.C}
              note="MHCflurry HLA-C coverage is sparser than HLA-A/B. These alleles fall back to affinity-only mode when the presentation model lacks training data."
            />
          </div>

          <div className="flex flex-wrap gap-3 text-xs">
            <span className="px-2.5 py-1 rounded-full bg-sky-500/10 text-sky-400 border border-sky-500/20 font-mono">
              HLA-A — {hlaAlleles.A.length} alleles
            </span>
            <span className="px-2.5 py-1 rounded-full bg-violet-500/10 text-violet-400 border border-violet-500/20 font-mono">
              HLA-B — {hlaAlleles.B.length} alleles
            </span>
            <span className="px-2.5 py-1 rounded-full bg-amber-500/10 text-amber-400 border border-amber-500/20 font-mono">
              HLA-C — {hlaAlleles.C.length} alleles
            </span>
            <span className="px-2.5 py-1 rounded-full bg-purple-500/10 text-purple-400 border border-purple-500/20 font-mono">
              Total — {total} alleles
            </span>
          </div>
        </section>

        {/* Priority ranking (v2) */}
        <section className="space-y-8">
          <h2 className="text-2xl font-bold text-foreground">Priority Ranking (P1–P4)</h2>
          <p className="leading-relaxed">
            Candidates are ranked on <strong>immunological evidence</strong>, not on splice/PTM
            origin. The headline <strong className="text-foreground">priority tier</strong> is a
            transparent conjunction of two orthogonal axes: multi-tool MHC{" "}
            <strong className="text-foreground">presentation</strong> consensus and BigMHC{" "}
            <strong className="text-foreground">immunogenicity</strong>. Splice origin, PTM
            disruption, expression and fold quality are reported as <em>decorating flags</em> — so
            we can state, for example, <em>"of P1 epitopes, X % are splice-enabled"</em> rather than
            gating the tier on it.
          </p>

          <div className="grid gap-4">
            {priorityTiers.map((t) => (
              <div key={t.name} className="flex flex-col sm:flex-row sm:items-center gap-4 p-4 rounded-xl border border-border/40 bg-muted/10">
                <Badge className={`${t.bg} ${t.color} hover:${t.bg} text-[10px] uppercase font-bold tracking-widest px-3 py-1 sm:w-56 justify-center ring-1 ring-border/50`}>
                  {t.name}
                </Badge>
                <p className="text-sm leading-relaxed text-muted-foreground font-medium">
                  {t.desc}
                </p>
              </div>
            ))}
          </div>

          <h3 className="text-lg font-bold text-foreground pt-2">Presentation consensus (the backbone)</h3>
          <p className="leading-relaxed text-sm text-muted-foreground">
            Three predictors are run per peptide × allele pair — MHCflurry 2.x, NetMHCpan-4.2 (EL)
            and BigMHC EL — and vote at two stringencies (strong / weak) using field-standard
            eluted-ligand %rank cutoffs. The count of votes sets the presentation tier. A predictor
            that did not score a row would not vote; since October 2026 no such rows exist, so every
            tier rests on the full three-tool vote and the same peptide × allele pair carries the same
            tier in every isoform it occurs in.
          </p>
          <div className="grid gap-3">
            {presentationTiers.map((t) => (
              <div key={t.name} className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
                <Badge className={`${t.bg} ${t.color} hover:${t.bg} text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap sm:w-28 justify-center`}>
                  {t.name}
                </Badge>
                <p className="text-sm leading-relaxed text-muted-foreground font-medium">{t.desc}</p>
              </div>
            ))}
          </div>

          <div className="rounded-xl border border-border/40 bg-muted/10 p-4 not-prose">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              Neoantigen table · October 2026 release
            </p>
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
              {neoCounts.map((c) => (
                <div key={c.label}>
                  <div className="text-xl font-bold font-mono text-foreground">{c.value}</div>
                  <div className="text-xs text-muted-foreground">{c.label}</div>
                </div>
              ))}
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              In the July 2026 release 18 of 40 compute shards had run without NetMHCpan, so 40 % of
              rows rested on a two-tool vote. Those 27,593 isoforms were re-predicted with all three
              tools in October 2026 and the whole table re-tiered; P1 grew from 92,790 to 102,547
              peptides because the missing votes had erred downward.
            </p>
          </div>
        </section>

        {/* Proteome uniqueness */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Proteome uniqueness</h2>
          <p className="leading-relaxed">
            Every neoepitope is novel relative to its gene's <em>canonical</em> protein by
            construction. That alone does not make it tumour-specific: the same 8–11-mer can occur in
            another annotated isoform of the same gene, or in a different gene's protein, in which
            case it is a self-peptide. After every load, each candidate peptide is therefore looked up
            in every GENCODE v47 reference protein (the md5-deduplicated protein-coding set; long-read
            isoforms are not part of this background) and classified:
          </p>
          <div className="grid gap-3">
            {uniquenessClasses.map((u) => (
              <div key={u.chip} className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
                <span className={`shrink-0 inline-flex items-center justify-center w-14 px-2 py-0.5 rounded border text-[10px] font-bold font-mono ${u.color}`}>
                  {u.chip}
                </span>
                <div className="space-y-1">
                  <p className="text-sm leading-relaxed text-muted-foreground font-medium">
                    <span className="font-mono text-xs text-foreground/80 mr-2">{u.name}</span>{u.desc}
                  </p>
                  <p className="text-xs text-muted-foreground"><span className="font-mono text-foreground/80">{u.p1}</span> of the 102,547 P1 peptides</p>
                </div>
              </div>
            ))}
          </div>
          <p className="leading-relaxed text-sm text-muted-foreground">
            The class is shown as a <strong>UNQ / SELF / REF</strong> chip on every row of the atlas
            and the per-isoform neoantigen panel; a companion <strong>n/3</strong> chip appears only
            when fewer than three predictors scored a row (currently none). It never changes the
            priority tier — it is reported alongside, so a reader can ask "of the P1 peptides, which
            are truly absent from the normal proteome?" (59 %). In the public API,{" "}
            <code className="font-mono text-xs">/neoantigens/top?flag=proteome_unique</code> restricts
            the atlas to unique peptides and <code className="font-mono text-xs">flag=all_tools</code>{" "}
            to rows scored by all three predictors.
          </p>
        </section>

        {/* Decorating flags */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Decorating flags</h2>
          <p className="leading-relaxed">
            Orthogonal annotations layered on top of any tier. They describe <em>why</em> a
            candidate is interesting (splice origin, PTM disruption) or add supporting evidence
            (expression, fold quality) — but they never change the priority tier.
          </p>
          <div className="grid gap-3">
            {decoratingFlags.map((t) => (
              <div key={t.name} className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
                <Badge className={`${t.bg} ${t.color} hover:${t.bg} text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap sm:w-40 justify-center`}>
                  {t.name}
                </Badge>
                <p className="text-sm leading-relaxed text-muted-foreground font-medium">{t.desc}</p>
              </div>
            ))}
          </div>

          {/* Defensibility caveats */}
          <div className="flex items-start gap-3 p-4 rounded-xl border border-amber-500/25 bg-amber-500/5">
            <AlertTriangle className="w-4 h-4 text-amber-400 mt-0.5 shrink-0" />
            <p className="text-sm leading-relaxed text-muted-foreground">
              <strong className="text-amber-300">Interpretation.</strong> All tiers are{" "}
              <strong>predicted</strong>, not experimentally validated. Presentation (the P-tier
              backbone) is well benchmarked; immunogenicity relies on a single model (BigMHC-IM)
              and should be read as supporting evidence. Because neoepitopes are junction-derived
              they are novel relative to the <em>canonical</em> protein, but a peptide may still
              occur elsewhere in the normal proteome — check the proteome-uniqueness chip (above)
              before reading a candidate as tumour-specific. About 39 % of isoforms with predictions
              are GENCODE reference isoforms, not long-read novel ones.
            </p>
          </div>
        </section>

        {/* Epitope Atlas */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Epitope Atlas</h2>
          <p className="text-sm text-slate-400 leading-relaxed">
            The <strong className="text-foreground/80">Epitope Atlas</strong> (accessible via the
            Results page toggle) provides a proteome-wide view of all predicted neoantigens in a
            single searchable, filterable table — complementing the per-isoform neoantigen panel.
          </p>
          <div className="grid gap-4 sm:grid-cols-2">
            <div className="p-4 rounded-xl border border-border/40 bg-muted/10 space-y-1.5">
              <h4 className="font-bold text-sm text-foreground">Deduplication</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Each peptide × isoform pair is reduced to a single representative row — its
                best-scoring allele hit — so the same junction peptide can't inflate counts across
                every allele in the panel.
              </p>
            </div>
            <div className="p-4 rounded-xl border border-border/40 bg-muted/10 space-y-1.5">
              <h4 className="font-bold text-sm text-foreground">Sortable Columns</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Sort by priority tier, presentation (multi-tool consensus), immunogenicity, IC50,
                pLDDT, junction type, gene, or isoform. Default order is priority ascending (P1 first).
                Column headers show sort direction and carry inline tooltips explaining each metric.
              </p>
            </div>
            <div className="p-4 rounded-xl border border-border/40 bg-muted/10 space-y-1.5">
              <h4 className="font-bold text-sm text-foreground">Filters</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Filter simultaneously by gene symbol, priority tier (P1–P4), decorating flag
                (immunogenic, proteome-unique, all-three-tools, splice-enabled, expressed,
                PTM-disrupting, high-fidelity fold), binding category, HLA allele (dynamically
                loaded), junction type, and minimum pLDDT. Each row carries its flags as chips
                (IMM · SPL · EXP · PTM · UNQ / SELF / REF).
              </p>
            </div>
            <div className="p-4 rounded-xl border border-border/40 bg-muted/10 space-y-1.5">
              <h4 className="font-bold text-sm text-foreground">Summary Statistics</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Header chips show real-time counts per priority tier (P1–P3), immunogenic and
                expressed peptides, plus unique peptides and unique genes matching the active
                filter set.
              </p>
            </div>
          </div>
        </section>

        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "MHCflurry 2.2", note: "MHC-I binding & presentation", href: "https://github.com/openvax/mhcflurry" },
              { label: "NetMHCpan-4.2", note: "eluted-ligand presentation prediction", href: "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/" },
              { label: "BigMHC", note: "presentation (EL) & immunogenicity (IM)", href: "https://github.com/KarchinLab/bigmhc" },
              { label: "IEDB", note: "HLA reference allele set", href: "https://www.iedb.org/" },
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
          <Link href="/documentation/nmd-transcript-fate" className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2">
            <span className="opacity-60">&larr;</span> Prev: NMD &amp; Transcript Fate
          </Link>
          <div className="flex flex-col items-end gap-1">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest text-right">Next Section</span>
            <Link href="/documentation/ptm-accessibility" className="text-lg font-bold text-primary hover:underline flex items-center gap-2">
              PTM Accessibility <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
