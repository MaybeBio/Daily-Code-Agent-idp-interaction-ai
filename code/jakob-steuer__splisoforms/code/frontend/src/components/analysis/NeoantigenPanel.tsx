"use client";
import { useEffect, useState } from "react";
import {
  Syringe, ArrowUpDown, HelpCircle, CheckCircle2,
  AlertTriangle, Microscope, FlaskConical, Dna, Activity,
} from "lucide-react";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { cn } from "@/lib/utils";
import { getApiUrl } from "@/lib/api";

// ─── Types ────────────────────────────────────────────────────────────────────

interface NeoantigenPred {
  id: number;
  isoform_id: number;
  peptide_sequence: string;
  peptide_start: number;
  peptide_end: number;
  hla_allele: string;
  binding_affinity_nm: number;
  binding_category: string;
  percentile_rank?: number | null;
  junction_type: string;
  // MHCflurry 2.2
  presentation_score?: number | null;
  processing_score?:   number | null;
  // PTM
  ptm_impact_status?: string | null;
  avg_plddt?:         number | null;
  is_ptm_sensitive?:  boolean;
  // Multi-tool ensemble (v5)
  netmhcpan_affinity_nm?: number | null;
  netmhcpan_rank_el?:     number | null;
  netmhcpan_rank_ba?:     number | null;
  bigmhc_el_score?:       number | null;
  bigmhc_im_score?:       number | null;
  consensus_rank?:        number | null;
  is_ba_rescue?:          boolean;
  // Audit-3 additions
  tools_agreeing?:        number | null;
  agreement_tier?:        string | null;  // 'high' | 'medium' | 'low'
  ptm_canonical_residue_lost?: boolean;
  // ranking v2
  presentation_tier?: string | null;      // STRONG | MODERATE | WEAK | NONE
  priority_tier?:     string | null;       // P1 | P2 | P3 | P4
  is_immunogenic?:    boolean;
  is_expressed?:      boolean;
  min_plddt?:         number | null;
  // proteome-wide novelty vs every GENCODE reference protein (NULL = not computed)
  proteome_uniqueness?: "unique" | "self_gene" | "other_gene" | null;
  is_proteome_unique?:  boolean | null;
}

type SortKey =
  | "consensus_rank"
  | "presentation_score"
  | "bigmhc_im_score"
  | "processing_score"
  | "binding_affinity_nm"
  | "peptide_start"
  | "hla_allele";

// ranking v2 — priority tier badge styles
const PRIORITY_BADGE: Record<string, string> = {
  P1: "bg-emerald-500/15 text-emerald-300 border-emerald-500/30",
  P2: "bg-sky-500/15 text-sky-300 border-sky-500/30",
  P3: "bg-amber-500/15 text-amber-300 border-amber-500/30",
  P4: "bg-muted/40 text-muted-foreground border-border/50",
};
const PRIORITY_TITLE: Record<string, string> = {
  P1: "P1 · High-confidence — strong presentation + immunogenic",
  P2: "P2 · Strong presentation",
  P3: "P3 · Supported — moderate presentation + immunogenic",
  P4: "P4 · Candidate",
};

// ─── Shared helpers ───────────────────────────────────────────────────────────

function InfoTip({ text }: { text: string }) {
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger>
          <HelpCircle className="w-2.5 h-2.5 text-muted-foreground/50 hover:text-muted-foreground transition-colors cursor-help" />
        </TooltipTrigger>
        <TooltipContent className="bg-popover border-border text-foreground text-xs max-w-xs shadow-xl">
          {text}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}

function ScoreBar({ value, colorClass }: { value: number; colorClass: string }) {
  return (
    <div className="flex items-center gap-1.5 min-w-[72px]">
      <div className="flex-1 h-1.5 rounded-full bg-muted/50 overflow-hidden">
        <div className={`h-full rounded-full ${colorClass}`} style={{ width: `${Math.round(value * 100)}%` }} />
      </div>
      <span className="text-[10px] font-mono tabular-nums w-7 text-right">{(value * 100).toFixed(0)}%</span>
    </div>
  );
}

/** Colour-coded %rank chip — lower rank = better */
function RankChip({ rank, className }: { rank: number; className?: string }) {
  const { bg, text } =
    rank <= 0.5  ? { bg: "bg-emerald-500/20", text: "text-emerald-400" } :
    rank <= 2.0  ? { bg: "bg-amber-500/15",   text: "text-amber-400"   } :
    rank <= 10   ? { bg: "bg-orange-500/10",  text: "text-orange-400"  } :
                   { bg: "bg-muted/30",        text: "text-muted-foreground" };
  return (
    <span className={cn(`inline-flex items-center font-mono text-[10px] px-1.5 py-0.5 rounded border border-white/5 ${bg} ${text}`, className)}>
      {rank.toFixed(2)}%
    </span>
  );
}

/** Three-dot agreement indicator for a (peptide, allele) row.
 *
 *  Uses the persisted tools_agreeing / agreement_tier when available
 *  (audit-3, commit f5c2b8d6a3e4) so UI and SQL agree on what "agreement"
 *  means.  Thresholds:
 *      MHCflurry  %rank < 0.5
 *      NetMHCpan  %rank_el < 0.5
 *      BigMHC     EL score > 0.5
 */
function AgreementDots({ pred }: { pred: NeoantigenPred }) {
  const hasNetMHC = pred.netmhcpan_rank_el != null;
  const hasBigMHC = pred.bigmhc_el_score   != null;
  const hasFlurry = pred.percentile_rank   != null;

  const flurryOn = hasFlurry && (pred.percentile_rank!   < 0.5);
  const netmhcOn = hasNetMHC && (pred.netmhcpan_rank_el! < 0.5);
  const bigmhcOn = hasBigMHC && (pred.bigmhc_el_score!   > 0.5);

  const total = [hasFlurry, hasNetMHC, hasBigMHC].filter(Boolean).length;
  const agree = pred.tools_agreeing ?? [flurryOn, netmhcOn, bigmhcOn].filter(Boolean).length;
  const tier  = pred.agreement_tier ?? (agree >= 3 ? "high" : agree === 2 ? "medium" : "low");

  if (total === 0) return <span className="text-muted-foreground/30 text-[10px]">—</span>;

  const tierBadge =
    tier === "high"   ? "bg-emerald-500/15 text-emerald-300 border-emerald-500/30" :
    tier === "medium" ? "bg-amber-500/15   text-amber-300   border-amber-500/30"   :
                        "bg-muted/30        text-muted-foreground border-border/40";

  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger>
          <span className="inline-flex items-center gap-1.5 cursor-default">
            <span className="inline-flex items-center gap-0.5">
              {hasFlurry && <span className={cn("w-2 h-2 rounded-full", flurryOn ? "bg-sky-400"     : "bg-muted/40")} />}
              {hasNetMHC && <span className={cn("w-2 h-2 rounded-full", netmhcOn ? "bg-violet-400"  : "bg-muted/40")} />}
              {hasBigMHC && <span className={cn("w-2 h-2 rounded-full", bigmhcOn ? "bg-emerald-400" : "bg-muted/40")} />}
            </span>
            <span className={cn("text-[9px] font-bold uppercase tracking-wider px-1.5 py-0.5 rounded border", tierBadge)}>
              {tier}
            </span>
          </span>
        </TooltipTrigger>
        {/* flex-col overrides the base TooltipContent's `inline-flex items-center`, which
            otherwise lays these rows out horizontally and pushes everything past ~320px
            (BigMHC onward) outside the box. */}
        <TooltipContent className="flex flex-col items-stretch bg-popover border-border text-foreground text-xs shadow-xl space-y-1 p-2 min-w-[180px]">
          <div className="font-bold text-[10px] uppercase tracking-widest text-muted-foreground mb-1">Tool agreement</div>
          {hasFlurry && <div className="flex justify-between gap-4"><span className="text-sky-400">MHCflurry %rank</span><span>{flurryOn ? `✓ ${pred.percentile_rank!.toFixed(2)}%` : `– ${pred.percentile_rank!.toFixed(2)}%`}</span></div>}
          {hasNetMHC && <div className="flex justify-between gap-4"><span className="text-violet-400">NetMHCpan EL</span><span>{netmhcOn ? `✓ ${pred.netmhcpan_rank_el!.toFixed(2)}%` : `– ${pred.netmhcpan_rank_el!.toFixed(2)}%`}</span></div>}
          {hasBigMHC && <div className="flex justify-between gap-4"><span className="text-emerald-400">BigMHC EL</span><span>{bigmhcOn ? `✓ ${(pred.bigmhc_el_score! * 100).toFixed(0)}%` : `– ${(pred.bigmhc_el_score! * 100).toFixed(0)}%`}</span></div>}
          <div className="border-t border-border/30 pt-1 mt-1 text-muted-foreground">{agree}/{total} tools agree → {tier} confidence</div>
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}

/** Tool availability status dots in the header */
function ToolBadge({ label, available, color }: { label: string; available: boolean; color: string }) {
  return (
    <span className={cn(
      "inline-flex items-center gap-1.5 text-[10px] font-semibold px-2 py-1 rounded-md border",
      available
        ? `border-${color}-500/30 bg-${color}-500/10 text-${color}-400`
        : "border-border/30 bg-muted/20 text-muted-foreground/40"
    )}>
      <span className={cn("w-1.5 h-1.5 rounded-full", available ? `bg-${color}-400` : "bg-muted/50")} />
      {label}
    </span>
  );
}

// ─── Sort header ──────────────────────────────────────────────────────────────

function SortHeader({
  k, label, sortKey, sortAsc, onSort, className = "", tip, alignRight = false,
}: {
  k: SortKey; label: string; sortKey: SortKey; sortAsc: boolean;
  onSort: (k: SortKey) => void; className?: string; tip?: string; alignRight?: boolean;
}) {
  const active = sortKey === k;
  return (
    <th className={`py-2 px-2 cursor-pointer select-none ${className}`} onClick={() => onSort(k)}>
      <span className={cn("flex items-center gap-1 text-[10px] uppercase tracking-wider font-bold", alignRight && "justify-end")}>
        {label}
        {tip && <InfoTip text={tip} />}
        <ArrowUpDown className={cn("w-2.5 h-2.5", active ? "text-primary opacity-100" : "opacity-40")} />
        {active && <span className="text-primary">{sortAsc ? "↑" : "↓"}</span>}
      </span>
    </th>
  );
}

// ─── Main component ───────────────────────────────────────────────────────────

export default function NeoantigenPanel({
  isoformId,
  isCanonical = false,
  hoveredResidues,
  selectedResidues,
  onHover,
  onSelect
}: {
  isoformId: number;
  isCanonical?: boolean;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelect?: (residues: number[]) => void;
}) {
  const [predictions, setPredictions] = useState<NeoantigenPred[]>([]);
  const [loading,     setLoading]     = useState(true);
  // Default sort: presentation_score (MHCflurry's calibrated metric) instead of
  // consensus_rank — the latter combines incommensurable scales (audit-3).
  const [sortKey,     setSortKey]     = useState<SortKey>("presentation_score");
  // presentation_score is higher-is-better → descending shows the strongest first
  // (matches toggleSort, which sets ascending only for %rank/affinity/position keys).
  const [sortAsc,     setSortAsc]     = useState(false);

  useEffect(() => {
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/isoform/${isoformId}/neoantigen`)
      .then(r => r.ok ? r.json() : [])
      .then((data: NeoantigenPred[]) => {
        // Presentation_score is the default; keep it that way regardless of
        // whether consensus_rank is populated (consensus is informational only).
        setPredictions(data);
      })
      .catch(() => setPredictions([]))
      .finally(() => setLoading(false));
  }, [isoformId]);

  if (loading) return (
    <div className="text-sm text-muted-foreground animate-pulse p-4">Loading neoantigen predictions…</div>
  );
  if (predictions.length === 0) return (
    <div className="p-6 text-center text-muted-foreground">
      <Syringe className="w-8 h-8 mx-auto mb-2 opacity-40" />
      {isCanonical ? (
        <>
          <p className="text-sm">No neoepitopes — this is the canonical isoform.</p>
          <p className="text-xs mt-1 opacity-60 max-w-md mx-auto">
            Neoepitopes are peptides absent from the canonical protein. The canonical isoform is
            the reference it would be compared against, so by definition it has none.
          </p>
        </>
      ) : (
        <>
          <p className="text-sm">No neoepitopes predicted for this isoform.</p>
          <p className="text-xs mt-1 opacity-60 max-w-md mx-auto">
            Either its sequence carries no region absent from the canonical protein, or the novel
            regions it does have yielded no predicted MHC binders.
          </p>
        </>
      )}
    </div>
  );

  // ── Derived stats ─────────────────────────────────────────────────────────
  const strong        = predictions.filter(p => p.binding_category === "strong_binder").length;
  const weak          = predictions.filter(p => p.binding_category === "weak_binder").length;
  const highConsensus = predictions.filter(p => p.consensus_rank != null && p.consensus_rank <= 2.0).length;
  const presented     = predictions.filter(p => p.presentation_score != null && p.presentation_score >= 0.5).length;
  const immunogenic   = predictions.filter(p => p.bigmhc_im_score != null && p.bigmhc_im_score >= 0.5).length;
  const baRescued     = predictions.filter(p => p.is_ba_rescue).length;

  const hasFlurry    = predictions.some(p => p.presentation_score != null);
  const hasNetMHC    = predictions.some(p => p.netmhcpan_rank_el != null);
  const hasBigMHC    = predictions.some(p => p.bigmhc_el_score != null);
  const hasConsensus = predictions.some(p => p.consensus_rank != null);
  const hasIM        = predictions.some(p => p.bigmhc_im_score != null);

  const toolCount = [hasFlurry, hasNetMHC, hasBigMHC].filter(Boolean).length;

  // ── Sort ──────────────────────────────────────────────────────────────────
  const sorted = [...predictions].sort((a, b) => {
    const nullLast = (v: number | null | undefined) => v ?? (sortAsc ? Infinity : -Infinity);
    let va: number, vb: number;
    switch (sortKey) {
      case "consensus_rank":      va = nullLast(a.consensus_rank);      vb = nullLast(b.consensus_rank);      break;
      case "presentation_score":  va = nullLast(a.presentation_score);  vb = nullLast(b.presentation_score);  break;
      case "bigmhc_im_score":     va = nullLast(a.bigmhc_im_score);     vb = nullLast(b.bigmhc_im_score);     break;
      case "processing_score":    va = nullLast(a.processing_score);    vb = nullLast(b.processing_score);    break;
      case "binding_affinity_nm": va = a.binding_affinity_nm;           vb = b.binding_affinity_nm;           break;
      case "peptide_start":       va = a.peptide_start;                 vb = b.peptide_start;                 break;
      case "hla_allele":
        return sortAsc
          ? a.hla_allele.localeCompare(b.hla_allele)
          : b.hla_allele.localeCompare(a.hla_allele);
      default: va = 0; vb = 0;
    }
    return sortAsc ? va - vb : vb - va;
  });

  const toggleSort = (k: SortKey) => {
    if (sortKey === k) { setSortAsc(!sortAsc); return; }
    setSortKey(k);
    // %rank and affinity: lower is better → ascending. Scores: higher is better → descending.
    setSortAsc(k === "consensus_rank" || k === "binding_affinity_nm" || k === "peptide_start");
  };

  // ── Render ────────────────────────────────────────────────────────────────
  return (
    <div className="space-y-5 p-1">

      {/* ── Tool status bar ─────────────────────────────────────────────── */}
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div className="flex flex-wrap gap-1.5">
          <span className="text-[10px] uppercase font-bold text-muted-foreground/60 tracking-widest self-center mr-1">
            Prediction tools
          </span>
          {/* Use static class names to avoid Tailwind purge issues */}
          <span className={cn(
            "inline-flex items-center gap-1.5 text-[10px] font-semibold px-2 py-1 rounded-md border",
            hasFlurry
              ? "border-sky-500/30 bg-sky-500/10 text-sky-400"
              : "border-border/30 bg-muted/20 text-muted-foreground/40"
          )}>
            <span className={cn("w-1.5 h-1.5 rounded-full", hasFlurry ? "bg-sky-400" : "bg-muted/50")} />
            MHCflurry 2.2
          </span>
          <span className={cn(
            "inline-flex items-center gap-1.5 text-[10px] font-semibold px-2 py-1 rounded-md border",
            hasNetMHC
              ? "border-violet-500/30 bg-violet-500/10 text-violet-400"
              : "border-border/30 bg-muted/20 text-muted-foreground/40"
          )}>
            <span className={cn("w-1.5 h-1.5 rounded-full", hasNetMHC ? "bg-violet-400" : "bg-muted/50")} />
            NetMHCpan-4.2
          </span>
          <span className={cn(
            "inline-flex items-center gap-1.5 text-[10px] font-semibold px-2 py-1 rounded-md border",
            hasBigMHC
              ? "border-emerald-500/30 bg-emerald-500/10 text-emerald-400"
              : "border-border/30 bg-muted/20 text-muted-foreground/40"
          )}>
            <span className={cn("w-1.5 h-1.5 rounded-full", hasBigMHC ? "bg-emerald-400" : "bg-muted/50")} />
            BigMHC EL+IM
          </span>
        </div>
        {toolCount < 3 && (
          <span className="text-[10px] text-muted-foreground/50 italic">
            {toolCount}/3 predictors scored this isoform
          </span>
        )}
      </div>

      {/* ── Summary chips ───────────────────────────────────────────────── */}
      <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-6 gap-2">
        {[
          { label: "Strong binders",    value: strong,        sub: "< 50 nM IC50",       color: "text-rose-400",    bg: "bg-rose-500/10",    border: "border-rose-500/20"    },
          { label: "Weak binders",      value: weak,          sub: "50–500 nM",           color: "text-amber-400",   bg: "bg-amber-500/10",   border: "border-amber-500/20"   },
          ...(hasConsensus ? [{ label: "High consensus", value: highConsensus, sub: "rank ≤2%, all tools",  color: "text-emerald-400", bg: "bg-emerald-500/10", border: "border-emerald-500/20" }] : []),
          ...(hasFlurry    ? [{ label: "Presented",      value: presented,    sub: "MHCflurry ≥50%",       color: "text-sky-400",     bg: "bg-sky-500/10",     border: "border-sky-500/20"     }] : []),
          ...(hasIM        ? [{ label: "Immunogenic",    value: immunogenic,  sub: "BigMHC-IM ≥50%",       color: "text-purple-400",  bg: "bg-purple-500/10",  border: "border-purple-500/20"  }] : []),
          ...(baRescued > 0 ? [{ label: "BA-rescued",   value: baRescued,    sub: "Cys-rich, EL-missed",  color: "text-orange-400",  bg: "bg-orange-500/10",  border: "border-orange-500/20"  }] : []),
        ].map((chip, i) => (
          <div key={i} className={`rounded-lg px-3 py-2 border ${chip.bg} ${chip.border} flex flex-col gap-0.5`}>
            <span className={`text-xl font-bold ${chip.color}`}>{chip.value}</span>
            <span className="text-[10px] font-bold text-foreground/70 leading-tight">{chip.label}</span>
            <span className="text-[9px] text-muted-foreground/60">{chip.sub}</span>
          </div>
        ))}
      </div>

      {/* ── BA-rescue notice ────────────────────────────────────────────── */}
      {baRescued > 0 && (
        <div className="flex items-start gap-2 text-[11px] bg-orange-500/5 border border-orange-500/15 rounded-lg px-3 py-2">
          <AlertTriangle className="w-3.5 h-3.5 mt-0.5 text-orange-400 shrink-0" />
          <span className="text-muted-foreground">
            <span className="font-semibold text-orange-400">{baRescued} BA-rescue candidate{baRescued > 1 ? "s" : ""}</span>
            {" "}— cysteine-rich peptide{baRescued > 1 ? "s" : ""} predicted by NetMHCpan binding-affinity mode
            but missed by eluted-ligand tools, possibly due to MS cysteine-depletion bias.
            Marked with <span className="font-mono bg-orange-500/20 text-orange-300 px-1 rounded">BA</span> in the table below.
          </span>
        </div>
      )}

      {/* ── Main table ──────────────────────────────────────────────────── */}
      <div className="overflow-x-auto max-h-[520px] overflow-y-auto rounded-lg border border-border/40">
        <table className="w-full text-xs">
          <thead className="sticky top-0 bg-card z-10 border-b border-border">
            <tr className="text-muted-foreground">
              <th className="text-left py-2 pl-3 pr-2 text-[10px] uppercase tracking-wider font-bold">Peptide</th>
              <SortHeader k="peptide_start" label="Pos" sortKey={sortKey} sortAsc={sortAsc} onSort={toggleSort} />
              <SortHeader k="hla_allele"    label="HLA" sortKey={sortKey} sortAsc={sortAsc} onSort={toggleSort} />

              {hasConsensus && (
                <SortHeader
                  k="consensus_rank" label="Presentation" sortKey={sortKey} sortAsc={sortAsc} onSort={toggleSort}
                  tip="Multi-tool presentation consensus: STRONG = ≥2 tools call a strong binder · MOD = 1 strong or ≥2 weak · WEAK = 1 weak. Hover for the underlying median %rank; tool agreement is the next column."
                />
              )}

              {/* Agreement dots — always show when ≥2 tools */}
              {toolCount >= 2 && (
                <th className="py-2 px-2 text-[10px] uppercase tracking-wider font-bold text-muted-foreground">
                  <span className="flex items-center gap-1">
                    Agreement
                    <InfoTip text="Per-tool signal: sky = MHCflurry, violet = NetMHCpan, green = BigMHC. Filled = above threshold. Hover for values." />
                  </span>
                </th>
              )}

              {hasFlurry && (
                <SortHeader
                  k="presentation_score" label="MHCflurry" sortKey={sortKey} sortAsc={sortAsc} onSort={toggleSort}
                  tip="MHCflurry 2.2 presentation (0–1): binding + processing. ≥0.5 = likely presented."
                />
              )}

              {hasNetMHC && (
                <th className="py-2 px-2 text-[10px] uppercase tracking-wider font-bold text-muted-foreground">
                  <span className="flex items-center gap-1">
                    NetMHCpan
                    <InfoTip text="NetMHCpan-4.2 EL %rank (lower = better). BA %rank shown as a rescue signal for cysteine-rich peptides." />
                  </span>
                </th>
              )}

              {hasIM && (
                <SortHeader
                  k="bigmhc_im_score" label="Immunog." sortKey={sortKey} sortAsc={sortAsc} onSort={toggleSort}
                  tip="BigMHC immunogenicity (0–1): probability a presented peptide activates T-cells. Orthogonal to available in MHCflurry or NetMHCpan."
                />
              )}

              <SortHeader
                k="binding_affinity_nm" label="IC50 (nM)" sortKey={sortKey} sortAsc={sortAsc} onSort={toggleSort}
                alignRight
                tip="MHCflurry 2.2 IC50 (nM). <50 = strong, <500 = weak."
              />
              <th className="text-left py-2 px-2 text-[10px] uppercase tracking-wider font-bold text-muted-foreground">Junction</th>
            </tr>
          </thead>

          <tbody>
            {sorted.map(p => {
              const isHighConf  = p.consensus_rank != null ? p.consensus_rank <= 2.0
                                : p.presentation_score    != null && p.presentation_score >= 0.5;
              const isImm       = p.is_immunogenic ?? (p.bigmhc_im_score != null && p.bigmhc_im_score >= 0.5);
              const spliceEnabled = (p.ptm_impact_status ?? "").startsWith("Splice-Enabled");
              const isCysRich   = p.peptide_sequence.split("C").length - 1 >= 2;

              const range = Array.from({length: p.peptide_end - p.peptide_start + 1}, (_, i) => p.peptide_start + i);
              const isHovered = hoveredResidues?.some(r => range.includes(r));
              const isSelected = selectedResidues?.some(r => range.includes(r));
              
              // We use bg colors or outlines to highlight the row
              const highlightClass = isSelected ? "bg-primary/15 ring-2 ring-primary ring-inset" : isHovered ? "bg-primary/5 ring-1 ring-primary/50 ring-inset" : "";

              return (
                <tr
                  key={p.id}
                  className={cn(
                    "border-b border-border/30 transition-colors cursor-pointer relative",
                    isSelected ? "bg-muted/60" : p.is_ba_rescue ? "bg-orange-500/5 hover:bg-orange-500/10" : "hover:bg-muted/20",
                    highlightClass
                  )}
                  onMouseEnter={() => onHover && onHover(range)}
                  onMouseLeave={() => onHover && onHover([])}
                  onClick={() => {
                    if (onSelect) {
                      onSelect(isSelected ? [] : range);
                    }
                  }}
                >
                  {/* Peptide sequence */}
                  <td className="py-2 pl-3 pr-2">
                    <div className="flex items-center gap-1.5 flex-wrap">
                      {p.priority_tier && (
                        <span title={PRIORITY_TITLE[p.priority_tier]}
                          className={cn("text-[9px] font-bold px-1 rounded border", PRIORITY_BADGE[p.priority_tier] ?? PRIORITY_BADGE.P4)}>
                          {p.priority_tier}
                        </span>
                      )}
                      <span className="font-mono font-medium text-foreground tracking-wide">
                        {p.peptide_sequence.split("").map((aa, i) => (
                          <span key={i} className={aa === "C" && p.is_ba_rescue ? "text-orange-400 font-bold" : ""}>{aa}</span>
                        ))}
                      </span>
                      {p.is_ba_rescue && (
                        <span className="text-[9px] font-bold bg-orange-500/20 text-orange-300 border border-orange-500/30 px-1 rounded">
                          BA
                        </span>
                      )}
                      {isImm && (
                        <span title="Predicted immunogenic (BigMHC-IM ≥ 0.5)" className="text-[9px] font-bold bg-violet-500/20 text-violet-300 border border-violet-500/30 px-1 rounded">
                          IMM
                        </span>
                      )}
                      {spliceEnabled && (
                        <span title="Splice-enabled: epitope exists only because of the splice event" className="text-[9px] font-bold bg-orange-500/10 text-orange-300 border border-orange-500/25 px-1 rounded">
                          SPL
                        </span>
                      )}
                      {p.is_expressed && (
                        <span title="Isoform detected in expression data" className="text-[9px] font-bold bg-emerald-500/10 text-emerald-300 border border-emerald-500/25 px-1 rounded">
                          EXP
                        </span>
                      )}
                      {p.proteome_uniqueness === "unique" && (
                        <span title="Proteome-unique: the peptide occurs in no GENCODE v47 reference protein" className="text-[9px] font-bold bg-rose-500/10 text-rose-300 border border-rose-500/25 px-1 rounded">
                          UNQ
                        </span>
                      )}
                      {(p.proteome_uniqueness === "self_gene" || p.proteome_uniqueness === "other_gene") && (
                        <span title={p.proteome_uniqueness === "self_gene"
                          ? "Also present in another annotated isoform of the same gene (reference self-peptide)"
                          : "Present in a reference protein of another gene (self-peptide; not tumour-specific)"}
                          className="text-[9px] font-bold bg-muted/30 text-muted-foreground border border-border/40 px-1 rounded">
                          {p.proteome_uniqueness === "self_gene" ? "SELF" : "REF"}
                        </span>
                      )}
                    </div>
                  </td>

                  {/* Position */}
                  <td className="py-2 px-2 font-mono text-foreground/60 whitespace-nowrap">
                    {p.peptide_start}–{p.peptide_end}
                  </td>

                  {/* HLA */}
                  <td className="py-2 px-2 font-mono text-foreground/80 whitespace-nowrap text-[10px]">
                    {p.hla_allele}
                  </td>

                  {/* Presentation tier (multi-tool consensus; median %rank on hover) */}
                  {hasConsensus && (
                    <td className="py-2 px-2">
                      {p.presentation_tier && p.presentation_tier !== "NONE" ? (
                        <span
                          title={p.consensus_rank != null ? `median %rank ${p.consensus_rank.toFixed(2)}%` : undefined}
                          className={cn(
                            "inline-flex items-center rounded px-1.5 py-0.5 text-[9px] font-bold uppercase tracking-wide border",
                            p.presentation_tier === "STRONG"   ? "bg-emerald-500/15 text-emerald-300 border-emerald-500/30" :
                            p.presentation_tier === "MODERATE" ? "bg-amber-500/15 text-amber-300 border-amber-500/30" :
                                                                 "bg-muted/40 text-muted-foreground border-border/50"
                          )}>
                          {p.presentation_tier === "STRONG" ? "Strong" : p.presentation_tier === "MODERATE" ? "Mod" : "Weak"}
                        </span>
                      ) : <span className="text-muted-foreground/30 text-[10px]">—</span>}
                    </td>
                  )}

                  {/* Agreement dots */}
                  {toolCount >= 2 && (
                    <td className="py-2 px-2">
                      <AgreementDots pred={p} />
                    </td>
                  )}

                  {/* MHCflurry presentation score */}
                  {hasFlurry && (
                    <td className="py-2 px-2">
                      {p.presentation_score != null ? (
                        <ScoreBar
                          value={p.presentation_score}
                          colorClass={
                            p.presentation_score >= 0.7 ? "bg-sky-500"
                            : p.presentation_score >= 0.5 ? "bg-sky-400/70"
                            : "bg-muted-foreground/40"
                          }
                        />
                      ) : <span className="text-muted-foreground/30 text-[10px]">n/a</span>}
                    </td>
                  )}

                  {/* NetMHCpan EL + BA ranks */}
                  {hasNetMHC && (
                    <td className="py-2 px-2 space-y-0.5">
                      {p.netmhcpan_rank_el != null ? (
                        <div className="flex items-center gap-1">
                          <span className="text-[9px] text-muted-foreground/50 w-4">EL</span>
                          <RankChip rank={p.netmhcpan_rank_el} />
                        </div>
                      ) : null}
                      {p.netmhcpan_rank_ba != null ? (
                        <div className="flex items-center gap-1">
                          <span className="text-[9px] text-muted-foreground/50 w-4">BA</span>
                          <RankChip rank={p.netmhcpan_rank_ba} />
                        </div>
                      ) : null}
                      {p.netmhcpan_rank_el == null && p.netmhcpan_rank_ba == null && (
                        <span className="text-muted-foreground/30 text-[10px]">—</span>
                      )}
                    </td>
                  )}

                  {/* BigMHC IM immunogenicity */}
                  {hasIM && (
                    <td className="py-2 px-2">
                      {p.bigmhc_im_score != null ? (
                        <ScoreBar
                          value={p.bigmhc_im_score}
                          colorClass={
                            p.bigmhc_im_score >= 0.7 ? "bg-purple-500"
                            : p.bigmhc_im_score >= 0.5 ? "bg-purple-400/70"
                            : "bg-muted-foreground/40"
                          }
                        />
                      ) : <span className="text-muted-foreground/30 text-[10px]">—</span>}
                    </td>
                  )}

                  {/* IC50 */}
                  <td className="py-2 px-2 text-right font-mono font-medium">
                    <span className={
                      p.binding_affinity_nm <= 50 ? "text-rose-400"
                      : p.binding_affinity_nm <= 500 ? "text-amber-400"
                      : "text-muted-foreground/60"
                    }>
                      {p.binding_affinity_nm.toFixed(1)}
                    </span>
                  </td>

                  {/* Junction */}
                  <td className="py-2 px-2 text-[10px] text-muted-foreground whitespace-nowrap">
                    {p.junction_type.replace(/_/g, " ")}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* ── Legend ──────────────────────────────────────────────────────── */}
      <div className="grid grid-cols-2 md:grid-cols-3 gap-x-6 gap-y-1.5 text-[10px] text-muted-foreground pt-1 border-t border-border/30">
        {hasConsensus && <>
          <span className="flex items-center gap-1.5"><span className="w-2 h-2 rounded bg-emerald-500/80" />Presentation STRONG — ≥2 tools call a strong binder</span>
          <span className="flex items-center gap-1.5"><span className="w-2 h-2 rounded bg-amber-500/80" />Presentation MOD — 1 strong or ≥2 weak</span>
          <span className="flex items-center gap-1.5"><span className="w-2 h-2 rounded bg-muted-foreground/50" />Presentation WEAK — 1 weak tool call</span>
        </>}
        {hasFlurry && <span className="flex items-center gap-1.5"><span className="w-2 h-2 rounded-full bg-sky-400" />MHCflurry: presentation score (blue bar)</span>}
        {hasIM     && <span className="flex items-center gap-1.5"><span className="w-2 h-2 rounded-full bg-purple-500" />BigMHC IM: T-cell immunogenicity (purple bar)</span>}
        {baRescued > 0 && <span className="flex items-center gap-1.5"><span className="font-mono text-[9px] bg-orange-500/20 text-orange-300 px-1 rounded">BA</span>Cys-rich peptide rescued by BA-mode; orange C in sequence</span>}
        {hasIM     && <span className="flex items-center gap-1.5"><span className="font-mono text-[9px] bg-purple-500/20 text-purple-300 px-1 rounded">IMM</span>BigMHC-IM ≥50%: predicted immunogenic</span>}
      </div>
    </div>
  );
}
