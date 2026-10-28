"use client";

import { useEffect, useState, useCallback } from "react";
import { useSearchParams } from "next/navigation";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import {
  Syringe, Filter, Search,
  Sparkles, Zap, Shield, FlaskConical, Activity, ShieldCheck,
  ChevronLeft, ChevronRight, Info,
  ArrowUp, ArrowDown, ArrowUpDown, HelpCircle,
} from "lucide-react";
import { cn } from "@/lib/utils";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { getApiUrl } from "@/lib/api";

// ─── Types ────────────────────────────────────────────────────────────────────

type SortCol = "priority" | "consensus" | "immunogenicity" | "affinity" | "plddt" | "gene" | "isoform" | "peptide" | "hla" | "junction";

interface NeoantigenRow {
  id: number;
  isoform_id: number;
  peptide_sequence: string;
  peptide_start: number;
  peptide_end: number;
  hla_allele: string;
  affinity_nm: number;
  binding_category: string;
  junction_type: string;
  ptm_impact_status?: string;
  presentation_score?: number | null;
  avg_plddt?: number | null;
  min_plddt?: number | null;
  is_ptm_sensitive?: boolean;
  ptm_canonical_residue_lost?: boolean;
  transcript_id: string;
  transcript_name?: string | null;
  gene_id: number;
  gene_symbol: string;
  // ranking v2
  presentation_tier?: string | null;
  priority_tier?: string | null;
  is_immunogenic?: boolean;
  is_expressed?: boolean;
  priority_rank: number;
  consensus_rank?: number | null;
  bigmhc_im_score?: number | null;
  netmhcpan_rank_el?: number | null;
  is_ba_rescue?: boolean;
  // proteome-wide novelty of the peptide vs every GENCODE reference protein
  proteome_uniqueness?: "unique" | "self_gene" | "other_gene" | null;
  is_proteome_unique?: boolean | null;
  // how many of the 3 predictors actually scored this row (NULL tools do not vote)
  tools_available?: number | null;
  tools_agreeing?: number | null;
}

interface BrowserStats {
  total: number;
  p1: number;
  p2: number;
  p3: number;
  p4: number;
  strong_presentation: number;
  immunogenic: number;
  expressed: number;
  splice_enabled: number;
  ba_rescue?: number;
  proteome_unique?: number;
  missing_netmhcpan?: number;
  unique_peptides: number;
  unique_genes: number;
}

interface ApiResponse {
  items: NeoantigenRow[];
  total_count: number;
  stats: BrowserStats;
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

// priority_tier P1–P4: the headline tier (presentation consensus × immunogenicity)
const PRIORITY_META: Record<string, { label: string; color: string; bg: string; border: string }> = {
  P1: { label: "High-confidence",     color: "text-emerald-400",      bg: "bg-emerald-500/15", border: "border-emerald-500/30" },
  P2: { label: "Strong presentation", color: "text-sky-400",          bg: "bg-sky-500/10",     border: "border-sky-500/25"     },
  P3: { label: "Supported",           color: "text-amber-400",        bg: "bg-amber-500/10",   border: "border-amber-500/25"   },
  P4: { label: "Candidate",           color: "text-muted-foreground", bg: "bg-muted/30",       border: "border-border/40"      },
};
function priorityMeta(t?: string | null) { return PRIORITY_META[t ?? ""] ?? PRIORITY_META.P4; }

const PRESENTATION_LABEL: Record<string, string> = {
  STRONG: "Strong", MODERATE: "Moderate", WEAK: "Weak", NONE: "None",
};

function AffinityBadge({ nm, category }: { nm: number; category: string }) {
  const cls = category === "strong_binder"
    ? "text-rose-400 bg-rose-500/15 border-rose-500/30"
    : "text-amber-400 bg-amber-500/10 border-amber-500/25";
  return (
    <span className={`inline-flex items-center font-mono text-[10px] px-1.5 py-0.5 rounded border ${cls}`}>
      {nm.toFixed(1)} nM
    </span>
  );
}

function PlddtDot({ value }: { value?: number | null }) {
  if (value == null) return <span className="text-muted-foreground/40 text-[10px]">—</span>;
  const dot =
    value >= 90 ? "bg-[#0053D6]"
    : value >= 70 ? "bg-[#65CBF3]"
    : value >= 50 ? "bg-[#FFDB13]"
    : "bg-[#FF7D45]";
  return (
    <span className="flex items-center gap-1">
      <span className={`w-1.5 h-1.5 rounded-full ${dot} inline-block`} />
      <span className="font-mono text-[10px]">{value.toFixed(0)}</span>
    </span>
  );
}

// decorating flags rendered as small chips — orthogonal to the tier
function FlagChips({ row }: { row: NeoantigenRow }) {
  const spliceEnabled = (row.ptm_impact_status ?? "").startsWith("Splice-Enabled");
  const chips: { label: string; title: string; cls: string }[] = [];
  if (row.is_immunogenic)
    chips.push({ label: "IMM", title: "Predicted immunogenic (BigMHC-IM ≥ 0.5)", cls: "text-violet-400 bg-violet-500/10 border-violet-500/25" });
  if (spliceEnabled)
    chips.push({ label: "SPL", title: "Splice-enabled: epitope exists only because of the splice event", cls: "text-orange-400 bg-orange-500/10 border-orange-500/25" });
  if (row.is_expressed)
    chips.push({ label: "EXP", title: "Isoform detected in expression data", cls: "text-emerald-400 bg-emerald-500/10 border-emerald-500/25" });
  if (row.ptm_canonical_residue_lost)
    chips.push({ label: "PTM", title: "Ablates a canonical PTM site", cls: "text-cyan-400 bg-cyan-500/10 border-cyan-500/25" });
  if (row.proteome_uniqueness === "unique")
    chips.push({ label: "UNQ", title: "Proteome-unique: the peptide occurs in no GENCODE v47 reference protein", cls: "text-rose-400 bg-rose-500/10 border-rose-500/25" });
  else if (row.proteome_uniqueness === "self_gene")
    chips.push({ label: "SELF", title: "Also present in another annotated isoform of the same gene (reference self-peptide)", cls: "text-muted-foreground bg-muted/30 border-border/40" });
  else if (row.proteome_uniqueness === "other_gene")
    chips.push({ label: "REF", title: "Present in a reference protein of another gene (self-peptide; not tumour-specific)", cls: "text-muted-foreground bg-muted/30 border-border/40" });
  if (row.tools_available != null && row.tools_available < 3)
    chips.push({ label: `${row.tools_available}/3`, title: `Only ${row.tools_available} of 3 predictors scored this row (NetMHCpan missing); the tier is a ${row.tools_available}-tool vote`, cls: "text-amber-400 bg-amber-500/10 border-amber-500/25" });
  if (!chips.length) return <span className="text-muted-foreground/30 text-[10px]">—</span>;
  return (
    <span className="flex flex-wrap gap-1">
      {chips.map(c => (
        <span key={c.label} title={c.title}
          className={cn("inline-flex items-center px-1 py-0.5 rounded text-[8px] font-bold border", c.cls)}>
          {c.label}
        </span>
      ))}
    </span>
  );
}

function StatChip({ label, value, sub, color }: {
  label: string; value: string | number; sub?: string; color: string;
}) {
  return (
    <div className="flex flex-col gap-0.5 bg-card/40 border border-border/40 rounded-lg px-4 py-3 min-w-[110px]">
      <span className={`text-xl font-bold tabular-nums ${color}`}>
        {typeof value === "number" ? value.toLocaleString() : value}
      </span>
      <span className="text-[9px] font-bold uppercase tracking-widest text-muted-foreground">{label}</span>
      {sub && <span className="text-[9px] text-muted-foreground/60">{sub}</span>}
    </div>
  );
}

// ─── Inline tooltip ───────────────────────────────────────────────────────────

function InfoTip({ text }: { text: string }) {
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger>
          <HelpCircle className="w-2.5 h-2.5 text-muted-foreground/40 hover:text-muted-foreground transition-colors cursor-help" />
        </TooltipTrigger>
        <TooltipContent className="bg-popover border-border text-foreground text-xs max-w-xs shadow-xl font-normal normal-case tracking-normal">
          {text}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}

// ─── Column header with sort indicator ───────────────────────────────────────

function ColHeader({
  col, label, sortCol, sortOrder, onSort, className = "", tip, alignRight = false, alignCenter = false,
}: {
  col: SortCol;
  label: string;
  sortCol: SortCol;
  sortOrder: "asc" | "desc";
  onSort: (col: SortCol) => void;
  className?: string;
  tip?: string;
  alignRight?: boolean;
  alignCenter?: boolean;
}) {
  const active = sortCol === col;
  return (
    <th
      onClick={() => onSort(col)}
      className={cn(
        "py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] cursor-pointer select-none whitespace-nowrap",
        "hover:text-foreground transition-colors",
        active ? "text-primary" : "text-muted-foreground",
        className,
      )}
    >
      <span className={cn("flex items-center gap-1", alignRight && "justify-end", alignCenter && "justify-center")}>
        {label}
        {tip && <InfoTip text={tip} />}
        {active
          ? sortOrder === "asc"
            ? <ArrowUp   className="w-2.5 h-2.5" />
            : <ArrowDown className="w-2.5 h-2.5" />
          : <ArrowUpDown className="w-2.5 h-2.5 opacity-30" />
        }
      </span>
    </th>
  );
}

// ─── Constants ────────────────────────────────────────────────────────────────

const PRIORITY_FILTERS = [
  { id: "",   label: "All tiers",           icon: Filter   },
  { id: "P1", label: "P1 · High-confidence", icon: Sparkles },
  { id: "P2", label: "P2 · Strong pres.",    icon: ShieldCheck },
  { id: "P3", label: "P3 · Supported",       icon: Activity },
  { id: "P4", label: "P4 · Candidate",       icon: Filter   },
];

const FLAG_FILTERS = [
  { id: "",              label: "Any",           icon: Filter       },
  { id: "immunogenic",   label: "Immunogenic",   icon: Activity     },
  { id: "splice_enabled", label: "Splice-enabled", icon: Zap        },
  { id: "expressed",     label: "Expressed",     icon: FlaskConical },
  { id: "ptm_disrupting", label: "PTM-disrupting", icon: Shield     },
  { id: "high_fidelity", label: "High-fidelity fold", icon: ShieldCheck },
  { id: "proteome_unique", label: "Proteome-unique", icon: ShieldCheck },
  { id: "all_tools",     label: "All 3 tools scored", icon: Activity  },
];

// HLA alleles are fetched dynamically from /api/neoantigens/alleles

const LIMIT = 50;

// ─── Default sort direction per column ────────────────────────────────────────
const DEFAULT_ORDER: Record<SortCol, "asc" | "desc"> = {
  priority:       "asc",   // P1 first
  consensus:      "asc",   // lower consensus rank = better
  immunogenicity: "desc",  // higher immunogenicity score = better
  affinity:       "asc",   // lower nM = better
  plddt:          "desc",
  gene:           "asc",
  isoform:        "asc",
  peptide:        "asc",
  hla:            "asc",
  junction:       "asc",
};

// ─── Main component ───────────────────────────────────────────────────────────

export default function GlobalNeoantigenBrowser() {
  const [data, setData]       = useState<ApiResponse | null>(null);
  const [loading, setLoading] = useState(true);
  const [hlaOptions, setHlaOptions] = useState<string[]>([]);

  useEffect(() => {
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/neoantigens/alleles`)
      .then(r => r.ok ? r.json() : [])
      .then((rows: { allele: string }[]) => setHlaOptions(rows.map(r => r.allele)))
      .catch(() => {});
  }, []);

  const initialGene = (useSearchParams().get("gene") || "").trim();
  const [geneSearch,   setGeneSearch]   = useState(initialGene);
  const [geneApplied,  setGeneApplied]  = useState(initialGene);
  const [bindingCat,   setBindingCat]   = useState<"" | "strong_binder" | "weak_binder">("");
  const [hlaAllele,    setHlaAllele]    = useState("");
  const [junctionType, setJunctionType] = useState("");
  const [priorityTier, setPriorityTier] = useState("");
  const [flag,         setFlag]         = useState("");
  const [minPlddt,     setMinPlddt]     = useState(0);

  const [sortCol,   setSortCol]   = useState<SortCol>("priority");
  const [sortOrder, setSortOrder] = useState<"asc" | "desc">("asc");

  const [page, setPage] = useState(0);

  const fetchData = useCallback(async () => {
    setLoading(true);
    try {
      const apiUrl = getApiUrl();
      const params = new URLSearchParams({
        limit:      String(LIMIT),
        offset:     String(page * LIMIT),
        sort_by:    sortCol,
        sort_order: sortOrder,
      });
      if (bindingCat)   params.set("binding_category", bindingCat);
      if (hlaAllele)    params.set("hla_allele", hlaAllele);
      if (junctionType) params.set("junction_type", junctionType);
      if (priorityTier) params.set("priority_tier", priorityTier);
      if (flag)         params.set("flag", flag);
      if (minPlddt > 0) params.set("min_plddt", String(minPlddt));
      if (geneApplied)  params.set("gene", geneApplied);

      const res = await fetch(`${apiUrl}/api/neoantigens/top?${params}`);
      if (res.ok) setData(await res.json());
    } catch (e) {
      console.error("Neoantigen fetch failed", e);
    } finally {
      setLoading(false);
    }
  }, [page, sortCol, sortOrder, bindingCat, hlaAllele, junctionType, priorityTier, flag, minPlddt, geneApplied]);

  useEffect(() => { fetchData(); }, [fetchData]);

  const handleSort = (col: SortCol) => {
    if (col === sortCol) {
      setSortOrder(o => o === "asc" ? "desc" : "asc");
    } else {
      setSortCol(col);
      setSortOrder(DEFAULT_ORDER[col]);
    }
    setPage(0);
  };

  const applyFilter = (fn: () => void) => { fn(); setPage(0); };

  const stats = data?.stats;

  return (
    <div className="flex flex-col gap-6">

      {/* ── Summary stats ───────────────────────────────────────────────────── */}
      {stats && (
        <div className="flex flex-wrap gap-3">
          <StatChip label="Unique peptides" value={stats.unique_peptides} color="text-foreground" />
          <StatChip label="Unique genes"    value={stats.unique_genes}    color="text-foreground" />
          <StatChip label="P1"              value={stats.p1}  sub="High-confidence"     color="text-emerald-400" />
          <StatChip label="P2"              value={stats.p2}  sub="Strong presentation" color="text-sky-400"     />
          <StatChip label="P3"              value={stats.p3}  sub="Supported"           color="text-amber-400"   />
          <StatChip label="Immunogenic"     value={stats.immunogenic} sub="BigMHC-IM ≥ 0.5" color="text-violet-400" />
          {stats.expressed != null && (
            <StatChip label="Expressed"     value={stats.expressed}   sub="Detected in data" color="text-teal-400" />
          )}
        </div>
      )}

      {/* ── Info note ───────────────────────────────────────────────────────── */}
      <div className="flex items-start gap-2 text-[11px] text-muted-foreground bg-muted/30 border border-border/40 rounded-lg px-3 py-2">
        <Info className="w-3 h-3 mt-0.5 shrink-0" />
        <span>
          <strong className="text-foreground">Priority tier (P1–P4)</strong> ranks on immunological evidence:
          multi-tool MHC <strong className="text-foreground">presentation</strong> consensus (MHCflurry, NetMHCpan-4.2 EL, BigMHC EL)
          combined with BigMHC <strong className="text-foreground">immunogenicity</strong>. Splice origin, PTM disruption, expression
          and fold quality are shown as orthogonal flags, not ranking gates.
        </span>
      </div>

      <div className="flex flex-col lg:flex-row gap-6">
        {/* ── Filter sidebar ─────────────────────────────────────────────────── */}
        <div className="w-full lg:w-52 shrink-0">
          <div className="bg-card/40 border border-border/40 rounded-xl p-4 sticky top-24 space-y-5">
            <h4 className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground flex items-center gap-1.5">
              <Filter className="w-3 h-3" /> Filters
            </h4>

            {/* Gene search */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Gene</span>
              <div className="relative">
                <Search className="absolute left-2 top-2 h-3 w-3 text-muted-foreground" />
                <Input
                  placeholder="e.g. TP53"
                  value={geneSearch}
                  onChange={e => setGeneSearch(e.target.value)}
                  onKeyDown={e => { if (e.key === "Enter") applyFilter(() => setGeneApplied(geneSearch.trim())); }}
                  className="pl-6 h-7 text-xs bg-background border-border/60"
                />
                {geneApplied && (
                  <button onClick={() => applyFilter(() => { setGeneSearch(""); setGeneApplied(""); })}
                    className="absolute right-2 top-2 text-muted-foreground hover:text-foreground text-xs">&times;</button>
                )}
              </div>
              {geneSearch && geneSearch !== geneApplied && (
                <button onClick={() => applyFilter(() => setGeneApplied(geneSearch.trim()))}
                  className="text-[10px] text-primary hover:underline">Apply →</button>
              )}
            </div>

            {/* Priority tier */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Priority Tier</span>
              <div className="flex flex-col gap-1">
                {PRIORITY_FILTERS.map(f => (
                  <button key={f.id}
                    onClick={() => applyFilter(() => setPriorityTier(f.id))}
                    className={cn(
                      "flex items-center gap-2 text-left px-2 py-1.5 rounded text-[11px] transition-colors",
                      priorityTier === f.id
                        ? "bg-primary/20 text-primary border border-primary/30"
                        : "text-muted-foreground hover:bg-muted hover:text-foreground"
                    )}>
                    <f.icon className="w-3 h-3 shrink-0" />
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            {/* Decorating flags */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Flag</span>
              <div className="flex flex-col gap-1">
                {FLAG_FILTERS.map(f => (
                  <button key={f.id}
                    onClick={() => applyFilter(() => setFlag(f.id))}
                    className={cn(
                      "flex items-center gap-2 text-left px-2 py-1.5 rounded text-[11px] transition-colors",
                      flag === f.id
                        ? "bg-violet-500/20 text-violet-300 border border-violet-500/30"
                        : "text-muted-foreground hover:bg-muted hover:text-foreground"
                    )}>
                    <f.icon className="w-3 h-3 shrink-0" />
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            {/* Binding category */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Binding</span>
              <div className="flex flex-col gap-1">
                {[
                  { id: "" as const,             label: "Any"               },
                  { id: "strong_binder" as const, label: "Strong (< 50 nM)" },
                  { id: "weak_binder"   as const, label: "Weak (< 500 nM)"  },
                ].map(f => (
                  <button key={f.id}
                    onClick={() => applyFilter(() => setBindingCat(f.id))}
                    className={cn(
                      "text-left px-2 py-1 rounded text-[11px] transition-colors",
                      bindingCat === f.id
                        ? "bg-rose-500/20 text-rose-400 border border-rose-500/30"
                        : "text-muted-foreground hover:bg-muted hover:text-foreground"
                    )}>
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            {/* HLA allele */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">HLA Allele</span>
              <select
                value={hlaAllele}
                onChange={e => applyFilter(() => setHlaAllele(e.target.value))}
                className="w-full bg-background border border-border/60 rounded px-2 py-1.5 text-xs text-foreground focus:border-primary outline-none"
              >
                <option value="">All alleles</option>
                {hlaOptions.map(a => (
                  <option key={a} value={a}>{a}</option>
                ))}
              </select>
            </div>

            {/* Junction type */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Junction Type</span>
              <div className="flex flex-col gap-1">
                {[
                  { id: "",                    label: "All"                  },
                  { id: "novel_insertion",     label: "Novel insertion"      },
                  { id: "novel_substitution",  label: "Novel substitution"   },
                ].map(f => (
                  <button key={f.id}
                    onClick={() => applyFilter(() => setJunctionType(f.id))}
                    className={cn(
                      "text-left px-2 py-1 rounded text-[11px] transition-colors",
                      junctionType === f.id
                        ? "bg-primary/20 text-primary border border-primary/30"
                        : "text-muted-foreground hover:bg-muted hover:text-foreground"
                    )}>
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            {/* Min pLDDT */}
            <div className="space-y-1.5">
              <div className="flex justify-between">
                <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Min pLDDT</span>
                <span className="text-[10px] text-primary">{minPlddt || "—"}</span>
              </div>
              <input type="range" min="0" max="90" step="5"
                value={minPlddt}
                onChange={e => applyFilter(() => setMinPlddt(parseInt(e.target.value)))}
                className="w-full accent-primary h-1 bg-muted rounded-lg cursor-pointer"
              />
            </div>
          </div>
        </div>

        {/* ── Results table ──────────────────────────────────────────────────── */}
        <div className="flex-1 space-y-4 min-w-0">
          {loading ? (
            <div className="space-y-2">
              {[...Array(8)].map((_, i) => (
                <div key={i} className="h-10 bg-card/40 rounded animate-pulse border border-border/20" />
              ))}
            </div>
          ) : !data?.items.length ? (
            <div className="flex flex-col items-center justify-center py-20 text-muted-foreground">
              <Syringe className="w-10 h-10 mb-3 opacity-30" />
              <p className="text-sm">No peptides match the current filters.</p>
              <button
                onClick={() => {
                  setPriorityTier(""); setFlag(""); setBindingCat(""); setHlaAllele(""); setJunctionType("");
                  setMinPlddt(0); setGeneApplied(""); setGeneSearch(""); setPage(0);
                }}
                className="mt-3 text-xs text-primary hover:underline"
              >
                Reset all filters
              </button>
            </div>
          ) : (
            <>
              <div className="overflow-x-auto rounded-xl border border-border/40">
                <table className="w-full text-xs min-w-[860px]">
                  <thead className="bg-card/60 border-b border-border/40 sticky top-0 z-10">
                    <tr>
                      <ColHeader col="priority" label="Tier"
                        sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        tip="Priority tier P1–P4: multi-tool presentation consensus × BigMHC immunogenicity. P1 = strong presentation + immunogenic."
                      />
                      <ColHeader col="gene"     label="Gene"     sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort} />
                      <ColHeader col="isoform"  label="Isoform"  sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort} className="text-left" />
                      <ColHeader col="peptide"  label="Peptide"  sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort} />
                      <ColHeader col="hla"      label="HLA"      sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort} className="pl-3 pr-1" />
                      <ColHeader col="consensus" label="Presentation"
                        sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-center pr-3" alignCenter
                        tip="Multi-tool presentation consensus (MHCflurry, NetMHCpan-4.2 EL, BigMHC EL): STRONG = ≥2 tools call a strong binder · MOD = 1 strong or ≥2 weak · WEAK = 1 weak. Hover a cell for the underlying median %rank."
                      />
                      <ColHeader col="immunogenicity" label="IMM"
                        sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-center pr-3" alignCenter
                        tip="BigMHC immunogenicity (0–1): predicted T-cell recognition. ≥0.5 = immunogenic."
                      />
                      <ColHeader col="affinity" label="IC50 (nM)"
                        sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-left pl-1 pr-3"
                        tip="MHC binding affinity (IC50, nM). Lower = stronger. <50 strong, <500 weak."
                      />
                      <th className="py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] text-muted-foreground text-left">
                        <span className="flex items-center gap-1">Flags
                          <InfoTip text="Orthogonal annotations: IMM immunogenic · SPL splice-enabled · EXP expressed · PTM disrupts a canonical PTM site · UNQ proteome-unique / SELF or REF = peptide also in a reference protein · n/3 = fewer than three predictors scored the row." />
                        </span>
                      </th>
                      <ColHeader col="plddt"    label="pLDDT"
                        sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-center"
                        tip="AlphaFold confidence over the epitope window (0–100): ≥90 very high, ≥70 confident, ≥50 low, <50 very low."
                      />
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-border/30">
                    {data.items.map((row, idx) => {
                      const meta  = priorityMeta(row.priority_tier);
                      const isP1  = row.priority_tier === "P1";
                      return (
                        <tr
                          key={`${row.id}-${idx}`}
                          className={cn(
                            "hover:bg-muted/30 transition-colors",
                            isP1 && "bg-emerald-500/5"
                          )}
                        >
                          {/* Priority tier */}
                          <td className="py-2 px-3">
                            <span title={row.presentation_tier ? `${PRESENTATION_LABEL[row.presentation_tier] ?? row.presentation_tier} presentation` : undefined}
                              className={cn(
                                "inline-flex items-center gap-1 px-1.5 py-0.5 rounded text-[9px] font-semibold border whitespace-nowrap",
                                meta.bg, meta.color, meta.border
                              )}>
                              {isP1 && <Sparkles className="w-2.5 h-2.5" />}
                              {row.priority_tier ?? "—"} · {meta.label}
                            </span>
                          </td>
                          {/* Gene */}
                          <td className="py-2 px-3">
                            <Link href={`/gene/${row.gene_id}`}
                              className="font-semibold text-foreground/90 hover:text-primary transition-colors">
                              {row.gene_symbol}
                            </Link>
                          </td>
                          {/* Isoform */}
                          <td className="py-2 px-3">
                            <Link href={`/isoform/${row.isoform_id}`}
                              className="font-mono text-[10px] text-foreground/60 hover:text-primary transition-colors">
                              {row.transcript_name || row.transcript_id}
                            </Link>
                          </td>
                          {/* Peptide */}
                          <td className="py-2 px-3">
                            <div className="flex items-center gap-1.5">
                              <span className="font-mono tracking-wider text-foreground font-medium">
                                {row.peptide_sequence}
                              </span>
                              <span className="text-muted-foreground/50 text-[9px]">
                                {row.peptide_start}–{row.peptide_end}
                              </span>
                            </div>
                          </td>
                          {/* HLA allele */}
                          <td className="py-2 pl-3 pr-1 font-mono text-foreground/70 whitespace-nowrap">
                            {row.hla_allele}
                          </td>
                          {/* Presentation tier (multi-tool consensus; median %rank on hover) */}
                          <td className="py-2 pr-3 text-center">
                            {row.presentation_tier && row.presentation_tier !== "NONE" ? (
                              <span
                                title={row.consensus_rank != null ? `median %rank ${row.consensus_rank.toFixed(2)}%` : undefined}
                                className={cn(
                                  "inline-flex items-center rounded px-1.5 py-0.5 text-[9px] font-bold uppercase tracking-wide border",
                                  row.presentation_tier === "STRONG"   ? "bg-emerald-500/15 text-emerald-300 border-emerald-500/30" :
                                  row.presentation_tier === "MODERATE" ? "bg-amber-500/15 text-amber-300 border-amber-500/30" :
                                                                         "bg-muted/40 text-muted-foreground border-border/50"
                                )}>
                                {row.presentation_tier === "STRONG" ? "Strong" : row.presentation_tier === "MODERATE" ? "Mod" : "Weak"}
                                {row.is_ba_rescue && (
                                  <span className="ml-1 text-[8px] text-orange-400" title="BA rescue: cysteine-rich peptide missed by EL tools">BA</span>
                                )}
                              </span>
                            ) : (
                              <span className="text-muted-foreground/30">—</span>
                            )}
                          </td>
                          {/* Immunogenicity */}
                          <td className="py-2 pr-3 text-center font-mono text-[10px]">
                            {row.bigmhc_im_score != null ? (
                              <span className={cn(
                                "tabular-nums",
                                row.bigmhc_im_score >= 0.7 ? "text-violet-400 font-semibold" :
                                row.bigmhc_im_score >= 0.5 ? "text-sky-400" :
                                "text-muted-foreground"
                              )}>
                                {row.bigmhc_im_score.toFixed(3)}
                              </span>
                            ) : (
                              <span className="text-muted-foreground/30">—</span>
                            )}
                          </td>
                          {/* IC50 */}
                          <td className="py-2 pl-1 pr-3 text-left">
                            <AffinityBadge nm={row.affinity_nm} category={row.binding_category} />
                          </td>
                          {/* Flags */}
                          <td className="py-2 px-3">
                            <FlagChips row={row} />
                          </td>
                          {/* pLDDT */}
                          <td className="py-2 px-3 text-center">
                            <PlddtDot value={row.avg_plddt} />
                          </td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>

              {/* Pagination */}
              <div className="flex items-center justify-between pt-2">
                <span className="text-xs text-muted-foreground">
                  Showing{" "}
                  <span className="text-foreground font-medium">
                    {page * LIMIT + 1}–{Math.min((page + 1) * LIMIT, data.total_count)}
                  </span>
                  {" "}of{" "}
                  <span className="text-foreground font-medium">
                    {data.total_count.toLocaleString()}
                  </span>
                  {" "}unique peptide × isoform entries
                </span>
                <div className="flex gap-2">
                  <Button variant="outline" size="sm" disabled={page === 0}
                    onClick={() => setPage(p => p - 1)} className="h-7 text-xs gap-1">
                    <ChevronLeft className="w-3 h-3" /> Prev
                  </Button>
                  <Button variant="outline" size="sm"
                    disabled={(page + 1) * LIMIT >= data.total_count}
                    onClick={() => setPage(p => p + 1)} className="h-7 text-xs gap-1">
                    Next <ChevronRight className="w-3 h-3" />
                  </Button>
                </div>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
