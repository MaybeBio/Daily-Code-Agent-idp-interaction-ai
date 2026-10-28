"use client";

import { useEffect, useState, useCallback, Suspense } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import Link from "next/link";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import {
  ArrowLeft, Filter, Search, Info, Microscope,
  ChevronLeft, ChevronRight, ArrowUp, ArrowDown, ArrowUpDown, HelpCircle,
  LayoutGrid, Layers, ShieldAlert, CircleDot, AlertTriangle,
} from "lucide-react";
import { cn } from "@/lib/utils";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { getApiUrl } from "@/lib/api";

// ─── Types ────────────────────────────────────────────────────────────────────

// Server-sortable Discovery Index columns.
type SortCol = "impact_score" | "length" | "isoform_count" | "plddt" | "nmd" | "ptm" | "neoantigen" | "neoantigen_p1";

interface BrowseItem {
  id: number;
  transcript_id: string;
  transcript_name?: string | null;
  gene_id: number;
  gene_symbol: string;
  plddt_mean?: number | null;
  delta_length?: number | null;
  source_name?: string | null;
  has_neoantigen?: boolean;
  neoantigen_count?: number;
  neoantigen_p1_count?: number;
  num_domains_lost?: number;
  num_domains_gained?: number;
  nmd_status?: string | null;
  has_ptm_exposure_change?: boolean;
  has_ptm_buried_change?: boolean;
  has_ptm_loss?: boolean;
  isoform_count?: number;
  impact_score?: number | null;
}

interface GlobalStats {
  impact_categories: { category: string; count: number }[];
  nmd_distribution: { status: string; count: number }[];
  ptm_distribution: { impact: string; count: number }[];
  avg_delta_length?: number;
  avg_domains_lost?: number;
}

interface BrowseResponse {
  items: BrowseItem[];
  total_count: number;
  stats: GlobalStats;
}

interface SourceOpt { id: number; name: string; }

// ─── Small UI helpers (mirror the Neoepitope Atlas table) ──────────────────────

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

function StatChip({ label, value, sub, color }: { label: string; value: string | number; sub?: string; color: string }) {
  return (
    <div className="flex flex-col gap-0.5 bg-card/40 border border-border/40 rounded-lg px-4 py-3 min-w-[120px]">
      <span className={`text-xl font-bold tabular-nums ${color}`}>
        {typeof value === "number" ? value.toLocaleString() : value}
      </span>
      <span className="text-[9px] font-bold uppercase tracking-widest text-muted-foreground">{label}</span>
      {sub && <span className="text-[9px] text-muted-foreground/60">{sub}</span>}
    </div>
  );
}

function ColHeader({
  col, label, sortCol, sortOrder, onSort, className = "", tip, alignRight = false, alignCenter = false,
}: {
  col: SortCol; label: string; sortCol: SortCol; sortOrder: "asc" | "desc";
  onSort: (col: SortCol) => void; className?: string; tip?: string; alignRight?: boolean; alignCenter?: boolean;
}) {
  const active = sortCol === col;
  return (
    <th
      onClick={() => onSort(col)}
      className={cn(
        "py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] cursor-pointer select-none whitespace-nowrap hover:text-foreground transition-colors",
        active ? "text-primary" : "text-muted-foreground",
        className,
      )}
    >
      <span className={cn("flex items-center gap-1", alignRight && "justify-end", alignCenter && "justify-center")}>
        {label}
        {tip && <InfoTip text={tip} />}
        {active
          ? sortOrder === "asc" ? <ArrowUp className="w-2.5 h-2.5" /> : <ArrowDown className="w-2.5 h-2.5" />
          : <ArrowUpDown className="w-2.5 h-2.5 opacity-30" />}
      </span>
    </th>
  );
}

function PlddtDot({ value }: { value?: number | null }) {
  if (value == null) return <span className="text-muted-foreground/40 text-[10px]">—</span>;
  const dot = value >= 90 ? "bg-[#0053D6]" : value >= 70 ? "bg-[#65CBF3]" : value >= 50 ? "bg-[#FFDB13]" : "bg-[#FF7D45]";
  return (
    <span className="flex items-center justify-center gap-1">
      <span className={`w-1.5 h-1.5 rounded-full ${dot} inline-block`} />
      <span className="font-mono text-[10px]">{value.toFixed(0)}</span>
    </span>
  );
}

function impactBand(score: number) {
  if (score >= 60) return { color: "text-rose-400", bar: "bg-rose-500" };
  if (score >= 40) return { color: "text-amber-400", bar: "bg-amber-500" };
  if (score >= 20) return { color: "text-sky-400", bar: "bg-sky-500" };
  return { color: "text-muted-foreground", bar: "bg-muted-foreground/40" };
}

function shortSource(name?: string | null) {
  if (!name) return "—";
  if (name.includes("GENCODE")) return "GENCODE";
  if (name.includes("ccRCC")) return "ccRCC";
  if (name.includes("Breast")) return "BC";
  return name;
}

// ─── Constants ────────────────────────────────────────────────────────────────

const IMPACT_FILTERS = [
  { id: "all",         label: "All perturbations", icon: LayoutGrid },
  { id: "domain_loss", label: "Domain loss",       icon: Layers },
  { id: "neoantigen",  label: "Neoantigens",       icon: ShieldAlert },
  { id: "ptm_impact",  label: "PTM alteration",    icon: CircleDot },
  { id: "nmd_target",  label: "NMD target",        icon: AlertTriangle },
];

const LIMIT = 50;

const DEFAULT_ORDER: Record<SortCol, "asc" | "desc"> = {
  impact_score:  "desc",   // highest impact first
  length:        "desc",   // largest length delta first
  isoform_count: "desc",   // most isoforms per gene first
  plddt:         "desc",   // most confident first
  nmd:           "desc",   // NMD targets first
  ptm:           "desc",   // PTM loss first
  neoantigen:    "desc",   // most neoepitopes first
  neoantigen_p1: "desc",   // most P1 first
};

function nmdMeta(s?: string | null) {
  if (s === "nmd_target") return { label: "Target", cls: "bg-rose-500/15 text-rose-400 border-rose-500/30" };
  if (s === "nmd_escape") return { label: "Escape", cls: "bg-emerald-500/15 text-emerald-400 border-emerald-500/30" };
  if (s === "unknown") return { label: "Unknown", cls: "bg-muted/40 text-muted-foreground border-border/50" };
  return null;
}

// ─── Main explorer ─────────────────────────────────────────────────────────────

function DashboardExplorer() {
  const searchParams = useSearchParams();
  const router = useRouter();

  const [data, setData] = useState<BrowseResponse | null>(null);
  const [loading, setLoading] = useState(true);
  const [page, setPage] = useState(0);

  const [sources, setSources] = useState<SourceOpt[]>([]);

  // Filters
  const [minPlddt, setMinPlddt] = useState(0);
  const [minIsoforms, setMinIsoforms] = useState(1);
  const [impactFilter, setImpactFilter] = useState("all");
  const [sourceId, setSourceId] = useState("");

  // Sort — driven by the column headers
  const [sortCol, setSortCol] = useState<SortCol>("impact_score");
  const [sortOrder, setSortOrder] = useState<"asc" | "desc">("desc");

  // Search is kept in the URL so deep-links (e.g. from the homepage) work.
  const searchQuery = searchParams.get("query") || "";
  const [searchInput, setSearchInput] = useState(searchQuery);
  useEffect(() => { setSearchInput(searchParams.get("query") || ""); }, [searchParams]);

  useEffect(() => {
    fetch(`${getApiUrl()}/api/sources`)
      .then(r => (r.ok ? r.json() : []))
      .then(setSources)
      .catch(() => {});
  }, []);

  const handleSearch = () => {
    setPage(0);
    router.push(searchInput.trim() ? `/explore/impact?query=${encodeURIComponent(searchInput.trim())}` : "/explore/impact");
  };

  const fetchData = useCallback(async () => {
    setLoading(true);
    try {
      const params = new URLSearchParams({
        limit: String(LIMIT),
        offset: String(page * LIMIT),
        sort_by: sortCol,
        sort_order: sortOrder,
      });
      if (searchQuery) params.set("search", searchQuery);
      if (minPlddt > 0) params.set("min_plddt", String(minPlddt));
      if (minIsoforms > 1) params.set("min_isoforms", String(minIsoforms));
      if (impactFilter !== "all") params.set("impact_filter", impactFilter);
      if (sourceId) params.set("source_id", sourceId);

      const res = await fetch(`${getApiUrl()}/api/isoforms/browse?${params}`);
      if (res.ok) setData(await res.json());
    } catch (e) {
      console.error("Discovery load failed", e);
    } finally {
      setLoading(false);
    }
  }, [page, sortCol, sortOrder, searchQuery, minPlddt, minIsoforms, impactFilter, sourceId]);

  useEffect(() => { fetchData(); }, [fetchData]);

  const handleSort = (col: SortCol) => {
    if (col === sortCol) setSortOrder(o => (o === "asc" ? "desc" : "asc"));
    else { setSortCol(col); setSortOrder(DEFAULT_ORDER[col]); }
    setPage(0);
  };

  const applyFilter = (fn: () => void) => { fn(); setPage(0); };

  const s = data?.stats;
  const lost   = s?.impact_categories.find(c => c.category === "Domain Loss")?.count || 0;
  const gained = s?.impact_categories.find(c => c.category === "Domain Gained")?.count || 0;
  const neo    = s?.impact_categories.find(c => c.category === "Neoantigens")?.count || 0;
  const ptmTot = s?.ptm_distribution.reduce((a, c) => a + c.count, 0) || 0;
  const nmdTgt = s?.nmd_distribution.find(n => n.status === "Target")?.count || 0;

  return (
    <div className="flex flex-col gap-6">
      {/* Summary stat chips (dataset-wide context) */}
      {s && (
        <div className="flex flex-wrap gap-3">
          <StatChip label="Domain disruptions" value={lost + gained} sub="lost + gained" color="text-orange-400" />
          <StatChip label="Neoantigenic"        value={neo}          sub="novel binders" color="text-rose-400" />
          <StatChip label="PTM alterations"     value={ptmTot}       sub="exposed / buried / lost" color="text-cyan-400" />
          <StatChip label="NMD targets"         value={nmdTgt}       sub="mRNA decay" color="text-primary" />
          {s.avg_delta_length != null && (
            <StatChip label="Mean |Δlen|"       value={`${Math.abs(s.avg_delta_length).toFixed(0)} aa`} sub="vs canonical" color="text-amber-400" />
          )}
        </div>
      )}

      {/* Info note */}
      <div className="flex items-start gap-2 text-[11px] text-muted-foreground bg-muted/30 border border-border/40 rounded-lg px-3 py-2">
        <Info className="w-3 h-3 mt-0.5 shrink-0" />
        <span>
          The <strong className="text-foreground">Impact Score</strong> is a bounded 0–100 composite of five
          axes (domain loss, length change, neoantigen priority, NMD fate, PTM accessibility), tempered by
          structural confidence — see the{" "}
          <Link href="/documentation/discovery" className="text-primary hover:underline">Discovery guide</Link>.
          Rows are novel (non-canonical) isoforms; click a column header to sort.
        </span>
      </div>

      <div className="flex flex-col lg:flex-row gap-6">
        {/* ── Filters sidebar (filters only — sorting lives in the column headers) ── */}
        <div className="w-full lg:w-56 shrink-0">
          <div className="bg-card/40 border border-border/40 rounded-xl p-4 sticky top-24 space-y-5">
            <h4 className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground flex items-center gap-1.5">
              <Filter className="w-3 h-3" /> Filters
            </h4>

            {/* Search */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Search</span>
              <div className="relative">
                <Search className="absolute left-2 top-2 h-3 w-3 text-muted-foreground" />
                <Input
                  placeholder="Gene, transcript, or ID"
                  value={searchInput}
                  onChange={e => setSearchInput(e.target.value)}
                  onKeyDown={e => { if (e.key === "Enter") handleSearch(); }}
                  className="pl-6 h-7 text-xs bg-background border-border/60"
                />
                {searchInput && (
                  <button onClick={() => { setSearchInput(""); router.push("/explore/impact"); }}
                    className="absolute right-2 top-2 text-muted-foreground hover:text-foreground text-xs">&times;</button>
                )}
              </div>
              {searchInput.trim() && searchInput.trim() !== searchQuery && (
                <button onClick={handleSearch} className="text-[10px] text-primary hover:underline">Apply →</button>
              )}
            </div>

            {/* Impact profile */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Impact Profile</span>
              <div className="flex flex-col gap-1">
                {IMPACT_FILTERS.map(f => (
                  <button key={f.id}
                    onClick={() => applyFilter(() => setImpactFilter(f.id))}
                    className={cn(
                      "flex items-center gap-2 text-left px-2 py-1.5 rounded text-[11px] transition-colors",
                      impactFilter === f.id
                        ? "bg-primary/20 text-primary border border-primary/30"
                        : "text-muted-foreground hover:bg-muted hover:text-foreground"
                    )}>
                    <f.icon className="w-3 h-3 shrink-0" />
                    {f.label}
                  </button>
                ))}
              </div>
            </div>

            {/* Dataset */}
            <div className="space-y-1.5">
              <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Dataset</span>
              <select
                value={sourceId}
                onChange={e => applyFilter(() => setSourceId(e.target.value))}
                className="w-full bg-background border border-border/60 rounded px-2 py-1.5 text-xs text-foreground focus:border-primary outline-none"
              >
                <option value="">All datasets</option>
                {sources.map(src => (
                  <option key={src.id} value={src.id}>{shortSource(src.name)}</option>
                ))}
              </select>
            </div>

            {/* Min pLDDT */}
            <div className="space-y-1.5">
              <div className="flex justify-between">
                <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Min pLDDT</span>
                <span className="text-[10px] text-primary">{minPlddt || "—"}</span>
              </div>
              <input type="range" min="0" max="100" step="5" value={minPlddt}
                onChange={e => applyFilter(() => setMinPlddt(parseInt(e.target.value)))}
                className="w-full accent-primary h-1 bg-muted rounded-lg cursor-pointer" />
            </div>

            {/* Min isoforms per gene */}
            <div className="space-y-1.5">
              <div className="flex justify-between">
                <span className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground">Min isoforms / gene</span>
                <span className="text-[10px] text-secondary">{minIsoforms}</span>
              </div>
              <input type="range" min="1" max="50" value={minIsoforms}
                onChange={e => applyFilter(() => setMinIsoforms(parseInt(e.target.value)))}
                className="w-full accent-secondary h-1 bg-muted rounded-lg cursor-pointer" />
            </div>
          </div>
        </div>

        {/* ── Results table ──────────────────────────────────────────────────── */}
        <div className="flex-1 space-y-4 min-w-0">
          {loading ? (
            <div className="space-y-2">
              {[...Array(10)].map((_, i) => (
                <div key={i} className="h-10 bg-card/40 rounded animate-pulse border border-border/20" />
              ))}
            </div>
          ) : !data?.items.length ? (
            <div className="flex flex-col items-center justify-center py-20 text-muted-foreground">
              <Microscope className="w-10 h-10 mb-3 opacity-30" />
              <p className="text-sm">No isoforms match the current filters.</p>
              <button
                onClick={() => { setImpactFilter("all"); setSourceId(""); setMinPlddt(0); setMinIsoforms(1); setSearchInput(""); router.push("/explore/impact"); }}
                className="mt-3 text-xs text-primary hover:underline"
              >
                Reset all filters
              </button>
            </div>
          ) : (
            <>
              <div className="overflow-x-auto rounded-xl border border-border/40">
                <table className="w-full text-xs min-w-[880px]">
                  <thead className="bg-card/60 border-b border-border/40 sticky top-0 z-10">
                    <tr>
                      <ColHeader col="impact_score" label="Impact" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        tip="Structural Impact Score (0–100): a bounded composite of domain loss, length change, neoantigen priority, NMD fate and PTM accessibility, tempered by pLDDT." />
                      <th className="py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] text-muted-foreground text-left">Gene</th>
                      <th className="py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] text-muted-foreground text-left">Isoform</th>
                      <ColHeader col="length" label="Δ len" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-right" alignRight tip="Protein length change vs the canonical isoform (aa)." />
                      <th className="py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] text-muted-foreground text-center">
                        <span className="flex items-center justify-center gap-1">Domains
                          <InfoTip text="Pfam domains lost (−) or gained (+) relative to canonical." />
                        </span>
                      </th>
                      <ColHeader col="neoantigen" label="Neo" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-right" alignRight tip="Total ranked neoepitopes (representative peptides) for this isoform." />
                      <ColHeader col="neoantigen_p1" label="P1" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-right" alignRight tip="High-confidence (P1) neoepitopes: strong multi-tool presentation and predicted immunogenic." />
                      <ColHeader col="nmd" label="NMD" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-center" alignCenter tip="Nonsense-mediated decay prediction: Target (likely degraded), Escape (likely translated), or Unknown." />
                      <ColHeader col="ptm" label="PTM" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-center" alignCenter tip="PTM accessibility change vs canonical: Loss (site deleted), Exp (more exposed), Bur (more buried)." />
                      <ColHeader col="isoform_count" label="Isoforms" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-right" alignRight tip="Number of novel isoforms for this gene." />
                      <ColHeader col="plddt" label="pLDDT" sortCol={sortCol} sortOrder={sortOrder} onSort={handleSort}
                        className="text-center" alignCenter tip="Mean AlphaFold confidence over the model (0–100)." />
                      <th className="py-2.5 px-3 font-bold uppercase tracking-widest text-[9px] text-muted-foreground text-left">Source</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-border/30">
                    {data.items.map((row, idx) => {
                      const score = row.impact_score ?? 0;
                      const band = impactBand(score);
                      const dLost = row.num_domains_lost ?? 0;
                      const dGained = row.num_domains_gained ?? 0;
                      const dl = row.delta_length ?? 0;
                      return (
                        <tr key={`${row.id}-${idx}`} className="hover:bg-muted/30 transition-colors">
                          {/* Impact */}
                          <td className="py-2 px-3">
                            <div className="flex items-center gap-2">
                              <div className="w-14 bg-muted h-1.5 rounded-full overflow-hidden shrink-0">
                                <div className={`h-full ${band.bar}`} style={{ width: `${Math.min(Math.max(score, 0), 100)}%` }} />
                              </div>
                              <span className={cn("font-mono text-[11px] font-semibold tabular-nums", band.color)}>
                                {row.impact_score != null ? score.toFixed(0) : "—"}
                              </span>
                            </div>
                          </td>
                          {/* Gene */}
                          <td className="py-2 px-3">
                            <Link href={`/gene/${row.gene_id}`} className="font-semibold text-foreground/90 hover:text-primary transition-colors">
                              {row.gene_symbol}
                            </Link>
                          </td>
                          {/* Isoform */}
                          <td className="py-2 px-3">
                            <Link href={`/isoform/${row.id}`} className="font-mono text-[10px] text-foreground/60 hover:text-primary transition-colors">
                              {row.transcript_name || row.transcript_id}
                            </Link>
                          </td>
                          {/* Δ length */}
                          <td className="py-2 px-3 text-right font-mono text-[10px] text-foreground/80 tabular-nums">
                            {dl > 0 ? `+${dl}` : dl} aa
                          </td>
                          {/* Domains */}
                          <td className="py-2 px-3">
                            <div className="flex items-center justify-center gap-1">
                              {dLost > 0 && (
                                <span title="Pfam domains lost" className="inline-flex items-center rounded px-1 py-0.5 text-[9px] font-bold border bg-rose-500/15 text-rose-400 border-rose-500/30">−{dLost}</span>
                              )}
                              {dGained > 0 && (
                                <span title="Pfam domains gained" className="inline-flex items-center rounded px-1 py-0.5 text-[9px] font-bold border bg-emerald-500/15 text-emerald-400 border-emerald-500/30">+{dGained}</span>
                              )}
                              {dLost === 0 && dGained === 0 && <span className="text-muted-foreground/30">—</span>}
                            </div>
                          </td>
                          {/* Neoantigen total */}
                          <td className="py-2 px-3 text-right font-mono text-[10px] tabular-nums">
                            {row.neoantigen_count
                              ? <span className="text-foreground/80">{row.neoantigen_count.toLocaleString()}</span>
                              : <span className="text-muted-foreground/30">—</span>}
                          </td>
                          {/* Neoantigen P1 */}
                          <td className="py-2 px-3 text-right font-mono text-[10px] tabular-nums">
                            {row.neoantigen_p1_count
                              ? <span className="text-rose-400 font-semibold">{row.neoantigen_p1_count.toLocaleString()}</span>
                              : <span className="text-muted-foreground/30">—</span>}
                          </td>
                          {/* NMD */}
                          <td className="py-2 px-3 text-center">
                            {(() => {
                              const m = nmdMeta(row.nmd_status);
                              return m
                                ? <span className={cn("inline-flex items-center rounded px-1 py-0.5 text-[9px] font-bold border", m.cls)}>{m.label}</span>
                                : <span className="text-muted-foreground/30">—</span>;
                            })()}
                          </td>
                          {/* PTM */}
                          <td className="py-2 px-3">
                            <div className="flex items-center justify-center gap-1">
                              {row.has_ptm_loss && <span title="PTM site deleted" className="inline-flex items-center rounded px-1 py-0.5 text-[9px] font-bold border bg-rose-500/15 text-rose-400 border-rose-500/30">Loss</span>}
                              {row.has_ptm_exposure_change && <span title="PTM site more exposed" className="inline-flex items-center rounded px-1 py-0.5 text-[9px] font-bold border bg-cyan-500/15 text-cyan-400 border-cyan-500/30">Exp</span>}
                              {row.has_ptm_buried_change && <span title="PTM site more buried" className="inline-flex items-center rounded px-1 py-0.5 text-[9px] font-bold border bg-blue-500/15 text-blue-400 border-blue-500/30">Bur</span>}
                              {!row.has_ptm_loss && !row.has_ptm_exposure_change && !row.has_ptm_buried_change && <span className="text-muted-foreground/30">—</span>}
                            </div>
                          </td>
                          {/* Isoforms per gene */}
                          <td className="py-2 px-3 text-right font-mono text-[10px] text-foreground/70 tabular-nums">
                            {row.isoform_count ?? "—"}
                          </td>
                          {/* pLDDT */}
                          <td className="py-2 px-3 text-center">
                            <PlddtDot value={row.plddt_mean} />
                          </td>
                          {/* Source */}
                          <td className="py-2 px-3 text-[10px] text-muted-foreground whitespace-nowrap">
                            {shortSource(row.source_name)}
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
                  </span>{" "}
                  of <span className="text-foreground font-medium">{data.total_count.toLocaleString()}</span> isoforms
                </span>
                <div className="flex gap-2">
                  <Button variant="outline" size="sm" disabled={page === 0} onClick={() => setPage(p => p - 1)} className="h-7 text-xs gap-1">
                    <ChevronLeft className="w-3 h-3" /> Prev
                  </Button>
                  <Button variant="outline" size="sm" disabled={(page + 1) * LIMIT >= data.total_count} onClick={() => setPage(p => p + 1)} className="h-7 text-xs gap-1">
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

// ─── Page shell ────────────────────────────────────────────────────────────────

function ImpactPageInner() {
  const router = useRouter();
  const activeQuery = useSearchParams().get("query") || "";

  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-[1700px] mx-auto p-6 md:p-10">
        <header className="flex flex-col md:flex-row justify-between items-start md:items-center gap-4 mb-8">
          <div className="flex gap-4 items-center">
            <Button variant="ghost" size="icon" onClick={() => router.push("/")} className="hover:bg-primary/10">
              <ArrowLeft className="w-5 h-5" />
            </Button>
            <div>
              <h1 className="text-4xl font-bold tracking-tight bg-clip-text text-transparent bg-gradient-to-r from-primary via-foreground/80 to-secondary">
                Splicing Impact Explorer
              </h1>
              <p className="text-sm text-muted-foreground mt-1">
                {activeQuery
                  ? <>Results for <span className="text-primary font-medium">&ldquo;{activeQuery}&rdquo;</span></>
                  : "Novel isoforms ranked by structural impact across the cancer isoform landscape."}
              </p>
            </div>
          </div>
        </header>

        <DashboardExplorer />
      </div>
    </div>
  );
}

export default function ImpactPage() {
  return (
    <Suspense fallback={<div className="text-center py-24 animate-pulse text-primary tracking-widest uppercase text-xs">Initializing Discovery Engine...</div>}>
      <ImpactPageInner />
    </Suspense>
  );
}
