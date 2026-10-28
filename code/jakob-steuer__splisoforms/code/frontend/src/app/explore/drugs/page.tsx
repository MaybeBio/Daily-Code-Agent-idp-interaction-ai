"use client";

import { useEffect, useMemo, useState } from "react";
import Link from "next/link";
import { useRouter } from "next/navigation";
import { ArrowLeft, Pill, Search, TrendingDown, ShieldCheck, FlaskConical, ArrowUp, ArrowDown } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { getApiUrl } from "@/lib/api";
import { cn } from "@/lib/utils";

interface Stats {
  n_drugs: number; n_complexes: number; n_drugged_genes: number;
  fate_counts: Record<string, number>;
  n_subtractive_genes: number; n_total_scored_genes: number;
}
interface GeneRow {
  gene_id: number; symbol: string; n_drugs: number; n_complexes: number;
  n_isoforms_scored: number; n_subtractive: number; n_preserved: number; subtractive_frac: number;
}

type SortKey = "symbol" | "n_drugs" | "n_isoforms_scored" | "n_subtractive" | "subtractive_frac";

// Shared grid template so the header row and data rows line up exactly.
const GRID = "grid grid-cols-[1fr_84px_92px_104px_150px] gap-2 items-center";

// Sortable column header — click to sort, arrow shows the active column/direction.
function Th({ label, k, align = "right", sortKey, sortDir, onSort }: {
  label: string; k: SortKey; align?: "left" | "right";
  sortKey: SortKey; sortDir: "asc" | "desc"; onSort: (k: SortKey) => void;
}) {
  const active = sortKey === k;
  return (
    <button
      onClick={() => onSort(k)}
      className={cn(
        "flex items-center gap-1 hover:text-foreground transition-colors uppercase tracking-wider",
        align === "right" ? "justify-end" : "justify-start",
        active && "text-emerald-400"
      )}
    >
      <span>{label}</span>
      {active && (sortDir === "asc" ? <ArrowUp className="w-3 h-3" /> : <ArrowDown className="w-3 h-3" />)}
    </button>
  );
}

function StatCard({ label, value, sub, icon }: { label: string; value: string; sub?: string; icon?: React.ReactNode }) {
  return (
    <Card>
      <CardContent className="pt-5">
        <div className="flex items-center gap-2 text-muted-foreground text-xs mb-1">{icon}{label}</div>
        <div className="text-2xl font-bold">{value}</div>
        {sub && <div className="text-xs text-muted-foreground mt-1">{sub}</div>}
      </CardContent>
    </Card>
  );
}

export default function DrugsLensPage() {
  const router = useRouter();
  const [stats, setStats] = useState<Stats | null>(null);
  const [genes, setGenes] = useState<GeneRow[]>([]);
  const [q, setQ] = useState("");
  const [loading, setLoading] = useState(true);
  const [sortKey, setSortKey] = useState<SortKey>("subtractive_frac");
  const [sortDir, setSortDir] = useState<"asc" | "desc">("desc");

  useEffect(() => {
    const api = getApiUrl();
    setLoading(true);
    Promise.all([
      fetch(`${api}/api/drugs/stats`).then(r => r.json()),
      fetch(`${api}/api/drugs/genes?limit=10000`).then(r => r.json()),
    ]).then(([s, g]) => { setStats(s); setGenes(g); }).finally(() => setLoading(false));
  }, []);

  const toggleSort = (k: SortKey) => {
    if (k === sortKey) setSortDir(d => (d === "asc" ? "desc" : "asc"));
    else { setSortKey(k); setSortDir(k === "symbol" ? "asc" : "desc"); }
  };

  const filtered = useMemo(() => {
    const f = genes.filter(g => g.symbol.toLowerCase().includes(q.toLowerCase()));
    const dir = sortDir === "asc" ? 1 : -1;
    return [...f].sort((a, b) => {
      const av = a[sortKey], bv = b[sortKey];
      if (typeof av === "string" || typeof bv === "string")
        return String(av).localeCompare(String(bv)) * dir;
      return ((av as number) - (bv as number)) * dir;
    });
  }, [genes, q, sortKey, sortDir]);

  const decidable = stats
    ? (stats.fate_counts.lost ?? 0) + (stats.fate_counts.disrupted ?? 0) + (stats.fate_counts.preserved ?? 0)
    : 0;
  const subtractive = stats ? (stats.fate_counts.lost ?? 0) + (stats.fate_counts.disrupted ?? 0) : 0;

  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-[1500px] mx-auto p-6 md:p-10">
        <header className="flex items-center gap-4 mb-8">
          <Button variant="ghost" size="icon" onClick={() => router.push("/")} className="hover:bg-emerald-500/10">
            <ArrowLeft className="w-5 h-5" />
          </Button>
          <div className="flex items-center gap-3">
            <Pill className="w-7 h-7 text-emerald-500" />
            <div>
              <h1 className="text-3xl font-bold tracking-tight bg-clip-text text-transparent bg-gradient-to-r from-emerald-400 via-foreground/80 to-secondary">
                Drug–Isoform Atlas
              </h1>
              <p className="text-sm text-muted-foreground">How alternative splicing reshapes druggable pockets, from co-folded AF3 drug complexes.</p>
            </div>
          </div>
        </header>

        {/* Headline stats */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
          <StatCard label="Drugs" value={stats ? stats.n_drugs.toLocaleString() : "…"} icon={<FlaskConical className="w-3.5 h-3.5" />} />
          <StatCard label="Co-folded complexes" value={stats ? stats.n_complexes.toLocaleString() : "…"} icon={<Pill className="w-3.5 h-3.5" />} />
          <StatCard label="Drugged genes" value={stats ? stats.n_drugged_genes.toLocaleString() : "…"} icon={<ShieldCheck className="w-3.5 h-3.5" />} />
          <StatCard
            label="Genes with pocket loss"
            value={stats ? `${stats.n_subtractive_genes.toLocaleString()}` : "…"}
            sub={stats ? `of ${stats.n_total_scored_genes.toLocaleString()} scored (${Math.round(100 * stats.n_subtractive_genes / (stats.n_total_scored_genes || 1))}%)` : ""}
            icon={<TrendingDown className="w-3.5 h-3.5" />}
          />
        </div>

        {/* Subtractive-vs-preserved headline bar */}
        {stats && decidable > 0 && (
          <Card className="mb-6">
            <CardContent className="pt-5">
              <div className="flex justify-between text-xs text-muted-foreground mb-2">
                <span>Splicing is subtractive on druggability</span>
                <span>{subtractive.toLocaleString()} subtractive · {(stats.fate_counts.preserved ?? 0).toLocaleString()} preserved · {(stats.fate_counts.no_call ?? 0).toLocaleString()} no-call</span>
              </div>
              <div className="flex h-3 rounded-full overflow-hidden">
                <div className="bg-red-500" style={{ width: `${100 * (stats.fate_counts.lost ?? 0) / decidable}%` }} title="lost" />
                <div className="bg-amber-500" style={{ width: `${100 * (stats.fate_counts.disrupted ?? 0) / decidable}%` }} title="disrupted" />
                <div className="bg-emerald-500" style={{ width: `${100 * (stats.fate_counts.preserved ?? 0) / decidable}%` }} title="preserved" />
              </div>
            </CardContent>
          </Card>
        )}

        {/* Controls */}
        <div className="flex items-center gap-3 mb-3">
          <div className="relative flex-1 max-w-sm">
            <Search className="w-4 h-4 absolute left-3 top-1/2 -translate-y-1/2 text-muted-foreground" />
            <Input placeholder="Filter genes…" value={q} onChange={e => setQ(e.target.value)} className="pl-9" />
          </div>
        </div>

        {/* Gene table */}
        <Card>
          <CardHeader className="py-3">
            <CardTitle className="text-sm text-muted-foreground font-normal">
              {loading ? "Loading…" : `${filtered.length.toLocaleString()} drugged genes`}
            </CardTitle>
          </CardHeader>
          <CardContent className="pt-0">
            <div className={cn(GRID, "text-[11px] font-medium text-muted-foreground px-3 pb-2 border-b border-border/50 select-none")}>
              <Th label="Gene" k="symbol" align="left" sortKey={sortKey} sortDir={sortDir} onSort={toggleSort} />
              <Th label="Drugs" k="n_drugs" sortKey={sortKey} sortDir={sortDir} onSort={toggleSort} />
              <Th label="Isoforms" k="n_isoforms_scored" sortKey={sortKey} sortDir={sortDir} onSort={toggleSort} />
              <Th label="Subtractive" k="n_subtractive" sortKey={sortKey} sortDir={sortDir} onSort={toggleSort} />
              <Th label="Burden" k="subtractive_frac" sortKey={sortKey} sortDir={sortDir} onSort={toggleSort} />
            </div>
            <div className="max-h-[52vh] overflow-y-auto">
              {filtered.map(g => (
                <Link key={g.gene_id} href={`/drugs/${g.symbol}`}
                  className={cn(GRID, "px-3 py-2 text-sm hover:bg-muted/50 rounded-md border-b border-border/30")}>
                  <span className="font-semibold text-emerald-300 truncate">{g.symbol}</span>
                  <span className="text-right tabular-nums">{g.n_drugs}</span>
                  <span className="text-right tabular-nums text-muted-foreground">{g.n_isoforms_scored}</span>
                  <span className="text-right tabular-nums">{g.n_subtractive}/{g.n_subtractive + g.n_preserved}</span>
                  <span className="flex items-center gap-2 justify-end">
                    <span className="w-16 h-1.5 rounded-full bg-muted overflow-hidden">
                      <span className="block h-full bg-gradient-to-r from-amber-500 to-red-500" style={{ width: `${g.subtractive_frac * 100}%` }} />
                    </span>
                    <span className="tabular-nums text-xs w-9 text-right">{Math.round(g.subtractive_frac * 100)}%</span>
                  </span>
                </Link>
              ))}
              {!loading && filtered.length === 0 && <p className="text-sm text-muted-foreground px-3 py-6 text-center">No genes match “{q}”.</p>}
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  );
}
