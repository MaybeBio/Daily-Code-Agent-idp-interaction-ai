"use client";

import { useMemo } from "react";
import Link from "next/link";
import { Pill, ArrowRight, FlaskConical } from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { cn } from "@/lib/utils";

// data = /api/gene/{id}/drugs response; gene carries the isoform id→transcript map.
interface Fate {
  isoform_id: number; drug_name?: string; fate?: string;
  mechanism_class?: string; frac_retained?: number; decidable?: boolean;
}
interface Complex { id: number; drug: { name: string; pubchem_cid?: string }; holo_grade?: string; }

const FATE = {
  lost:      { label: "Lost",      cls: "bg-red-500/15 text-red-400 border-red-500/30",       bar: "bg-red-500" },
  disrupted: { label: "Disrupted", cls: "bg-amber-500/15 text-amber-400 border-amber-500/30", bar: "bg-amber-500" },
  preserved: { label: "Preserved", cls: "bg-emerald-500/15 text-emerald-400 border-emerald-500/30", bar: "bg-emerald-500" },
  no_call:   { label: "No call",   cls: "bg-slate-500/15 text-slate-400 border-slate-500/30",  bar: "bg-slate-500" },
} as const;
type FateKey = keyof typeof FATE;

const classify = (frac: number): FateKey =>
  frac < 0.25 ? "lost" : frac < 0.75 ? "disrupted" : "preserved";

export default function GeneDrugPanel({ data, gene }: { data: any; gene: any }) {
  const complexes: Complex[] = data.complexes ?? [];
  const fates: Fate[] = data.fates ?? [];
  const isoMap = useMemo(
    () => new Map<number, any>((gene.isoforms ?? []).map((i: any) => [i.id, i])),
    [gene.isoforms]);

  const nDrugs = useMemo(() => new Set(complexes.map(c => c.drug.pubchem_cid ?? c.drug.name)).size, [complexes]);

  // overall fate mix across every drug × isoform pair
  const mix = useMemo(() => {
    const m: Record<string, number> = { lost: 0, disrupted: 0, preserved: 0, no_call: 0 };
    for (const f of fates) if (f.fate && f.fate in m) m[f.fate]++;
    return m;
  }, [fates]);
  const total = mix.lost + mix.disrupted + mix.preserved + mix.no_call;

  // per-isoform rollup → mean retention (decidable) + worst fate, ranked worst-first
  const rows = useMemo(() => {
    const by = new Map<number, { sum: number; n: number; decidable: number }>();
    for (const f of fates) {
      const e = by.get(f.isoform_id) ?? { sum: 0, n: 0, decidable: 0 };
      if (f.decidable) { e.sum += f.frac_retained ?? 0; e.decidable++; }
      e.n++;
      by.set(f.isoform_id, e);
    }
    return [...by.entries()].map(([id, e]) => {
      const iso = isoMap.get(id);
      const meanFrac = e.decidable > 0 ? e.sum / e.decidable : null;
      const fate: FateKey = meanFrac == null ? "no_call" : classify(meanFrac);
      return {
        id,
        // Readable name when we have one; the raw source id otherwise.
        transcript_id: iso?.transcript_name ?? iso?.transcript_id ?? `#${id}`,
        len: iso?.sequence?.length,
        mane: iso?.is_mane_select,
        meanFrac, fate, nDrugs: e.n,
      };
    }).sort((a, b) => (a.meanFrac ?? 1.1) - (b.meanFrac ?? 1.1));
  }, [fates, isoMap]);

  return (
    <div className="space-y-4">
      {/* Summary + CTA */}
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div className="flex items-center gap-2 text-sm">
          <FlaskConical className="w-4 h-4 text-emerald-500" />
          <span className="font-semibold">{nDrugs}</span> co-folded drugs ·
          <span className="font-semibold">{complexes.length}</span> complexes ·
          <span className="font-semibold">{rows.length}</span> isoforms scored
        </div>
        <Link href={`/drugs/${gene.symbol}`}
          className="inline-flex items-center gap-1.5 text-sm font-semibold text-emerald-400 hover:text-emerald-300 hover:underline">
          <Pill className="w-4 h-4" /> Open Drug–Isoform viewer <ArrowRight className="w-4 h-4" />
        </Link>
      </div>

      {/* Pocket-fate mix across all drug × isoform pairs */}
      {total > 0 && (
        <div>
          <div className="flex justify-between text-xs text-muted-foreground mb-1.5">
            <span>Drug-pocket fate across isoforms</span>
            <span>{mix.lost + mix.disrupted} subtractive · {mix.preserved} preserved{mix.no_call ? ` · ${mix.no_call} no-call` : ""}</span>
          </div>
          <div className="flex h-2.5 rounded-full overflow-hidden">
            {(["lost", "disrupted", "preserved", "no_call"] as const).map(k =>
              mix[k] > 0 ? <div key={k} className={FATE[k].bar} style={{ width: `${100 * mix[k] / total}%` }} title={`${FATE[k].label}: ${mix[k]}`} /> : null
            )}
          </div>
        </div>
      )}

      {/* Per-isoform rollup, worst-affected first */}
      <div className="space-y-1">
        {rows.slice(0, 8).map(r => (
          <Link key={r.id} href={`/drugs/${gene.symbol}`}
            className="grid grid-cols-[1fr_auto_auto] gap-3 items-center px-3 py-1.5 rounded-md hover:bg-muted/50 text-sm">
            <span className="flex items-center gap-2 truncate">
              <span className={cn("w-2 h-2 rounded-full shrink-0", FATE[r.fate].bar)} />
              <span className="font-mono truncate">{r.transcript_id}</span>
              {r.mane && <Badge variant="secondary" className="text-[10px]">MANE</Badge>}
              {r.len && <span className="text-xs text-muted-foreground">{r.len} aa</span>}
            </span>
            <span className="text-xs text-muted-foreground tabular-nums">
              {r.meanFrac == null ? "—" : `${Math.round(r.meanFrac * 100)}% retained`}
            </span>
            <Badge className={cn("border", FATE[r.fate].cls)}>{FATE[r.fate].label}</Badge>
          </Link>
        ))}
        {rows.length > 8 && (
          <Link href={`/drugs/${gene.symbol}`} className="block px-3 py-1.5 text-xs text-muted-foreground hover:text-foreground">
            + {rows.length - 8} more isoforms →
          </Link>
        )}
      </div>
    </div>
  );
}
