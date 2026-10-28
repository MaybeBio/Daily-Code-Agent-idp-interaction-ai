"use client";

import { useMemo } from "react";
import Link from "next/link";
import { Pill, ArrowRight, FlaskConical } from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { cn } from "@/lib/utils";

interface Fate {
  id: number;
  drug_name?: string;
  pubchem_cid?: string;
  n_pocket_residues?: number;
  n_retained?: number;
  frac_retained?: number;
  fate?: string;
  mechanism_class?: string;
  decidable?: boolean;
  pocket_plddt?: number;
  retained_residues_json?: string;
}

const FATE = {
  lost:      { label: "Lost",      cls: "bg-red-500/15 text-red-400 border-red-500/30",       bar: "bg-red-500" },
  disrupted: { label: "Disrupted", cls: "bg-amber-500/15 text-amber-400 border-amber-500/30", bar: "bg-amber-500" },
  preserved: { label: "Preserved", cls: "bg-emerald-500/15 text-emerald-400 border-emerald-500/30", bar: "bg-emerald-500" },
  no_call:   { label: "No call",   cls: "bg-slate-500/15 text-slate-400 border-slate-500/30",  bar: "bg-slate-500" },
} as const;
type FateKey = keyof typeof FATE;
const fateOf = (f?: string): FateKey => (f && f in FATE ? (f as FateKey) : "no_call");

export default function IsoformDrugPanel({
  fates, geneSymbol, onSelect,
}: {
  fates: Fate[];
  geneSymbol?: string;
  onSelect?: (residues: number[]) => void;
}) {
  const mix = useMemo(() => {
    const m: Record<string, number> = { lost: 0, disrupted: 0, preserved: 0, no_call: 0 };
    for (const f of fates) m[fateOf(f.fate)]++;
    return m;
  }, [fates]);
  const total = fates.length;

  // worst-retention first (undecidable / no-call sink to the bottom)
  const rows = useMemo(
    () => [...fates].sort((a, b) => (a.frac_retained ?? 1.1) - (b.frac_retained ?? 1.1)),
    [fates]);

  const highlight = (f: Fate) => {
    if (!onSelect) return;
    try { onSelect(JSON.parse(f.retained_residues_json ?? "[]")); } catch { onSelect([]); }
  };

  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center justify-between gap-3">
        <div className="flex items-center gap-2 text-sm">
          <FlaskConical className="w-4 h-4 text-emerald-500" />
          Fate of <span className="font-semibold">{total}</span> drug pocket{total === 1 ? "" : "s"} in this isoform
        </div>
        {geneSymbol && (
          <Link href={`/drugs/${geneSymbol}`}
            className="inline-flex items-center gap-1.5 text-sm font-semibold text-emerald-400 hover:text-emerald-300 hover:underline">
            <Pill className="w-4 h-4" /> Open Drug–Isoform viewer <ArrowRight className="w-4 h-4" />
          </Link>
        )}
      </div>

      {total > 0 && (
        <div>
          <div className="flex justify-between text-xs text-muted-foreground mb-1.5">
            <span>Pocket fate across drugs</span>
            <span>{mix.lost + mix.disrupted} subtractive · {mix.preserved} preserved{mix.no_call ? ` · ${mix.no_call} no-call` : ""}</span>
          </div>
          <div className="flex h-2.5 rounded-full overflow-hidden">
            {(["lost", "disrupted", "preserved", "no_call"] as const).map(k =>
              mix[k] > 0 ? <div key={k} className={FATE[k].bar} style={{ width: `${100 * mix[k] / total}%` }} title={`${FATE[k].label}: ${mix[k]}`} /> : null
            )}
          </div>
        </div>
      )}

      <div className="space-y-1">
        {rows.map(f => {
          const st = FATE[fateOf(f.fate)];
          return (
            <button key={f.id} onClick={() => highlight(f)}
              className="w-full grid grid-cols-[1fr_auto_auto] gap-3 items-center px-3 py-1.5 rounded-md hover:bg-muted/50 text-sm text-left"
              title={onSelect ? "Highlight retained pocket residues in the viewer" : undefined}>
              <span className="flex items-center gap-2 truncate">
                <span className={cn("w-2 h-2 rounded-full shrink-0", st.bar)} />
                <span className="truncate font-medium">{f.drug_name ?? "drug"}</span>
                {f.mechanism_class && <span className="text-xs text-muted-foreground truncate">{f.mechanism_class.replace(/_/g, " ")}</span>}
              </span>
              <span className="text-xs text-muted-foreground tabular-nums whitespace-nowrap">
                {f.n_retained}/{f.n_pocket_residues} · {Math.round((f.frac_retained ?? 0) * 100)}%
              </span>
              <Badge className={cn("border", st.cls)}>{st.label}</Badge>
            </button>
          );
        })}
      </div>

      <p className="text-xs text-muted-foreground leading-relaxed">
        Each drug's binding pocket is defined from real contacts on the co-folded canonical structure, then carried
        onto this isoform by sequence alignment. Calls in low-confidence regions (pocket pLDDT &lt; 70) are reported as
        <span className="italic"> no&nbsp;call</span> rather than guessed. Click a drug to highlight its retained pocket residues.
      </p>
    </div>
  );
}
