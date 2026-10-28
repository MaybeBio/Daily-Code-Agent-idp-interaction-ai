"use client";

import { LayoutGrid } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { cn } from "@/lib/utils";

interface FateLike {
  isoform_id: number;
  fate?: string;
  decidable?: boolean;
}

function Stat({ value, label }: { value: number; label: string }) {
  return (
    <div className="flex flex-col items-center justify-center rounded-lg border border-border/40 bg-muted/20 py-2.5 px-1">
      <span className="text-xl font-bold text-foreground leading-none">{value.toLocaleString()}</span>
      <span className="text-[10px] uppercase tracking-widest text-muted-foreground mt-1 text-center leading-tight">
        {label}
      </span>
    </div>
  );
}

/**
 * Pooled counts for the gene's drug layer, plus a one-glance read on what splicing does to
 * its druggable pockets.
 *
 * The impact bar is computed over DECIDABLE pocket calls only. Undecidable ones (`no_call`
 * — the pocket sits in a low-confidence region, so retention cannot be judged) are reported
 * separately rather than folded in: counting them as "preserved" would overstate
 * druggability, and as "lost" would overstate the damage.
 */
export default function DrugGeneOverviewPanel({
  nComplexes,
  nDrugs,
  nIsoforms,
  nPockets,
  fates,
}: {
  nComplexes: number;
  nDrugs: number;
  nIsoforms: number;
  nPockets: number;
  fates: FateLike[];
}) {
  const lost = fates.filter(f => f.fate === "lost").length;
  const disrupted = fates.filter(f => f.fate === "disrupted").length;
  const preserved = fates.filter(f => f.fate === "preserved").length;
  const noCall = fates.filter(f => (f.fate ?? "no_call") === "no_call").length;

  const decidable = lost + disrupted + preserved;
  const subtractive = lost + disrupted;
  const pct = (n: number) => (decidable ? (n / decidable) * 100 : 0);

  // Isoforms that lose or disrupt at least one pocket — the headline biological claim.
  const hitIsoforms = new Set(
    fates.filter(f => f.fate === "lost" || f.fate === "disrupted").map(f => f.isoform_id),
  ).size;
  const judgedIsoforms = new Set(
    fates.filter(f => ["lost", "disrupted", "preserved"].includes(f.fate ?? "")).map(f => f.isoform_id),
  ).size;

  return (
    <Card>
      <CardHeader className="py-4">
        <CardTitle className="text-lg font-semibold flex items-center gap-2">
          <LayoutGrid className="w-4 h-4 text-emerald-500" />
          Overview
        </CardTitle>
      </CardHeader>

      <CardContent className="pt-0 space-y-4">
        <div className="grid grid-cols-4 gap-2">
          <Stat value={nComplexes} label="Co-folded" />
          <Stat value={nDrugs} label="Drugs" />
          <Stat value={nPockets} label="Pockets" />
          <Stat value={nIsoforms} label="Isoforms" />
        </div>

        <div className="space-y-2">
          <div className="flex items-baseline justify-between gap-2">
            <span className="text-[10px] uppercase tracking-widest text-muted-foreground font-semibold">
              Splicing impact on drug binding
            </span>
            {decidable > 0 && (
              <span className="text-[10px] text-muted-foreground">
                {decidable.toLocaleString()} scored
              </span>
            )}
          </div>

          {decidable === 0 ? (
            <p className="text-xs text-muted-foreground leading-relaxed">
              No pocket could be scored for this gene — every isoform–drug pair falls in a
              low-confidence region of the structure.
            </p>
          ) : (
            <>
              <div className="flex h-2 w-full overflow-hidden rounded-full bg-muted">
                {preserved > 0 && (
                  <div className="bg-emerald-500" style={{ width: `${pct(preserved)}%` }}
                       title={`Preserved: ${preserved}`} />
                )}
                {disrupted > 0 && (
                  <div className="bg-amber-500" style={{ width: `${pct(disrupted)}%` }}
                       title={`Disrupted: ${disrupted}`} />
                )}
                {lost > 0 && (
                  <div className="bg-red-500" style={{ width: `${pct(lost)}%` }}
                       title={`Lost: ${lost}`} />
                )}
              </div>

              <div className="flex flex-wrap gap-x-4 gap-y-1 text-[11px]">
                <Legend cls="bg-emerald-500" label="Preserved" n={preserved} pct={pct(preserved)} />
                <Legend cls="bg-amber-500" label="Disrupted" n={disrupted} pct={pct(disrupted)} />
                <Legend cls="bg-red-500" label="Lost" n={lost} pct={pct(lost)} />
              </div>

              <p className="text-xs text-muted-foreground leading-relaxed">
                <span className="text-foreground font-medium">
                  {hitIsoforms} of {judgedIsoforms} isoform{judgedIsoforms === 1 ? "" : "s"}
                </span>{" "}
                lose or disrupt at least one drug pocket, and{" "}
                <span className="text-foreground font-medium">
                  {Math.round(pct(subtractive))}%
                </span>{" "}
                of scored isoform–drug pairs are subtractive.
                {noCall > 0 && (
                  <>
                    {" "}
                    A further {noCall.toLocaleString()} pair{noCall === 1 ? " was" : "s were"} left
                    unscored: the pocket lies in a low-confidence region, so retention cannot be
                    judged.
                  </>
                )}
              </p>
            </>
          )}
        </div>
      </CardContent>
    </Card>
  );
}

function Legend({ cls, label, n, pct }: { cls: string; label: string; n: number; pct: number }) {
  return (
    <span className="inline-flex items-center gap-1.5 text-muted-foreground">
      <span className={cn("w-2 h-2 rounded-full shrink-0", cls)} />
      {label}
      <span className="text-foreground/80 font-medium">{n}</span>
      <span className="text-muted-foreground/70">({Math.round(pct)}%)</span>
    </span>
  );
}
