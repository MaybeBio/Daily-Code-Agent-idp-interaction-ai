"use client";

import React, { useMemo, useState, useEffect } from "react";
import { Grid3X3, Microscope } from "lucide-react";
import TissueVisualizer from "@/components/TissueVisualizer";
import { CellTypeSection, type CellTypeEntry } from "@/components/analysis/ExpressionInfoPanel";
import { getApiUrl } from "@/lib/api";

interface SourceMembership {
  source_id: number;
  source_name: string;
  in_normal: boolean;
  in_primary_tumor: boolean;
  in_metastasis: boolean;
}

interface Props {
  gene: any;
  selectedIsoformIds: number[];
  sortedIsoformIds?: number[];
}

const TISSUE_COLS = [
  { key: "in_normal"        as const, label: "Normal",        color: "rgba(34, 197, 94, "  },
  { key: "in_primary_tumor" as const, label: "Primary Tumor", color: "rgba(239, 68, 68, "  },
  { key: "in_metastasis"    as const, label: "Metastasis",    color: "rgba(168, 85, 247, " },
];

// Shorten long source names for compact display
function shortSource(name: string): string {
  if (name.includes("Karakulak") || name.toLowerCase().includes("ccrcc")) return "ccRCC";
  if (name.toLowerCase().includes("breast"))  return "BC";
  if (name.toLowerCase().includes("afdb") || name.toLowerCase().includes("alphafold")) return "AFDB";
  // Fallback: first word
  return name.split(/[\s(]/)[0];
}

export default function GeneExpressionHeatmap({ gene, selectedIsoformIds, sortedIsoformIds }: Props) {
  const displayIsoforms = useMemo(() => {
    if (selectedIsoformIds.length > 0) {
      const isoById = new Map((gene.isoforms ?? []).map((i: any) => [i.id, i]));
      if (sortedIsoformIds && sortedIsoformIds.length > 0) {
        return sortedIsoformIds
          .filter(id => selectedIsoformIds.includes(id))
          .map(id => isoById.get(id))
          .filter(Boolean);
      }
      return gene.isoforms?.filter((iso: any) => selectedIsoformIds.includes(iso.id)) || [];
    }
    return [...(gene.isoforms || [])].sort((a: any, b: any) => {
      if (a.is_canonical !== b.is_canonical) return a.is_canonical ? -1 : 1;
      if (a.is_flagged !== b.is_flagged) return a.is_flagged ? -1 : 1;
      return 0;
    }).slice(0, 5);
  }, [gene, selectedIsoformIds, sortedIsoformIds]);

  if (!displayIsoforms.length) return null;

  const allSourceNames = useMemo(() => {
    const names = new Set<string>();
    for (const iso of displayIsoforms) {
      const mbs: SourceMembership[] = iso.dataset_presence ?? [];
      mbs.forEach(m => names.add(m.source_name));
    }
    return Array.from(names);
  }, [displayIsoforms]);
  const multiSource = allSourceNames.length > 1;

  const [cellTypeData, setCellTypeData] = useState<CellTypeEntry[]>([]);
  
  useEffect(() => {
    if (gene?.id) {
      const apiUrl = getApiUrl();
      fetch(`${apiUrl}/api/gene/${gene.id}/cell-type`)
        .then(r => r.ok ? r.json() : [])
        .then(setCellTypeData)
        .catch(() => setCellTypeData([]));
    }
  }, [gene?.id]);

  const { cohortTissues, hpaTissues } = useMemo(() => {
    const cohort = new Set<string>();
    const hpa = new Set<string>();
    
    for (const iso of displayIsoforms) {
      if (iso.dataset_presence) {
        iso.dataset_presence.forEach((dp: any) => {
          const name = dp.source_name.toLowerCase();
          if (name.includes('ccrcc') || name.includes('kidney')) cohort.add('kidney');
          if (name.includes('brca') || name.includes('bc') || name.includes('breast')) cohort.add('breast');
        });
      }

      if (iso.expression_data) {
        iso.expression_data.forEach((expr: any) => {
          if (expr.sample) hpa.add(expr.sample);
        });
      }
    }

    // Also check gene-level expression_data just in case
    if (gene.expression_data) {
      gene.expression_data.forEach((expr: any) => {
        if (expr.sample) hpa.add(expr.sample);
      });
    }
    
    // Add all cell types from the HPA API endpoint
    cellTypeData.forEach(entry => {
      if (entry.cell_type) hpa.add(entry.cell_type);
    });

    return { cohortTissues: Array.from(cohort), hpaTissues: Array.from(hpa) };
  }, [displayIsoforms, gene, cellTypeData]);
  return (
    <div className="flex flex-col xl:flex-row gap-6">
      {/* Visualizer Panel */}
      <div className="flex-shrink-0 w-full xl:w-72 bg-card rounded-xl border border-border/50 p-4 flex flex-col items-center justify-start shadow-sm">
         <TissueVisualizer cohortTissues={cohortTissues} hpaTissues={hpaTissues} className="w-full flex justify-center mb-6" />
         
         <div className="w-full space-y-2 mt-auto">
           <div className="flex items-center gap-2 pb-1 border-b border-border/40">
             <Microscope className="w-3.5 h-3.5 text-cyan-400 shrink-0" />
             <h4 className="text-[11px] font-semibold uppercase tracking-widest text-foreground/80">
               Cell-Type Context — HPA (Canonical Isoform)
             </h4>
           </div>
           <CellTypeSection data={cellTypeData} />
         </div>
      </div>

      {/* Heatmap Panel */}
      <div className="flex-1 bg-card rounded-xl border border-border overflow-hidden p-4 shadow-sm">
        <div className="mb-4">
          <h3 className="font-semibold text-foreground flex items-center gap-2">
            <Grid3X3 className="w-4 h-4 text-primary" />
            Tumor/Normal Detection (Long-Read)
          </h3>
          <p className="text-xs text-muted-foreground mt-1">
            Indicates which sample types this transcript was detected in.
          </p>
        </div>

      <div className="overflow-x-auto pb-2">
        <table className="w-full text-left text-xs whitespace-nowrap">
          <thead>
            <tr>
              <th className="pb-3 px-2 font-medium text-muted-foreground">Transcript</th>
              <th className="pb-3 px-2 font-medium text-muted-foreground">Status</th>
              {multiSource && (
                <th className="pb-3 px-2 font-medium text-muted-foreground">Dataset</th>
              )}
              {TISSUE_COLS.map(s => (
                <th key={s.key} className="pb-3 px-2 font-medium text-muted-foreground text-center">
                  {s.label}
                </th>
              ))}
            </tr>
          </thead>
          <tbody className="divide-y divide-border/30">
            {displayIsoforms.map((iso: any) => {
              const mbs: SourceMembership[] = iso.dataset_presence ?? [];
              // If membership data available, render one sub-row per dataset;
              // otherwise fall back to the legacy boolean flags as a single row.
              const rows: Array<{ label: string | null; in_normal: boolean; in_primary_tumor: boolean; in_metastasis: boolean }> =
                mbs.length > 0
                  ? mbs.map(m => ({
                      label: m.source_name,
                      in_normal: m.in_normal,
                      in_primary_tumor: m.in_primary_tumor,
                      in_metastasis: m.in_metastasis,
                    }))
                  : [{
                      label: null,
                      in_normal: iso.in_normal,
                      in_primary_tumor: iso.in_primary_tumor,
                      in_metastasis: iso.in_metastasis,
                    }];

              return rows.map((row, rowIdx) => (
                <tr
                  key={`${iso.id}-${rowIdx}`}
                  className={`transition-colors ${rowIdx === 0 ? "hover:bg-muted/30" : "hover:bg-muted/20"} ${rowIdx > 0 ? "border-t-0" : ""}`}
                >
                  {/* Transcript + Status only on first sub-row, spanning all dataset rows */}
                  {rowIdx === 0 ? (
                    <>
                      <td className="py-2 px-2 font-mono text-primary/80" rowSpan={rows.length}>
                        {iso.transcript_name || iso.transcript_id}
                      </td>
                      <td className="py-2 px-2" rowSpan={rows.length}>
                        {iso.is_canonical
                          ? <span className="px-1.5 py-0.5 rounded-full bg-blue-500/20 text-blue-400 border border-blue-500/30 text-[10px]">Canonical</span>
                          : <span className="px-1.5 py-0.5 rounded-full bg-purple-500/20 text-purple-400 border border-purple-500/30 text-[10px]">Novel</span>
                        }
                      </td>
                    </>
                  ) : null}

                  {/* Dataset label column (only when multiple sources) */}
                  {multiSource && (
                    <td className="py-2 px-2 text-[10px] text-muted-foreground/70">
                      {row.label ? shortSource(row.label) : "—"}
                    </td>
                  )}

                  {TISSUE_COLS.map(s => {
                    const active = row[s.key];
                    return (
                      <td key={s.key} className="py-2 px-2 text-center w-24">
                        <div
                          className={`w-6 h-6 rounded mx-auto transition-all ${active ? "scale-110 shadow-lg" : "opacity-10 scale-90"}`}
                          style={{
                            backgroundColor: active ? `${s.color}0.8)` : "var(--muted)",
                            border: active ? `1px solid ${s.color}1.0)` : "1px solid var(--border)",
                          }}
                          title={`${row.label ?? ""} — ${active ? `Present in ${s.label}` : `Absent in ${s.label}`}`}
                        />
                      </td>
                    );
                  })}
                </tr>
              ));
            })}
          </tbody>
        </table>
      </div>

      <div className="mt-4 flex gap-4 text-[10px] text-muted-foreground uppercase tracking-widest font-bold">
        <div className="flex items-center gap-1"><div className="w-2 h-2 rounded bg-green-500" /> Normal</div>
        <div className="flex items-center gap-1"><div className="w-2 h-2 rounded bg-red-500" /> Primary Tumor</div>
        <div className="flex items-center gap-1"><div className="w-2 h-2 rounded bg-purple-500" /> Metastasis</div>
      </div>
      </div>
    </div>
  );
}
