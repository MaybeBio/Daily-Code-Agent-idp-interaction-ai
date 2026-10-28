"use client";
import { useEffect, useState, useMemo } from "react";
import { Layers, ChevronDown, ChevronUp, AlertTriangle, HelpCircle } from "lucide-react";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { getApiUrl } from "@/lib/api";

function InfoTip({ text }: { text: string }) {
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger>
          <HelpCircle className="w-2.5 h-2.5 text-muted-foreground/50 hover:text-muted-foreground transition-colors cursor-help ml-0.5" />
        </TooltipTrigger>
        <TooltipContent className="bg-popover border-border text-foreground text-xs max-w-xs shadow-xl font-normal normal-case tracking-normal">
          {text}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}

interface ResidueImpact {
  id: number;
  isoform_id: number;
  residue_number: number;
  dssp_assignment: string;
  phi_angle?: number;
  psi_angle?: number;
  residue_plddt?: number;
  exon_number?: number;
  canonical_dssp?: string;
  structural_change?: string;
  canonical_isoform_id?: number;
}

const DSSP_COLORS: Record<string, string> = {
  H: "#f43f5e",   // rose-500  — α-Helix
  E: "#3b82f6",   // blue-500  — β-Sheet
  C: "#475569",   // slate-600 — Coil
};

const DSSP_BG: Record<string, string> = {
  H: "bg-rose-500",
  E: "bg-blue-500",
  C: "bg-slate-600",
};

const DSSP_LABELS: Record<string, string> = {
  H: "α-Helix",
  E: "β-Sheet",
  C: "Coil / Unstructured",
};

/** pLDDT → tailwind colour for confidence indication */
function plddt_color(v: number | undefined): string {
  if (v === undefined || v === null) return "#475569";
  if (v >= 90) return "#22c55e";  // green-500
  if (v >= 70) return "#3b82f6";  // blue-500
  if (v >= 50) return "#f59e0b";  // amber-500
  return "#ef4444";               // red-500
}

/** Compress continuous runs of the same DSSP code into segments for efficient rendering */
function toSegments(residues: ResidueImpact[]): Array<{ code: string; count: number; avgPlddt?: number }> {
  if (!residues.length) return [];
  const segs: Array<{ code: string; count: number; plddt_sum: number; plddt_n: number }> = [];
  for (const r of residues) {
    const last = segs[segs.length - 1];
    if (last && last.code === r.dssp_assignment) {
      last.count++;
      if (r.residue_plddt != null) { last.plddt_sum += r.residue_plddt; last.plddt_n++; }
    } else {
      segs.push({ code: r.dssp_assignment, count: 1, plddt_sum: r.residue_plddt ?? 0, plddt_n: r.residue_plddt != null ? 1 : 0 });
    }
  }
  return segs.map(s => ({ code: s.code, count: s.count, avgPlddt: s.plddt_n > 0 ? s.plddt_sum / s.plddt_n : undefined }));
}

export default function StructuralImpactPanel({ isoformId }: { isoformId: number }) {
  const [residues, setResidues] = useState<ResidueImpact[]>([]);
  const [loading, setLoading] = useState(true);
  const [showTable, setShowTable] = useState(false);

  useEffect(() => {
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/isoform/${isoformId}/structural-impact`)
      .then(r => r.ok ? r.json() : [])
      .then(setResidues)
      .catch(() => setResidues([]))
      .finally(() => setLoading(false));
  }, [isoformId]);

  const segments = useMemo(() => toSegments(residues), [residues]);
  const allCoil = useMemo(() => residues.length > 0 && residues.every(r => r.dssp_assignment === "C"), [residues]);

  if (loading) return <div className="text-sm text-muted-foreground animate-pulse p-4">Loading structural impact…</div>;
  if (residues.length === 0) return (
    <div className="p-6 text-center text-muted-foreground">
      <Layers className="w-8 h-8 mx-auto mb-2 opacity-40" />
      <p className="text-sm">No structural impact data available.</p>
      <p className="text-xs mt-1 opacity-60 max-w-md mx-auto">
        Structural impact is measured per residue against the canonical protein, so it requires a
        predicted structure and a canonical to compare with. The canonical isoform itself has none.
      </p>
    </div>
  );

  // Stats
  const counts: Record<string, number> = { H: 0, E: 0, C: 0 };
  const changeCounts: Record<string, number> = {};
  residues.forEach(r => {
    counts[r.dssp_assignment] = (counts[r.dssp_assignment] || 0) + 1;
    if (r.structural_change && r.structural_change !== "same") {
      changeCounts[r.structural_change] = (changeCounts[r.structural_change] || 0) + 1;
    }
  });
  const total = residues.length;

  return (
    <div className="space-y-4 p-1">

      {/* Stale-data warning: if all assignments are C, likely from pre-fix backfill */}
      {allCoil && (
        <div className="flex items-start gap-2 bg-amber-500/8 border border-amber-500/25 rounded-lg p-3 text-xs text-amber-400">
          <AlertTriangle className="w-4 h-4 shrink-0 mt-0.5" />
          <div>
            <span className="font-semibold">All residues show as Coil</span> — this isoform was likely
            analyzed with an earlier version of the DSSP module that had a silent assignment bug.
            Re-run with <code className="bg-amber-500/10 px-1 rounded">--force</code> to get corrected secondary structure.
            Showing pLDDT confidence coloring as fallback.
          </div>
        </div>
      )}

      {/* Structure strip — segment-compressed for performance */}
      <div>
        <span className="text-xs text-muted-foreground mb-1 block">
          Secondary Structure Assignment
          <span className="ml-2 opacity-60">({total.toLocaleString()} residues, {segments.length} segments)</span>
        </span>
        <div
          className="flex h-5 rounded overflow-hidden border border-slate-700"
          title="Hover individual segments for details"
        >
          {allCoil
            /* When all-coil: show pLDDT confidence colors as a useful fallback */
            ? residues.map(r => (
                <div
                  key={r.residue_number}
                  className="flex-shrink-0"
                  style={{
                    width: `${100 / total}%`,
                    minWidth: "1px",
                    backgroundColor: plddt_color(r.residue_plddt),
                  }}
                  title={`Res ${r.residue_number}: pLDDT ${r.residue_plddt?.toFixed(0) ?? "N/A"}`}
                />
              ))
            /* Normal path: render compressed DSSP segments */
            : segments.map((seg, i) => (
                <div
                  key={i}
                  className="flex-shrink-0"
                  style={{
                    width: `${(seg.count / total) * 100}%`,
                    minWidth: "1px",
                    backgroundColor: DSSP_COLORS[seg.code] || DSSP_COLORS.C,
                  }}
                  title={`${DSSP_LABELS[seg.code] || "Coil"} × ${seg.count} residues${seg.avgPlddt != null ? `, avg pLDDT ${seg.avgPlddt.toFixed(0)}` : ""}`}
                />
              ))
          }
        </div>
      </div>

      {/* Legend */}
      <div className="flex gap-4 text-xs">
        {(["H", "E", "C"] as const).map(code => (
          <div key={code} className="flex items-center gap-1.5">
            <div className={`w-3 h-3 rounded-sm ${DSSP_BG[code]}`} />
            <span className="text-slate-400">{DSSP_LABELS[code]}</span>
            <span className="text-slate-300 font-medium">{((counts[code] / total) * 100).toFixed(0)}%</span>
          </div>
        ))}
        {allCoil && (
          <div className="flex items-center gap-1.5 ml-2 border-l border-slate-700 pl-2">
            <span className="text-slate-500 text-[10px]">↑ showing pLDDT instead:</span>
            {[["#22c55e","≥90"], ["#3b82f6","≥70"], ["#f59e0b","≥50"], ["#ef4444","<50"]].map(([c,l]) => (
              <span key={l} className="flex items-center gap-0.5">
                <span className="inline-block w-2.5 h-2.5 rounded-sm" style={{ backgroundColor: c }} />
                <span className="text-slate-500">{l}</span>
              </span>
            ))}
          </div>
        )}
      </div>

      {/* Structural changes */}
      {Object.keys(changeCounts).length > 0 && (
        <div className="bg-amber-500/5 border border-amber-500/20 rounded-lg p-3">
          <span className="text-xs text-amber-400 font-semibold block mb-1">Structural Changes vs Canonical</span>
          <div className="flex flex-wrap gap-2">
            {Object.entries(changeCounts).sort((a, b) => b[1] - a[1]).map(([change, count]) => (
              <span key={change} className="bg-amber-500/10 text-amber-300 px-2 py-0.5 rounded text-[10px]">
                {count}× {change.replace(/_/g, " → ")}
              </span>
            ))}
          </div>
        </div>
      )}

      {/* Table toggle */}
      <button onClick={() => setShowTable(!showTable)} className="flex items-center gap-1 text-xs text-muted-foreground hover:text-slate-300 transition">
        {showTable ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
        {showTable ? "Hide" : "Show"} per-residue data ({total.toLocaleString()} residues)
      </button>

      {showTable && (
        <div className="overflow-x-auto max-h-[320px] overflow-y-auto">
          <table className="w-full text-xs">
            <thead className="sticky top-0 bg-slate-900">
              <tr className="border-b border-slate-700 text-muted-foreground">
                <th className="text-left py-2 px-2">#</th>
                <th className="text-center py-2 px-2">
                  <span className="inline-flex items-center gap-0.5">
                    DSSP
                    <InfoTip text="Secondary structure (DSSP): H = helix, E = strand, C = coil." />
                  </span>
                </th>
                <th className="text-right py-2 px-2">
                  <span className="inline-flex items-center justify-end gap-0.5">
                    φ
                    <InfoTip text="Backbone torsion φ, around the N–Cα bond (−180° to 180°)." />
                  </span>
                </th>
                <th className="text-right py-2 px-2">
                  <span className="inline-flex items-center justify-end gap-0.5">
                    ψ
                    <InfoTip text="Backbone torsion ψ, around the Cα–C bond (−180° to 180°)." />
                  </span>
                </th>
                <th className="text-right py-2 px-2">
                  <span className="inline-flex items-center justify-end gap-0.5">
                    pLDDT
                    <InfoTip text="AlphaFold confidence (0–100): ≥90 very high, ≥70 confident, ≥50 low, <50 very low." />
                  </span>
                </th>
                <th className="text-center py-2 px-2">
                  <span className="inline-flex items-center gap-0.5">
                    Change
                    <InfoTip text="Structural change vs the canonical at this residue: gained, lost, or altered." />
                  </span>
                </th>
              </tr>
            </thead>
            <tbody>
              {residues.slice(0, 300).map(r => (
                <tr key={r.residue_number} className="border-b border-slate-800/50 hover:bg-slate-800/30">
                  <td className="py-1 px-2 font-mono text-slate-400">{r.residue_number}</td>
                  <td className="py-1 px-2 text-center">
                    <span
                      className={`inline-flex items-center justify-center w-5 h-5 rounded text-[10px] font-bold text-white ${DSSP_BG[r.dssp_assignment] || "bg-slate-600"}`}
                    >
                      {r.dssp_assignment}
                    </span>
                  </td>
                  <td className="py-1 px-2 text-right font-mono text-slate-400">{r.phi_angle?.toFixed(0) ?? "—"}</td>
                  <td className="py-1 px-2 text-right font-mono text-slate-400">{r.psi_angle?.toFixed(0) ?? "—"}</td>
                  <td className="py-1 px-2 text-right font-mono">
                    <span style={{ color: plddt_color(r.residue_plddt) }}>{r.residue_plddt?.toFixed(0) ?? "—"}</span>
                  </td>
                  <td className="py-1 px-2 text-center text-[10px] text-slate-500 capitalize">
                    {r.structural_change === "same" || !r.structural_change ? "—" : r.structural_change.replace(/_/g, " → ")}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          {total > 300 && (
            <p className="text-[10px] text-muted-foreground text-center py-2">
              Showing first 300 of {total.toLocaleString()} residues
            </p>
          )}
        </div>
      )}
    </div>
  );
}
