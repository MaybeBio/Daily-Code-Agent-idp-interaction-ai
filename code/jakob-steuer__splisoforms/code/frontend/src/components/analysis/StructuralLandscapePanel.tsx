"use client";
import React, { useEffect, useState, useMemo } from "react";
import { Layers, Activity, Zap, AlertTriangle } from "lucide-react";
import { getApiUrl } from "@/lib/api";

interface ResidueImpact {
  residue_number: number;
  dssp_assignment: string;
  residue_plddt?: number;
}

interface DisorderScore {
  disorder_scores: number[];
  anchor_scores: number[];
}

interface StructuralMetrics {
  plddt_series?: number[];
  pae_blocks?: string;
}

const DSSP_COLORS: Record<string, string> = {
  H: "#f43f5e",
  E: "#3b82f6",
  C: "#475569",
};

/** AlphaFold 3 dark theme palette */
function plddt_v4_color(v: number | undefined): string {
  if (v === undefined || v === null) return "#475569";
  if (v >= 90) return "#0053D6";
  if (v >= 70) return "#63C7F4";
  if (v >= 50) return "#FFD600";
  return "#FF8A00";
}

function disorder_class(score: number): string {
  if (score < 0.4) return "bg-slate-200/50 dark:bg-slate-600/50";
  if (score < 0.7) return "bg-amber-500/60";
  return "bg-amber-500";
}

function InteractiveTrack({
  total,
  onHover,
  onSelect,
  hoveredResidues,
  selectedResidues,
  className,
  children
}: {
  total: number;
  onHover?: (residues: number[]) => void;
  onSelect?: (residues: number[]) => void;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  className?: string;
  children: React.ReactNode;
}) {
  const [dragStart, setDragStart] = useState<number | null>(null);

  useEffect(() => {
    const handleMouseUp = () => setDragStart(null);
    window.addEventListener("mouseup", handleMouseUp);
    return () => window.removeEventListener("mouseup", handleMouseUp);
  }, []);

  const getIdx = (e: React.MouseEvent) => {
    const rect = e.currentTarget.getBoundingClientRect();
    return Math.floor(((e.clientX - rect.left) / rect.width) * total) + 1;
  };

  return (
    <div 
      className={`flex h-4 rounded overflow-hidden border cursor-crosshair relative ${className}`}
      onMouseDown={(e) => {
        if (!total) return;
        const idx = getIdx(e);
        setDragStart(idx);
        if (onSelect) onSelect([idx]);
      }}
      onMouseMove={(e) => {
        if (!total) return;
        const idx = getIdx(e);
        if (dragStart !== null) {
          const start = Math.min(dragStart, idx);
          const end = Math.max(dragStart, idx);
          const range = Array.from({length: end - start + 1}, (_, i) => start + i);
          if (onSelect) onSelect(range);
        } else {
          if (onHover) onHover([idx]);
        }
      }}
      onMouseLeave={() => {
        if (dragStart === null && onHover) onHover([]);
      }}
    >
      {children}
      {hoveredResidues?.map(r => (
         <div key={`h-${r}`} className="absolute top-0 bottom-0 bg-primary/80 pointer-events-none mix-blend-screen" style={{ left: `${((r - 1) / total) * 100}%`, width: `${Math.max(0.1, 100 / total)}%`, minWidth: '1px' }} />
      ))}
      {selectedResidues?.map(r => (
         <div key={`s-${r}`} className="absolute top-0 bottom-0 bg-orange-500 pointer-events-none mix-blend-screen" style={{ left: `${((r - 1) / total) * 100}%`, width: `${Math.max(0.1, 100 / total)}%`, minWidth: '1px' }} />
      ))}
    </div>
  );
}

export default function StructuralLandscapePanel({ 
  isoformId,
  hoveredResidues,
  selectedResidues,
  onHover,
  onSelect
}: { 
  isoformId: number;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelect?: (residues: number[]) => void;
}) {
  const [residues, setResidues] = useState<ResidueImpact[]>([]);
  const [disorder, setDisorder] = useState<DisorderScore | null>(null);
  const [metrics, setMetrics] = useState<StructuralMetrics | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const apiUrl = getApiUrl();
    Promise.all([
      fetch(`${apiUrl}/api/isoform/${isoformId}/structural-impact`).then(r => r.ok ? r.json() : []),
      fetch(`${apiUrl}/api/isoform/${isoformId}`).then(r => r.ok ? r.json() : null)
    ])
    .then(([impactData, isoformData]) => {
      setResidues(impactData);
      if (isoformData) {
        if (isoformData.disorder_score) setDisorder(isoformData.disorder_score);
        if (isoformData.structural_metrics) setMetrics(isoformData.structural_metrics);
      }
    })
    .catch(() => {})
    .finally(() => setLoading(false));
  }, [isoformId]);

  const segments = useMemo(() => {
    const segs: Array<{ code: string; count: number }> = [];
    for (const r of residues) {
      const last = segs[segs.length - 1];
      if (last && last.code === r.dssp_assignment) last.count++;
      else segs.push({ code: r.dssp_assignment, count: 1 });
    }
    return segs;
  }, [residues]);

  const total = metrics?.plddt_series?.length || residues.length || (disorder?.disorder_scores.length ?? 0);

  // Overall pLDDT confidence — used to surface a banner when the whole
  // structure is below the conventional usable cutoff.
  const meanPlddt = useMemo(() => {
    const series = metrics?.plddt_series;
    if (!series || series.length === 0) return null;
    const vals = series.filter(v => typeof v === "number") as number[];
    if (vals.length === 0) return null;
    return vals.reduce((s, v) => s + v, 0) / vals.length;
  }, [metrics?.plddt_series]);

  if (loading) return <div className="p-4 text-center text-slate-500 animate-pulse text-[10px]">LOADING_LANDSCAPE...</div>;
  if (total === 0) return <div className="p-8 text-center text-slate-600 border border-slate-800 rounded-xl">DATA_PENDING</div>;

  return (
    <div className="space-y-8">
      {/* Low-overall-confidence banner (audit-6).  Conventional AlphaFold
          cutoff is pLDDT 70; below this most downstream structural-impact
          calls (DSSP transitions, SASA-based PTM accessibility, neoantigen
          structural confidence) are unreliable. */}
      {meanPlddt !== null && meanPlddt < 70 && (
        <div className="flex items-start gap-2 bg-amber-500/10 border border-amber-500/30 rounded-lg p-3 text-xs text-amber-300">
          <AlertTriangle className="w-4 h-4 shrink-0 mt-0.5" />
          <div>
            <div className="font-bold mb-0.5">
              Low-confidence prediction (mean pLDDT {meanPlddt.toFixed(1)})
            </div>
            <p className="leading-relaxed text-amber-200/85">
              Downstream calls that depend on side-chain placement
              (DSSP secondary structure, SASA-based PTM accessibility,
              structural rigidity) are systematically less reliable for
              this isoform.  Audits 2/4/6 routed many such calls to
              <code className="font-mono mx-1">low_confidence</code> /
              <code className="font-mono mx-1">undetermined</code> states.
            </p>
          </div>
        </div>
      )}

      {/* Confidence */}
      <div className="space-y-4 pb-4 border-b border-slate-800">
        <div className="flex justify-between text-[10px] font-bold text-slate-500 uppercase tracking-widest">
          <div className="flex items-center gap-2"><Zap className="w-3 h-3 text-blue-400" /> Confidence (pLDDT)</div>
          <div className="flex gap-2">
             <span className="text-[#0053D6]">90+</span>
             <span className="text-[#63C7F4]">70+</span>
             <span className="text-[#FFD600]">50+</span>
             <span className="text-[#FF8A00]">{"<"}50</span>
          </div>
        </div>
        <div className="relative">
          <InteractiveTrack 
            total={total}
            onHover={onHover}
            onSelect={onSelect}
            hoveredResidues={hoveredResidues}
            selectedResidues={selectedResidues}
            className="border-border bg-background"
          >
            {metrics?.plddt_series ? metrics.plddt_series.map((s, i) => (
               <div key={i} style={{ width: (100 / total) + "%", backgroundColor: plddt_v4_color(s) }} />
            )) : residues.map((r, i) => (
               <div key={i} style={{ width: (100 / total) + "%", backgroundColor: plddt_v4_color(r.residue_plddt) }} />
            ))}
          </InteractiveTrack>
          
          {/* Rigid Blocks Sub-track */}
          {metrics?.pae_blocks && (
            <div className="mt-1 flex h-1.5 w-full bg-slate-200/50 dark:bg-slate-900/50 rounded-full overflow-hidden relative border border-slate-300/50 dark:border-slate-800/50">
              {JSON.parse(metrics.pae_blocks).map((block: any, i: number) => (
                <div 
                  key={i} 
                  className="absolute h-full bg-indigo-500/60"
                  style={{ 
                    left: `${((block.start - 1) / total) * 100}%`, 
                    width: `${((block.end - block.start + 1) / total) * 100}%` 
                  }}
                />
              ))}
            </div>
          )}
        </div>

        {metrics?.pae_blocks && (
          <div className="flex gap-3 p-3 bg-blue-500/5 border border-blue-500/20 rounded-lg">
            <AlertTriangle className="w-4 h-4 text-blue-400 shrink-0" />
            <div className="text-[10px] text-blue-400/80 font-bold uppercase tracking-tight">
              Rigid structural blocks detected. These regions are predicted to be stable and are prioritized for structural impact assessment.
            </div>
          </div>
        )}
      </div>

      {/* DSSP */}
      <div className="space-y-2">
        <div className="flex justify-between text-[10px] font-bold text-slate-500 uppercase tracking-widest">
          <div className="flex items-center gap-2"><Layers className="w-3 h-3 text-rose-500" /> Secondary Structure (DSSP)</div>
          <div className="flex gap-2">
            <span className="text-rose-500">Helix</span>
            <span className="text-blue-500">Sheet</span>
            <span className="text-slate-500">Loop</span>
          </div>
        </div>
        <InteractiveTrack 
          total={total}
          onHover={onHover}
          onSelect={onSelect}
          hoveredResidues={hoveredResidues}
          selectedResidues={selectedResidues}
          className="border-slate-300 dark:border-slate-700 bg-slate-100 dark:bg-slate-800"
        >
          {segments.map((seg, i) => (
            <div key={i} style={{ width: (seg.count / total * 100) + "%", backgroundColor: DSSP_COLORS[seg.code] || "#475569" }} className="border-0" />
          ))}
        </InteractiveTrack>
      </div>

      {/* Disorder */}
      <div className="space-y-2">
        <div className="flex justify-between text-[10px] font-bold text-slate-500 uppercase tracking-widest">
          <div className="flex items-center gap-2"><Activity className="w-3 h-3 text-amber-500" /> Disorder propensity (AIUPred)</div>
          <div className="flex gap-2">
            <span className="text-slate-500">Ordered</span>
            <span className="text-amber-500/60">Flexible</span>
            <span className="text-amber-500">Disordered</span>
          </div>
        </div>
        <InteractiveTrack 
          total={total}
          onHover={onHover}
          onSelect={onSelect}
          hoveredResidues={hoveredResidues}
          selectedResidues={selectedResidues}
          className="border-slate-300 dark:border-slate-700 bg-slate-100 dark:bg-slate-800"
        >
          {disorder ? disorder.disorder_scores.map((s, i) => (
            <div key={i} style={{ width: (100 / total) + "%" }} className={`border-0 ${disorder_class(s)}`} />
          )) : <div className="w-full bg-slate-100 dark:bg-slate-800" />}
        </InteractiveTrack>
      </div>

      {/* Binding */}
      <div className="space-y-2">
        <div className="flex justify-between text-[10px] font-bold text-slate-500 uppercase tracking-widest">
          <div className="flex items-center gap-2"><Zap className="w-3 h-3 text-emerald-500" /> Functional Binding (AIUPred-binding)</div>
          <div className="flex gap-2">
            <span className="text-slate-500">Low Propensity</span>
            <span className="text-emerald-500">Binding Site</span>
          </div>
        </div>
        <InteractiveTrack 
          total={total}
          onHover={onHover}
          onSelect={onSelect}
          hoveredResidues={hoveredResidues}
          selectedResidues={selectedResidues}
          className="border-slate-300 dark:border-slate-700 bg-slate-100 dark:bg-slate-800"
        >
          {disorder ? disorder.anchor_scores.map((s, i) => (
            <div key={i} style={{ width: (100 / total) + "%" }} className={`border-0 ${s > 0.5 ? "bg-emerald-500/60" : "bg-transparent"}`} />
          )) : <div className="w-full bg-slate-100 dark:bg-slate-800" />}
        </InteractiveTrack>
      </div>
    </div>
  );
}
