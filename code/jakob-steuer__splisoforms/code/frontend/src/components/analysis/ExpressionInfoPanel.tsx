"use client";

import { useEffect, useState } from "react";
import { Microscope, BarChart2, Info } from "lucide-react";
import { getApiUrl } from "@/lib/api";
import { sourceShortLabel } from "@/lib/sources";

// ── HPA gene-level cell-type data ─────────────────────────────────────────────
export interface CellTypeEntry {
  id: number;
  gene_id: number;
  cell_type: string;
  expression_ncpm: number;
  source: string;
  specificity_category?: string;
}

// ── MAS-seq isoform-level per-sample data ─────────────────────────────────────
interface ExpressionRecord {
  id: number;
  isoform_id: number;
  sample: string;
  condition: string;
  fl_count: number;
  cpm?: number | null;
  isoform_frac?: number | null;
  n_cells?: number | null;
  cell_frac?: number | null;
}

// ── Shared helpers ─────────────────────────────────────────────────────────────
function HBar({ value, max, color }: { value: number; max: number; color: string }) {
  const pct = max > 0 ? Math.min((value / max) * 100, 100) : 0;
  return (
    <div className="flex items-center gap-2 min-w-0">
      <div className="flex-1 h-2 rounded-full bg-muted/40 overflow-hidden">
        <div className={`h-full rounded-full ${color}`} style={{ width: `${pct}%` }} />
      </div>
      <span className="text-[10px] font-mono tabular-nums w-14 text-right shrink-0">
        {value >= 1000 ? `${(value / 1000).toFixed(1)}k` : value.toFixed(value < 10 ? 2 : 1)}
      </span>
    </div>
  );
}

function FracBar({ value, color = "bg-violet-500" }: { value: number; color?: string }) {
  return (
    <div className="flex items-center gap-2 min-w-0">
      <div className="flex-1 h-2 rounded-full bg-muted/40 overflow-hidden">
        <div className={`h-full rounded-full ${color}`} style={{ width: `${Math.min(value * 100, 100)}%` }} />
      </div>
      <span className="text-[10px] font-mono tabular-nums w-10 text-right shrink-0">
        {(value * 100).toFixed(1)}%
      </span>
    </div>
  );
}

// Tissue-flag display used when no quantitative expression data is available
// (e.g. breast cancer isoforms ingested without per-sample FL counts)
interface TissueFlags {
  in_normal?: boolean | null;
  in_primary_tumor?: boolean | null;
  in_metastasis?: boolean | null;
}

const TISSUE_DEFS = [
  { key: "in_normal"        as const, label: "Normal",        color: "bg-sky-400",  text: "text-sky-300"  },
  { key: "in_primary_tumor" as const, label: "Primary Tumor", color: "bg-rose-500", text: "text-rose-300" },
  { key: "in_metastasis"    as const, label: "Metastasis",    color: "bg-violet-500", text: "text-violet-300" },
];

function TissuePresenceSection({ flags }: { flags: TissueFlags }) {
  const detected = TISSUE_DEFS.filter(t => flags[t.key]);
  const absent   = TISSUE_DEFS.filter(t => !flags[t.key]);

  return (
    <div className="space-y-2">
      <div className="flex flex-wrap gap-2 text-[10px]">
        {TISSUE_DEFS.map(t => (
          <span
            key={t.key}
            className={`flex items-center gap-1.5 px-2 py-0.5 rounded border ${
              flags[t.key]
                ? `${t.color.replace("bg-", "bg-").replace("-500", "-500/10")} ${t.text} border-current/30`
                : "bg-muted/20 text-muted-foreground/40 border-border/30"
            }`}
          >
            <span className={`w-1.5 h-1.5 rounded-full shrink-0 ${flags[t.key] ? t.color : "bg-muted"}`} />
            {t.label}
            {flags[t.key] ? " ✓" : " —"}
          </span>
        ))}
      </div>
      <p className="text-[10px] text-muted-foreground">
        Detected in {detected.length}/{TISSUE_DEFS.length} tissue compartment{detected.length !== 1 ? "s" : ""}.
        {absent.length > 0 && ` Absent in: ${absent.map(t => t.label).join(", ")}.`}
        {" "}Quantitative per-sample FL/CPM data not available for this dataset.
      </p>
    </div>
  );
}

// ── Sub-section: isoform-level expression ─────────────────────────────────────
function IsoformExpressionSection({
  expressionData,
  tissueFlags,
}: {
  expressionData: ExpressionRecord[];
  tissueFlags?: TissueFlags;
}) {
  // No quantitative data — fall back to tissue presence flags if available
  if (expressionData.length === 0) {
    const hasTissueFlags = tissueFlags && TISSUE_DEFS.some(t => tissueFlags[t.key]);
    if (hasTissueFlags && tissueFlags) {
      return <TissuePresenceSection flags={tissueFlags} />;
    }
    return (
      <div className="py-4 text-center text-muted-foreground">
        <BarChart2 className="w-6 h-6 mx-auto mb-1.5 opacity-40" />
        <p className="text-xs">Not quantified in this dataset.</p>
        <p className="text-[10px] mt-0.5 opacity-60">
          Isoform was not retained through the per-sample quality filter.
        </p>
      </div>
    );
  }

  // Derive sample order dynamically: normals first, then tumors, alphabetical within group
  const sampleOrder = [...new Set(expressionData.map(r => r.sample))].sort((a, b) => {
    const aIsNormal = expressionData.find(r => r.sample === a)?.condition !== "tumor";
    const bIsNormal = expressionData.find(r => r.sample === b)?.condition !== "tumor";
    if (aIsNormal !== bIsNormal) return aIsNormal ? -1 : 1;
    return a.localeCompare(b);
  });

  const get = (s: string) => expressionData.find(r => r.sample === s);
  const hasCells   = expressionData.some(r => r.n_cells != null);
  const maxCpm     = Math.max(...expressionData.map(r => r.cpm ?? 0));
  const tumorRecs  = expressionData.filter(r => r.condition === "tumor");
  const normalRecs = expressionData.filter(r => r.condition !== "tumor");
  const nSamples   = expressionData.length;

  const meanTumorCpm  = tumorRecs.length
    ? tumorRecs.reduce((s, r) => s + (r.cpm ?? 0), 0) / tumorRecs.length : null;
  const meanTumorFrac = tumorRecs.some(r => r.isoform_frac != null)
    ? tumorRecs.reduce((s, r) => s + (r.isoform_frac ?? 0), 0) / tumorRecs.length : null;
  const normalCpm  = normalRecs.length
    ? normalRecs.reduce((s, r) => s + (r.cpm ?? 0), 0) / normalRecs.length : null;

  return (
    <div className="space-y-3">
      {/* Summary chips */}
      <div className="flex flex-wrap gap-1.5 text-[10px]">
        <span className="bg-blue-500/10 text-blue-400 border border-blue-500/20 px-2 py-0.5 rounded">
          {nSamples} sample{nSamples !== 1 ? "s" : ""}
        </span>
        {meanTumorCpm != null && (
          <span className="bg-rose-500/10 text-rose-400 border border-rose-500/20 px-2 py-0.5 rounded">
            {meanTumorCpm.toFixed(1)} CPM (tumor mean)
          </span>
        )}
        {normalCpm != null && (
          <span className="bg-slate-500/10 text-foreground/70 border border-slate-500/20 px-2 py-0.5 rounded">
            {normalCpm.toFixed(1)} CPM (normal mean)
          </span>
        )}
        {meanTumorFrac != null && (
          <span className="bg-violet-500/10 text-violet-400 border border-violet-500/20 px-2 py-0.5 rounded">
            ψ {(meanTumorFrac * 100).toFixed(1)}% (tumor mean)
          </span>
        )}
      </div>

      {/* Per-sample table */}
      <table className="w-full text-xs">
        <thead>
          <tr className="border-b border-border text-muted-foreground text-left">
            <th className="py-1.5 px-1.5 w-20">Sample</th>
            <th className="py-1.5 px-1.5 w-16 text-right">FL</th>
            <th className="py-1.5 px-1.5">CPM</th>
            <th className="py-1.5 px-1.5">ψ fraction</th>
            {hasCells && <th className="py-1.5 px-1.5">Cells</th>}
          </tr>
        </thead>
        <tbody>
          {sampleOrder.map(sample => {
            const rec = get(sample);
            const isTumor = rec?.condition === "tumor";
            return (
              <tr key={sample} className={`border-b border-border/20 ${rec ? "hover:bg-muted/20" : "opacity-35"}`}>
                <td className="py-1.5 px-1.5">
                  <div className="flex items-center gap-1.5">
                    <span className={`w-1.5 h-1.5 rounded-full shrink-0 ${isTumor ? "bg-rose-500" : "bg-sky-400"}`} />
                    <span className={`font-mono font-medium text-[10px] ${isTumor ? "text-rose-300" : "text-sky-300"}`}>
                      {sample}
                    </span>
                  </div>
                </td>
                <td className="py-1.5 px-1.5 text-right font-mono text-foreground/70 text-[10px]">
                  {rec ? rec.fl_count.toLocaleString() : "—"}
                </td>
                <td className="py-1.5 px-1.5 min-w-[100px]">
                  {rec?.cpm != null
                    ? <HBar value={rec.cpm} max={maxCpm} color={isTumor ? "bg-rose-500" : "bg-sky-400"} />
                    : <span className="text-muted-foreground/30 text-[10px]">—</span>}
                </td>
                <td className="py-1.5 px-1.5 min-w-[100px]">
                  {rec?.isoform_frac != null
                    ? <FracBar value={rec.isoform_frac} />
                    : <span className="text-muted-foreground/30 text-[10px]">—</span>}
                </td>
                {hasCells && (
                  <td className="py-1.5 px-1.5 min-w-[90px]">
                    {rec?.cell_frac != null
                      ? <FracBar value={rec.cell_frac} color="bg-amber-500" />
                      : <span className="text-muted-foreground/30 text-[10px]">—</span>}
                  </td>
                )}
              </tr>
            );
          })}
        </tbody>
      </table>

      <p className="text-[10px] text-muted-foreground">
        Long-read sequencing · CPM normalised for sequencing depth · ψ = isoform share of gene FL reads
        {hasCells && " · Cell % from per-cell UMI sparse matrix"}
      </p>
    </div>
  );
}

// ── Sub-section: HPA gene-level cell-type expression ──────────────────────────
export function CellTypeSection({ data }: { data: CellTypeEntry[] }) {
  if (data.length === 0) {
    return (
      <div className="py-4 text-center text-muted-foreground">
        <Microscope className="w-6 h-6 mx-auto mb-1.5 opacity-40" />
        <p className="text-xs">No cell-type expression data available.</p>
      </div>
    );
  }

  const sorted   = [...data].sort((a, b) => b.expression_ncpm - a.expression_ncpm);
  const top      = sorted.slice(0, 20);
  const maxVal   = top[0]?.expression_ncpm || 1;
  const specificity = data[0]?.specificity_category;

  return (
    <div className="space-y-2">
      {specificity && (
        <div className="flex items-center gap-2 text-[10px]">
          <span className="text-muted-foreground">Specificity:</span>
          <span className="bg-primary/15 text-primary px-2 py-0.5 rounded capitalize">
            {specificity.replace(/_/g, " ")}
          </span>
          <span className="text-muted-foreground ml-auto">
            Source:{" "}
            <a href="https://www.proteinatlas.org/" target="_blank" rel="noreferrer"
              className="text-primary hover:underline">
              {data[0]?.source ?? "HPA"}
            </a>
          </span>
        </div>
      )}

      <div className="space-y-0.5">
        {top.map(entry => (
          <div key={entry.id} className="flex items-center gap-2 group">
            <span className="text-[10px] text-foreground/70 w-32 truncate text-right shrink-0
                             group-hover:text-foreground transition" title={entry.cell_type}>
              {entry.cell_type}
            </span>
            <div className="flex-1 h-3.5 bg-muted/50 rounded overflow-hidden">
              <div
                className="h-full rounded bg-gradient-to-r from-cyan-500/60 to-cyan-400/80 transition-all"
                style={{ width: `${Math.max((entry.expression_ncpm / maxVal) * 100, 1)}%` }}
              />
            </div>
            <span className="text-[10px] font-mono text-foreground/80 w-14 text-right shrink-0">
              {entry.expression_ncpm >= 1000
                ? `${(entry.expression_ncpm / 1000).toFixed(1)}k`
                : entry.expression_ncpm.toFixed(0)}
            </span>
          </div>
        ))}
      </div>

      <p className="text-[10px] text-muted-foreground text-center">
        Top {top.length}/{data.length} cell types · nCPM · Gene-level
      </p>
    </div>
  );
}

// ── Multi-source tissue presence (when isoform has memberships) ────────────────
interface SourceMembership {
  source_id: number;
  source_name: string;
  in_normal: boolean;
  in_primary_tumor: boolean;
  in_metastasis: boolean;
}

function MultiSourceTissueSection({ memberships }: { memberships: SourceMembership[] }) {
  // Every dataset the isoform belongs to is listed, so this matches the gene page's
  // "Mapped Transcripts". A dataset can legitimately carry NO tissue flag:
  //   * a reference annotation (GENCODE) has no tissue semantics at all — it is not a
  //     sampled cohort, so all three flags are false for all 75,895 of its memberships;
  //   * a long-read cohort can detect the transcript without assigning it a tissue.
  // Previously such rows were dropped and the whole section returned null, leaving the
  // heading above an empty gap even though the isoform clearly belonged to 2 datasets.
  const hasTissue = (m: SourceMembership) =>
    m.in_normal || m.in_primary_tumor || m.in_metastasis;

  return (
    <div className="space-y-2">
      <table className="w-full text-[10px]">
        <thead>
          <tr className="border-b border-border/40 text-muted-foreground text-left">
            <th className="py-1 pr-3">Dataset</th>
            {TISSUE_DEFS.map(t => (
              <th key={t.key} className="py-1 px-1.5 text-center">{t.label}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {memberships.map(m => (
            <tr key={m.source_id} className="border-b border-border/20">
              <td className="py-1.5 pr-3 text-foreground/70 max-w-[160px] truncate" title={m.source_name}>
                <a href="/documentation" className="hover:text-primary hover:underline">
                  {sourceShortLabel(m.source_name)}
                </a>
              </td>
              {hasTissue(m) ? (
                TISSUE_DEFS.map(t => (
                  <td key={t.key} className="py-1.5 px-1.5 text-center">
                    {m[t.key]
                      ? <span className={`inline-block w-3 h-3 rounded-full ${t.color}`} title="Detected" />
                      : <span className="inline-block w-3 h-3 rounded-full bg-muted/30" title="Not detected" />
                    }
                  </td>
                ))
              ) : (
                <td
                  colSpan={TISSUE_DEFS.length}
                  className="py-1.5 px-1.5 text-center text-muted-foreground/60 italic"
                  title="Present in this dataset, but it carries no tissue annotation for this isoform"
                >
                  present · no tissue annotation
                </td>
              )}
            </tr>
          ))}
        </tbody>
      </table>
      <p className="text-[10px] text-muted-foreground">
        Detection across {memberships.length} dataset{memberships.length !== 1 ? "s" : ""}.
        Reference annotations (GENCODE) list the transcript but carry no tissue sampling.
        Quantitative per-sample FL/CPM data shown below when available.
      </p>
    </div>
  );
}

// ── Combined panel ─────────────────────────────────────────────────────────────
interface CellTypePanelProps {
  geneId: number;
  expressionData?: ExpressionRecord[];
  tissueFlags?: TissueFlags;
  sourceName?: string;
  sourceMemberships?: SourceMembership[];
}

export default function CellTypePanel({ geneId, expressionData = [], tissueFlags, sourceName, sourceMemberships }: CellTypePanelProps) {
  const [cellTypeData, setCellTypeData] = useState<CellTypeEntry[]>([]);
  const [loading, setLoading]           = useState(true);

  useEffect(() => {
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/gene/${geneId}/cell-type`)
      .then(r => r.ok ? r.json() : [])
      .then(setCellTypeData)
      .catch(() => setCellTypeData([]))
      .finally(() => setLoading(false));
  }, [geneId]);

  if (loading) return (
    <div className="text-sm text-muted-foreground animate-pulse p-4">Loading expression data…</div>
  );

  return (
    <div className="space-y-6 p-1">

      {/* ── Section 1a: multi-dataset tissue detection (when memberships available) ── */}
      {sourceMemberships && sourceMemberships.length > 1 && (
        <div className="space-y-2">
          <div className="flex items-center gap-2 pb-1 border-b border-border/40">
            <BarChart2 className="w-3.5 h-3.5 text-cyan-400 shrink-0" />
            <h4 className="text-[11px] font-semibold uppercase tracking-widest text-foreground/80">
              Tissue Detection — All Datasets
            </h4>
          </div>
          <MultiSourceTissueSection memberships={sourceMemberships} />
        </div>
      )}

      {/* ── Section 1b: isoform-level quantitative expression / single-dataset tissue detection ── */}
      <div className="space-y-2">
        <div className="flex items-center gap-2 pb-1 border-b border-border/40">
          <BarChart2 className="w-3.5 h-3.5 text-rose-400 shrink-0" />
          <h4 className="text-[11px] font-semibold uppercase tracking-widest text-foreground/80">
            Isoform Expression
          </h4>
          <span className="ml-auto text-[10px] text-muted-foreground" title={sourceName || undefined}>
            {sourceName ? sourceShortLabel(sourceName) : "per sample"}
          </span>
        </div>
        <IsoformExpressionSection
          expressionData={expressionData}
          tissueFlags={sourceMemberships && sourceMemberships.length > 1 ? undefined : tissueFlags}
        />
      </div>

      {/* ── Section 2: HPA gene-level cell-type expression ── */}
      <div className="space-y-2">
        <div className="flex items-center gap-2 pb-1 border-b border-border/40">
          <Microscope className="w-3.5 h-3.5 text-cyan-400 shrink-0" />
          <h4 className="text-[11px] font-semibold uppercase tracking-widest text-foreground/80">
            Cell-Type Context — HPA (gene-level)
          </h4>
          <span className="ml-auto text-[10px] text-muted-foreground" title="HPA single-cell RNA nCPM for the cell types in which this gene is enriched/enhanced; not isoform-resolved">gene-level · scRNA nCPM</span>
        </div>
        <CellTypeSection data={cellTypeData} />
      </div>

    </div>
  );
}
