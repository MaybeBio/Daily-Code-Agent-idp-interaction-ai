"use client";

import { useEffect, useMemo, useState, Fragment } from "react";
import Link from "next/link";
import { useRouter } from "next/navigation";
import dynamic from "next/dynamic";
import {
  ArrowRight, ShieldAlert, Loader2, X, Plus, Layers, Info,
} from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { getApiUrl } from "@/lib/api";
import { sourceShortLabel } from "@/lib/sources";

// Mol* multi-isoform overlay (loads + superposes each isoform on the canonical).
const GeneStructureViewerClient = dynamic(
  () => import("@/app/gene/[id]/components/GeneStructureViewerClient"),
  { ssr: false, loading: () => (
    <div className="flex items-center justify-center h-full min-h-[400px] text-muted-foreground gap-2">
      <Loader2 className="w-6 h-6 animate-spin text-primary" /> Loading 3D overlay…
    </div>
  ) }
);

const MAX_COMPARE = 5;

interface IsoformDetail {
  id: number;
  gene_id: number;
  transcript_id: string;
  transcript_name?: string | null;
  sequence?: string;
  delta_length?: number;
  is_canonical?: boolean;
  in_primary_tumor?: boolean;
  in_normal?: boolean;
  in_metastasis?: boolean;
  structures?: { af3_cif_url?: string }[];
  gene?: { symbol?: string; canonical_sequence?: string };
  source?: { name?: string };
  metadata_record?: { structural_category?: string; subcategory?: string } | null;
  pfam_domains?: { domain_status?: string; pfam_name?: string }[];
  ted_domains?: any[];
  neoantigen_predictions?: { binding_affinity_nm?: number; ptm_impact_status?: string }[];
  ptm_accessibility?: { accessibility_change?: string }[];
  nmd_prediction?: { nmd_status?: string; nmd_confidence?: number; orf_type?: string } | null;
  structural_metrics?: { plddt_mean?: number; plddt_min?: number; pae_mean?: number } | null;
  string_interactions?: unknown[];
  homodimer_analysis?: { disruption_score?: number } | null;
  heterodimer_analyses?: { disruption_score?: number }[];
  disorder_score?: { disorder_scores?: number[] } | null;
}

const NMD_STYLE: Record<string, string> = {
  nmd_target: "bg-rose-500/15 text-rose-400 border-rose-500/30",
  nmd_escape: "bg-emerald-500/15 text-emerald-400 border-emerald-500/30",
};
function nmdBadge(status?: string) {
  if (!status) return <span className="text-muted-foreground">—</span>;
  const cls = NMD_STYLE[status] || "bg-muted text-muted-foreground border-border";
  return <Badge variant="outline" className={`text-[10px] ${cls}`}>{status.replace(/_/g, " ")}</Badge>;
}

export default function CompareViewClient({ ids }: { ids: number[] }) {
  const router = useRouter();
  const [isoforms, setIsoforms] = useState<IsoformDetail[]>([]);
  const [geneIsoforms, setGeneIsoforms] = useState<any[]>([]);
  const [canonical, setCanonical] = useState<IsoformDetail | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    const apiUrl = getApiUrl();
    (async () => {
      setLoading(true);
      setError(null);
      try {
        const details: IsoformDetail[] = await Promise.all(
          ids.map((id) => fetch(`${apiUrl}/api/isoform/${id}`).then((r) => {
            if (!r.ok) throw new Error(`Failed to load isoform ${id}`);
            return r.json();
          }))
        );
        if (cancelled) return;
        setIsoforms(details);

        // Same-gene: resolve the gene's full isoform list (for canonical ref + add picker)
        const geneId = details[0]?.gene_id;
        if (geneId != null) {
          const list = await fetch(`${apiUrl}/api/genes/${geneId}/isoforms`).then((r) => r.ok ? r.json() : []);
          if (cancelled) return;
          setGeneIsoforms(Array.isArray(list) ? list : []);

          // Canonical reference for the overlay: prefer one already in the set,
          // else fetch the gene's canonical isoform detail (for cif + sequence).
          let canon = details.find((d) => d.is_canonical) || null;
          if (!canon) {
            const canonItem = (Array.isArray(list) ? list : []).find((i: any) => i.is_canonical);
            if (canonItem) {
              canon = await fetch(`${apiUrl}/api/isoform/${canonItem.id}`).then((r) => r.ok ? r.json() : null);
            }
          }
          if (!cancelled) setCanonical(canon);
        }
      } catch (err: any) {
        if (!cancelled) setError(err.message || "Failed to load comparison data.");
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, [ids]);

  const crossGene = useMemo(
    () => isoforms.length > 1 && new Set(isoforms.map((i) => i.gene_id)).size > 1,
    [isoforms]
  );

  // ── Self-selection helpers ────────────────────────────────────────────────
  const setIds = (next: number[]) => router.push(`/compare?ids=${next.join(",")}`);
  const removeId = (id: number) => { if (ids.length > 2) setIds(ids.filter((x) => x !== id)); };
  const addId = (id: number) => { if (ids.length < MAX_COMPARE && !ids.includes(id)) setIds([...ids, id]); };
  const addable = geneIsoforms.filter((i: any) => !ids.includes(i.id));

  if (loading) {
    return (
      <div className="flex flex-col items-center justify-center min-h-[40vh] gap-3">
        <Loader2 className="w-6 h-6 text-primary animate-spin" />
        <span className="text-muted-foreground">Synchronising structural states…</span>
      </div>
    );
  }
  if (error || isoforms.length < 2) {
    return (
      <Card className="bg-destructive/10 border-destructive/20 text-destructive p-6">
        <div className="flex items-center gap-3">
          <ShieldAlert className="w-6 h-6" />
          <div>
            <h3 className="font-bold">Comparison Error</h3>
            <p className="text-sm opacity-80">{error || "Could not load isoforms for comparison."}</p>
          </div>
        </div>
      </Card>
    );
  }

  // ── Metric helpers ─────────────────────────────────────────────────────────
  const DASH = <span className="text-muted-foreground">—</span>;
  const ZERO = <span className="text-muted-foreground">0</span>;
  const pfamCount = (i: IsoformDetail, status: string) =>
    (i.pfam_domains || []).filter((d) => d.domain_status === status).length;
  const disorderPct = (i: IsoformDetail): number | null => {
    const s = i.disorder_score?.disorder_scores;
    if (!s || s.length === 0) return null;
    return Math.round((s.filter((x) => x > 0.5).length / s.length) * 100);
  };
  const goldMedal = (i: IsoformDetail) =>
    (i.neoantigen_predictions || []).filter((n) => n.ptm_impact_status === "Splice-Enabled (Gold Medal)").length;
  const strongBinders = (i: IsoformDetail) =>
    (i.neoantigen_predictions || []).filter((n) => (n.binding_affinity_nm ?? Infinity) < 50).length;
  const ptmLost = (i: IsoformDetail) =>
    (i.ptm_accessibility || []).filter((p) => p.accessibility_change === "lost").length;
  const interfaceDisruption = (i: IsoformDetail): number | null => {
    const scores = [
      i.homodimer_analysis?.disruption_score,
      ...(i.heterodimer_analyses || []).map((h) => h.disruption_score),
    ].filter((x): x is number => x != null);
    if (scores.length === 0) return null;
    return Math.round(Math.max(...scores) * 100);
  };
  const posInt = (n: number, cls = "") =>
    n > 0 ? <span className={cls}>{n}</span> : ZERO;
  const tissues = (i: IsoformDetail) => {
    const t = [i.in_normal && "N", i.in_primary_tumor && "T", i.in_metastasis && "M"].filter(Boolean);
    return t.length ? t.join(" · ") : DASH;
  };

  // Grouped so the ~20 metrics stay scannable. Each render pulls straight from the
  // isoform-detail payload (already fetched) — no extra requests.
  const metricGroups: { title: string; rows: { label: string; render: (i: IsoformDetail) => React.ReactNode }[] }[] = [
    { title: "Sequence & annotation", rows: [
      { label: "Length (aa)", render: (i) => i.sequence?.length ?? DASH },
      { label: "Δ vs canonical", render: (i) => {
          const d = i.delta_length;
          if (d == null || d === 0) return ZERO;
          return <span className={d > 0 ? "text-emerald-400" : "text-rose-400"}>{d > 0 ? "+" : ""}{d} aa</span>;
        } },
      { label: "Category", render: (i) => {
          const c = i.metadata_record?.structural_category;
          return c ? <span className="capitalize">{c.replace(/_/g, " ")}</span> : DASH;
        } },
      { label: "Source", render: (i) => <span className="text-xs">{sourceShortLabel(i.source?.name) || DASH}</span> },
    ] },
    { title: "Structure confidence", rows: [
      { label: "Mean pLDDT", render: (i) => i.structural_metrics?.plddt_mean != null ? i.structural_metrics.plddt_mean.toFixed(1) : DASH },
      { label: "Min pLDDT", render: (i) => i.structural_metrics?.plddt_min != null ? i.structural_metrics.plddt_min.toFixed(1) : DASH },
      { label: "Mean PAE (Å)", render: (i) => i.structural_metrics?.pae_mean != null ? i.structural_metrics.pae_mean.toFixed(1) : DASH },
    ] },
    { title: "Domains", rows: [
      { label: "Pfam present", render: (i) => posInt(pfamCount(i, "present")) },
      { label: "Pfam lost", render: (i) => posInt(pfamCount(i, "lost"), "text-rose-400") },
      { label: "Pfam gained", render: (i) => posInt(pfamCount(i, "gained"), "text-emerald-400") },
      { label: "Pfam truncated", render: (i) => posInt(pfamCount(i, "truncated"), "text-amber-400") },
      { label: "TED / structural domains", render: (i) => posInt((i.ted_domains || []).length) },
    ] },
    { title: "Transcript fate", rows: [
      { label: "NMD status", render: (i) => nmdBadge(i.nmd_prediction?.nmd_status) },
      { label: "ORF type", render: (i) => {
          const o = i.nmd_prediction?.orf_type;
          return o ? <span className="capitalize text-xs">{o.replace(/_/g, " ")}</span> : DASH;
        } },
      { label: "Disorder (% residues)", render: (i) => { const p = disorderPct(i); return p == null ? DASH : `${p}%`; } },
    ] },
    { title: "Neoantigens & PTM", rows: [
      { label: "Neoantigens", render: (i) => posInt((i.neoantigen_predictions || []).length) },
      { label: "Gold-Medal", render: (i) => posInt(goldMedal(i), "text-amber-400") },
      { label: "Strong binders (<50 nM)", render: (i) => posInt(strongBinders(i), "text-rose-400") },
      { label: "PTM sites", render: (i) => posInt((i.ptm_accessibility || []).length) },
      { label: "PTM sites lost", render: (i) => posInt(ptmLost(i), "text-orange-400") },
    ] },
    { title: "Interactions & expression", rows: [
      { label: "STRING partners", render: (i) => posInt((i.string_interactions || []).length) },
      { label: "Interface disruption", render: (i) => { const d = interfaceDisruption(i); return d == null ? DASH : <span className={d > 0 ? "text-rose-400" : ""}>{d}%</span>; } },
      { label: "Detected in (N/T/M)", render: (i) => tissues(i) },
    ] },
  ];

  const geneSymbol = isoforms[0]?.gene?.symbol;

  return (
    <div className="space-y-8 animate-in fade-in duration-500">
      {/* Selection chips + add */}
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-xs uppercase font-bold tracking-widest text-muted-foreground mr-1">
          Comparing {isoforms.length}{geneSymbol ? ` · ${geneSymbol}` : ""}
        </span>
        {isoforms.map((i) => (
          <span key={i.id} className="flex items-center gap-1.5 pl-3 pr-1.5 py-1 rounded-full border border-border/60 bg-card/50 text-sm">
            <span className="font-mono text-xs">{i.transcript_name || i.transcript_id}</span>
            {i.is_canonical && <Badge variant="outline" className="text-[9px] bg-amber-500/10 text-amber-500 border-amber-500/20">canon</Badge>}
            {ids.length > 2 && (
              <button onClick={() => removeId(i.id)} className="w-4 h-4 rounded-full hover:bg-muted flex items-center justify-center text-muted-foreground hover:text-foreground" title="Remove">
                <X className="w-3 h-3" />
              </button>
            )}
          </span>
        ))}
        {ids.length < MAX_COMPARE && addable.length > 0 && (
          <div className="relative inline-flex items-center">
            <Plus className="w-3.5 h-3.5 text-muted-foreground absolute left-2.5 pointer-events-none" />
            <select
              value=""
              onChange={(e) => { const v = parseInt(e.target.value, 10); if (!isNaN(v)) addId(v); }}
              className="appearance-none pl-7 pr-3 py-1 rounded-full border border-dashed border-border/60 bg-transparent text-sm text-muted-foreground hover:text-foreground hover:border-primary/40 cursor-pointer focus:outline-none"
            >
              <option value="">Add isoform…</option>
              {addable.map((i: any) => (
                <option key={i.id} value={i.id}>{i.transcript_id}{i.is_canonical ? " (canonical)" : ""}</option>
              ))}
            </select>
          </div>
        )}
      </div>

      {crossGene && (
        <div className="rounded-xl bg-amber-500/5 border border-amber-500/20 p-4 text-sm flex gap-3 items-start">
          <Info className="w-4 h-4 text-amber-400 shrink-0 mt-0.5" />
          <p className="text-foreground/85 leading-relaxed">
            These isoforms span different genes. Canonical-relative metrics (Δlength, domain loss) and the
            3D overlay assume a shared canonical reference and may not be meaningful across genes.
          </p>
        </div>
      )}

      {/* Combined 3D overlay — GeneStructureViewerClient brings its own bordered container */}
      <div className="space-y-2">
        <div className="flex items-center gap-2 text-xs uppercase tracking-widest text-primary font-bold">
          <Layers className="w-3.5 h-3.5" /> Combined Structural Overlay
        </div>
        {canonical?.structures?.[0]?.af3_cif_url ? (
          <GeneStructureViewerClient
            canonicalCifUrl={canonical.structures[0].af3_cif_url}
            canonicalSequence={canonical.sequence}
            selectedIsoforms={isoforms.map((i) => ({
              id: i.id,
              transcript_id: i.transcript_id,
              sequence: i.sequence,
              cifUrl: i.structures?.[0]?.af3_cif_url,
            }))}
          />
        ) : (
          <div className="flex items-center justify-center h-[400px] rounded-xl border border-border/40 bg-card/40 text-muted-foreground text-sm">
            No canonical structure available for overlay.
          </div>
        )}
      </div>

      {/* Metrics comparison table */}
      <Card className="bg-card/40 border-border/40">
        <CardHeader className="py-3 border-b border-border/20">
          <CardTitle className="text-xs uppercase tracking-widest text-muted-foreground">Metric Comparison</CardTitle>
        </CardHeader>
        <CardContent className="p-0 overflow-x-auto">
          <table className="w-full text-sm border-collapse">
            <thead>
              <tr className="border-b border-border/40">
                <th className="text-left px-4 py-3 font-semibold text-muted-foreground sticky left-0 bg-card/40">Metric</th>
                {isoforms.map((i) => (
                  <th key={i.id} className="text-left px-4 py-3 min-w-[150px]">
                    <Link href={`/isoform/${i.id}`} className="font-mono text-xs text-primary hover:underline flex items-center gap-1">
                      {i.transcript_name || i.transcript_id} <ArrowRight className="w-3 h-3" />
                    </Link>
                    {i.is_canonical && <Badge variant="outline" className="mt-1 text-[9px] bg-amber-500/10 text-amber-500 border-amber-500/20">canonical</Badge>}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {metricGroups.map((group) => (
                <Fragment key={group.title}>
                  <tr className="bg-muted/30 border-t border-border/30">
                    <td colSpan={isoforms.length + 1}
                        className="px-4 py-1.5 text-[10px] font-bold uppercase tracking-widest text-muted-foreground/80 sticky left-0 bg-muted/30">
                      {group.title}
                    </td>
                  </tr>
                  {group.rows.map((row, ri) => (
                    <tr key={row.label} className={ri % 2 ? "bg-muted/10" : ""}>
                      <td className="px-4 py-2.5 font-medium text-muted-foreground sticky left-0 bg-card/40">{row.label}</td>
                      {isoforms.map((i) => (
                        <td key={i.id} className="px-4 py-2.5 text-foreground">{row.render(i)}</td>
                      ))}
                    </tr>
                  ))}
                </Fragment>
              ))}
            </tbody>
          </table>
        </CardContent>
      </Card>
    </div>
  );
}
