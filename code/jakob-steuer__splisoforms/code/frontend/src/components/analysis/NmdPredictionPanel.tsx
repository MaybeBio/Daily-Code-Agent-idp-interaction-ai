"use client";
import { useEffect, useState } from "react";
import { getApiUrl } from "@/lib/api";
import {
  AlertTriangle,
  CheckCircle,
  HelpCircle,
  Microscope,
  Sparkles,
  XCircle,
} from "lucide-react";

interface NmdPrediction {
  id: number;
  isoform_id: number;
  nmd_status: string;
  has_ptc: boolean;
  truncation_ratio?: number;
  orf_type?: string;
  ptc_position?: number;
  last_ejc_distance?: number;
  evidence: string;
  canonical_isoform_id?: number;
  truncation_severity?: string;
  nmd_confidence?: number;
  is_candidate_smorf?: boolean;
  is_candidate_uorf?: boolean;
}

const STATUS_CONFIG: Record<
  string,
  { icon: typeof AlertTriangle; color: string; bg: string; label: string; desc: string }
> = {
  nmd_target:             { icon: XCircle,        color: "text-red-400",       bg: "bg-red-500/10 border-red-500/30",       label: "NMD Target",             desc: "50-nt rule: stop codon >50 nt upstream of last EJC. Likely degraded by NMD." },
  nmd_escape:             { icon: CheckCircle,    color: "text-emerald-400",   bg: "bg-emerald-500/10 border-emerald-500/30", label: "NMD Escape",         desc: "50-nt rule: stop in last exon, near last EJC, or in the 5'-proximal escape zone." },
  non_translatable:       { icon: XCircle,        color: "text-foreground/80", bg: "bg-slate-500/10 border-slate-500/30",   label: "Non-Translatable",       desc: "No complete open reading frame detected — translation does not initiate." },
  predicted_translatable: { icon: CheckCircle,    color: "text-sky-400",       bg: "bg-sky-500/10 border-sky-500/30",       label: "Predicted Translatable", desc: "Heuristic-only: ORF looks complete; no PTC evidence found. NOT a positive NMD-escape call (we lack CDS-level data to rule it out)." },
  unknown:                { icon: HelpCircle,     color: "text-amber-400",     bg: "bg-amber-500/10 border-amber-500/30",   label: "Unknown",                desc: "Insufficient data to apply the 50-nt rule. Could be an NMD target, an escape, or neither." },
};

const TRUNC_CONFIG: Record<string, { color: string; bg: string; label: string }> = {
  intact:           { color: "text-emerald-400", bg: "bg-emerald-500/10 border-emerald-500/30", label: "Intact" },
  minor:            { color: "text-amber-400",   bg: "bg-amber-500/10 border-amber-500/30",     label: "Minor truncation" },
  severe:           { color: "text-rose-400",    bg: "bg-rose-500/10 border-rose-500/30",       label: "Severe truncation" },
  non_translatable: { color: "text-slate-400",   bg: "bg-slate-500/10 border-slate-500/30",     label: "Non-translatable" },
};

export default function NmdPredictionPanel({ isoformId }: { isoformId: number }) {
  const [pred, setPred] = useState<NmdPrediction | null>(null);
  const [loading, setLoading] = useState(true);
  const [notFound, setNotFound] = useState(false);

  useEffect(() => {
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/isoform/${isoformId}/nmd`)
      .then((r) => {
        if (r.status === 404) {
          setNotFound(true);
          return null;
        }
        return r.ok ? r.json() : null;
      })
      .then((d) => d && setPred(d))
      .catch(() => setNotFound(true))
      .finally(() => setLoading(false));
  }, [isoformId]);

  if (loading) return <div className="text-sm text-muted-foreground animate-pulse p-4">Loading NMD prediction…</div>;
  if (notFound || !pred)
    return (
      <div className="p-6 text-center text-muted-foreground">
        <HelpCircle className="w-8 h-8 mx-auto mb-2 opacity-40" />
        <p className="text-sm">No NMD prediction available.</p>
        <p className="text-xs mt-1 opacity-60 max-w-md mx-auto">
          NMD prediction needs the transcript&apos;s exon structure and stop-codon position; this
          isoform has no prediction stored.
        </p>
      </div>
    );

  const cfg = STATUS_CONFIG[pred.nmd_status] || STATUS_CONFIG.unknown;
  const Icon = cfg.icon;
  const truncCfg = pred.truncation_severity ? TRUNC_CONFIG[pred.truncation_severity] : null;

  // Candidate badges override the standard verdict visually
  const candidateBadge =
    pred.is_candidate_smorf
      ? { icon: Sparkles, label: "Candidate microprotein (peptidein)", color: "text-violet-400", bg: "bg-violet-500/10 border-violet-500/30" }
      : pred.is_candidate_uorf
      ? { icon: Microscope, label: "Candidate upstream ORF (uORF)", color: "text-fuchsia-400", bg: "bg-fuchsia-500/10 border-fuchsia-500/30" }
      : null;

  return (
    <div className="space-y-4 p-1">
      {/* Candidate banner */}
      {candidateBadge && (
        <div className={`rounded-xl border ${candidateBadge.bg} p-3 flex items-start gap-3`}>
          <candidateBadge.icon className={`w-5 h-5 shrink-0 mt-0.5 ${candidateBadge.color}`} />
          <div>
            <p className={`text-sm font-semibold ${candidateBadge.color}`}>{candidateBadge.label}</p>
            <p className="text-xs text-foreground/70 mt-0.5">
              The 50-nt rule is not validated for this class. NMD verdict deliberately shown as
              "Unknown".
            </p>
          </div>
        </div>
      )}

      {/* Status badge */}
      <div className={`rounded-xl border ${cfg.bg} p-5 flex items-start gap-4`}>
        <div className={`shrink-0 mt-0.5 ${cfg.color}`}>
          <Icon className="w-8 h-8" />
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2 flex-wrap">
            <h4 className={`text-lg font-bold ${cfg.color}`}>{cfg.label}</h4>
            {pred.nmd_confidence != null && (
              <span className="text-xs font-mono px-2 py-0.5 rounded bg-muted text-muted-foreground border border-border">
                confidence {(pred.nmd_confidence * 100).toFixed(0)}%
              </span>
            )}
          </div>
          <p className="text-sm text-foreground/80 mt-1">{cfg.desc}</p>
          {pred.nmd_confidence != null && (
            <div className="mt-3">
              <div className="w-full h-1.5 bg-muted rounded-full overflow-hidden">
                <div
                  className={`h-full ${cfg.color.replace("text-", "bg-")} transition-all`}
                  style={{ width: `${Math.min(pred.nmd_confidence * 100, 100)}%` }}
                />
              </div>
              <p className="text-[10px] text-muted-foreground/70 mt-1">
                Confidence ramps from 0% at distance=50 nt (rule threshold) to 100% at ≥200 nt.
              </p>
            </div>
          )}
        </div>
      </div>

      {/* Detail grid */}
      <div className="grid grid-cols-2 gap-3 text-xs">
        <div className="bg-muted/30 rounded-lg p-3 border border-border">
          <span className="text-muted-foreground block mb-1">Premature Stop Codon (CDS-confirmed)</span>
          <span className={`font-semibold ${pred.has_ptc ? "text-red-400" : "text-emerald-400"}`}>
            {pred.has_ptc ? "Yes" : "No"}
          </span>
        </div>
        <div className="bg-muted/30 rounded-lg p-3 border border-border">
          <span className="text-muted-foreground block mb-1">ORF Type</span>
          <span className="font-semibold text-foreground capitalize">{pred.orf_type ? pred.orf_type.replace(/_/g, " ") : "Unknown"}</span>
        </div>
      </div>

      {/* Orthogonal: truncation severity (NOT an NMD signal) */}
      {truncCfg && (
        <div className={`rounded-lg border ${truncCfg.bg} p-3`}>
          <div className="flex items-center justify-between text-xs mb-2">
            <span className={`font-semibold ${truncCfg.color}`}>
              Protein length vs canonical — {truncCfg.label}
            </span>
            {pred.truncation_ratio != null && (
              <span className="font-mono font-semibold text-foreground">
                {(pred.truncation_ratio * 100).toFixed(1)}%
              </span>
            )}
          </div>
          {pred.truncation_ratio != null && (
            <div className="w-full h-2 bg-slate-700 rounded-full overflow-hidden">
              <div
                className={`h-full rounded-full ${truncCfg.color.replace("text-", "bg-")}`}
                style={{ width: `${Math.min(pred.truncation_ratio * 100, 100)}%` }}
              />
            </div>
          )}
          <p className="text-[10px] text-muted-foreground mt-2 leading-snug">
            Independent of NMD. A short protein can have a legitimate alternative TSS/poly-A; a
            full-length protein can still carry a PTC. Both signals matter.
          </p>
        </div>
      )}

      {/* Evidence */}
      <div className="text-xs text-muted-foreground bg-muted/20 rounded-lg p-3 border border-border/50">
        <span className="font-semibold text-foreground/80 block mb-1">Evidence</span>
        {pred.evidence}
      </div>
    </div>
  );
}
