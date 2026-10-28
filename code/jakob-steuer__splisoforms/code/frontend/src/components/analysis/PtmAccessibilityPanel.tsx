"use client";
import { useEffect, useState, useMemo } from "react";
import { spliceAwareAlign } from "@/utils/alignment";
import { Droplets, ChevronDown, ChevronUp, Info, HelpCircle } from "lucide-react";
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

interface PtmSite {
  id: number;
  isoform_id: number;
  uniprot_id: string;
  residue_position: number;
  novel_residue_position?: number;
  residue_type: string;
  ptm_type: string;
  source_db: string;
  canonical_sasa?: number;
  novel_sasa?: number;
  delta_sasa?: number;
  accessibility_change: string;
  canonical_isoform_id?: number;
}

// ── Color helpers ─────────────────────────────────────────────────────────────

function ptmTypeColor(type: string): string {
  const t = type.toLowerCase();
  if (t.includes("phospho")) return "#3B82F6"; // blue-500
  if (t.includes("ubiquit")) return "#F59E0B"; // amber-500
  if (t.includes("acetyl"))  return "#8B5CF6"; // violet-500
  return "#6B7280"; // slate-500
}

// Keys match the accessibility_change values stored in the DB:
// exposed / buried / intact / undetermined / lost. ("unchanged" is a legacy
// alias for intact, kept so older rows still resolve to a colour.)
const ACC_COLOR: Record<string, string> = {
  exposed:      "#10B981", // emerald-500
  buried:       "#EF4444", // red-500
  intact:       "#6B7280", // slate-500 — retained, no meaningful change
  undetermined: "#F59E0B", // amber-500 — SASA unreliable / low pLDDT
  lost:         "#F97316", // orange-500
  unchanged:    "#6B7280", // legacy alias for intact
};

const CHANGE_COLORS: Record<string, { bg: string; text: string }> = {
  exposed:      { bg: "bg-emerald-500/15", text: "text-emerald-700 dark:text-emerald-400" },
  buried:       { bg: "bg-red-500/15",     text: "text-red-700 dark:text-red-400" },
  intact:       { bg: "bg-slate-500/15",   text: "text-slate-700 dark:text-foreground/80" },
  undetermined: { bg: "bg-amber-500/15",   text: "text-amber-700 dark:text-amber-400" },
  lost:         { bg: "bg-orange-500/15",  text: "text-orange-700 dark:text-orange-400" },
  unchanged:    { bg: "bg-slate-500/15",   text: "text-slate-700 dark:text-foreground/80" },
};

// ── Lollipop track ────────────────────────────────────────────────────────────

// Three stem-height tiers to stagger overlapping sites (px from baseline to circle top).
const TIER_H = [18, 32, 46]; // px
const HEAD_D = 10;            // circle diameter (px)
const TRACK_H = 62;           // total container height (px)
const BASELINE_B = 14;        // px from container bottom to baseline

// Minimum proximity (as fraction of barLen) that triggers stagger.
const STAGGER_THRESH = 0.025;

function LollipopTrack({
  sites,
  barLen,
  colorMode,
  showLost = true,
  hoveredResidues,
  selectedResidues,
  onHover,
  onSelect
}: {
  sites: PtmSite[];
  barLen: number;
  colorMode: "ptm-type" | "accessibility";
  /** When false, lost sites are hidden (use in the novel track). */
  showLost?: boolean;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelect?: (residues: number[]) => void;
}) {
  const visible = showLost
    ? sites
    : sites.filter(s => s.accessibility_change !== "lost");

  // Sort by canonical position and assign stagger tiers.
  const sorted = [...visible].sort((a, b) => a.residue_position - b.residue_position);
  const tiers = new Array<number>(sorted.length).fill(0);
  for (let i = 1; i < sorted.length; i++) {
    const pi = sorted[i].residue_position / barLen;
    const used = new Set<number>();
    for (let j = i - 1; j >= 0; j--) {
      if ((pi - sorted[j].residue_position / barLen) > STAGGER_THRESH) break;
      used.add(tiers[j]);
    }
    let t = 0;
    while (used.has(t)) t++;
    tiers[i] = Math.min(t, TIER_H.length - 1);
  }

  return (
    <div className="relative w-full overflow-visible" style={{ height: TRACK_H }}>
      {/* Backbone */}
      <div
        className="absolute left-0 right-0"
        style={{
          bottom: BASELINE_B,
          height: 2,
          backgroundColor: "rgba(148,163,184,0.35)",
          borderRadius: 1,
        }}
      />
      {/* aa-range labels */}
      <div
        className="absolute left-0 right-0 flex justify-between px-0.5"
        style={{ bottom: 2 }}
      >
        <span className="text-[9px] text-muted-foreground/50 select-none">aa 1</span>
        <span className="text-[9px] text-muted-foreground/50 select-none">aa {barLen}</span>
      </div>

      {/* Lollipops */}
      {sorted.map((site, i) => {
        const leftPct = Math.min(Math.max((site.residue_position / barLen) * 100, 0), 100);
        const stemH    = TIER_H[tiers[i]];
        const isLost   = site.accessibility_change === "lost";
        const color    =
          colorMode === "ptm-type"
            ? ptmTypeColor(site.ptm_type)
            : (ACC_COLOR[site.accessibility_change] ?? ACC_COLOR.unchanged);

        const tipLine1 = `${site.ptm_type.replace(/_/g, " ")} · ${site.residue_type}${site.residue_position}`;
        const tipLine2 = isLost
          ? "Lost in this isoform (spliced out)"
          : colorMode === "accessibility"
            ? `${site.accessibility_change}${
                site.novel_residue_position
                  ? ` · novel pos ${site.novel_residue_position}`
                  : ""
              }${
                site.delta_sasa != null
                  ? ` · ΔSASA ${site.delta_sasa > 0 ? "+" : ""}${site.delta_sasa.toFixed(1)} Å²`
                  : ""
              }`
            : site.source_db;

        const targetRes = site.novel_residue_position ?? site.residue_position;
        const isHovered = hoveredResidues?.includes(targetRes);
        const isSelected = selectedResidues?.includes(targetRes);
        const ringClass = isSelected ? "ring-2 ring-primary ring-offset-1" : isHovered ? "ring-1 ring-primary/80" : "";

        return (
          <div
            key={`${site.id}-${colorMode}`}
            title={`${tipLine1}\n${tipLine2}`}
            className={`group transition-all rounded-full ${ringClass}`}
            style={{
              position: "absolute",
              bottom: BASELINE_B + 2,
              left: `${leftPct}%`,
              transform: "translateX(-50%)",
              display: "flex",
              flexDirection: "column",
              alignItems: "center",
              cursor: "pointer",
            }}
            onMouseEnter={() => !isLost && onHover && onHover([targetRes])}
            onMouseLeave={() => !isLost && onHover && onHover([])}
            onClick={() => {
              if (!isLost && onSelect) {
                const isCurrentlySelected = selectedResidues && selectedResidues.length > 0 && selectedResidues[0] === targetRes;
                onSelect(isCurrentlySelected ? [] : [targetRes]);
              }
            }}
          >
            {/* ✕ only in accessibility mode; canonical track always shows circles */}
            {colorMode === "accessibility" && isLost ? (
              <svg
                width={HEAD_D}
                height={HEAD_D}
                viewBox="0 0 10 10"
                style={{ flexShrink: 0 }}
              >
                <line x1="1.5" y1="1.5" x2="8.5" y2="8.5" stroke={color} strokeWidth="2" strokeOpacity="0.80" />
                <line x1="8.5" y1="1.5" x2="1.5" y2="8.5" stroke={color} strokeWidth="2" strokeOpacity="0.80" />
              </svg>
            ) : (
              <div
                style={{
                  width: HEAD_D,
                  height: HEAD_D,
                  borderRadius: "50%",
                  backgroundColor: color,
                  opacity: 0.85,
                  flexShrink: 0,
                }}
              />
            )}
            {/* Stem */}
            <div
              style={{
                width: 2,
                height: stemH - HEAD_D,
                backgroundColor: color,
                opacity: isLost ? 0.40 : 0.65,
                flexShrink: 0,
              }}
            />
          </div>
        );
      })}
    </div>
  );
}

// ── Legend items ──────────────────────────────────────────────────────────────

function LegendCircle({ color, label }: { color: string; label: string }) {
  return (
    <div className="flex items-center gap-1">
      <div
        style={{ width: 8, height: 8, borderRadius: "50%", backgroundColor: color, opacity: 0.85 }}
      />
      <span className="text-[10px] text-muted-foreground">{label}</span>
    </div>
  );
}

function LegendX({ color, label }: { color: string; label: string }) {
  return (
    <div className="flex items-center gap-1">
      <svg width="9" height="9" viewBox="0 0 9 9">
        <line x1="1" y1="1" x2="8" y2="8" stroke={color} strokeWidth="2" strokeOpacity="0.85" />
        <line x1="8" y1="1" x2="1" y2="8" stroke={color} strokeWidth="2" strokeOpacity="0.85" />
      </svg>
      <span className="text-[10px] text-muted-foreground">{label}</span>
    </div>
  );
}

// ── Main component ────────────────────────────────────────────────────────────

export default function PtmAccessibilityPanel({
  isoformId,
  isCanonical = false,
  canonicalSeq,
  isoformSeq,
  hoveredResidues,
  selectedResidues,
  onHover,
  onSelect
}: {
  isoformId: number;
  isCanonical?: boolean;
  canonicalSeq?: string;
  isoformSeq?: string;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelect?: (residues: number[]) => void;
}) {
  const [sites, setSites]     = useState<PtmSite[]>([]);
  const [loading, setLoading] = useState(true);
  const [expanded, setExpanded] = useState(false);

  useEffect(() => {
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/isoform/${isoformId}/ptm`)
      .then(r => r.ok ? r.json() : [])
      .then(setSites)
      .catch(() => setSites([]))
      .finally(() => setLoading(false));
  }, [isoformId]);

  // ── Hooks must all be called before any conditional return ────────────────

  // Same NW-based coverage ribbon used in DomainsPanel
  const coverageBins = useMemo<("present" | "deleted" | "outside")[] | null>(() => {
    if (!canonicalSeq || !isoformSeq) return null;
    // spliceAwareAlign chains exact k-mer anchors, so deletion blocks land
    // cleanly at exon boundaries rather than being scattered by NW.
    const { canonAln, queryAln } = spliceAwareAlign(canonicalSeq, isoformSeq);
    const alnLen = canonAln.length;
    let firstIso = alnLen, lastIso = -1;
    for (let i = 0; i < alnLen; i++) {
      if (queryAln[i] !== "-") { if (firstIso === alnLen) firstIso = i; lastIso = i; }
    }
    const perResidue: ("present" | "deleted" | "outside")[] = [];
    for (let i = 0; i < alnLen; i++) {
      if (canonAln[i] === "-") continue;
      if (queryAln[i] !== "-") perResidue.push("present");
      else if (i < firstIso || i > lastIso) perResidue.push("outside");
      else perResidue.push("deleted");
    }
    const N = Math.min(200, perResidue.length || 1);
    const counts = Array.from({ length: N }, () => ({ present: 0, deleted: 0, outside: 0 }));
    const cLen = perResidue.length;
    perResidue.forEach((s, i) => {
      const bin = Math.min(Math.floor((i * N) / cLen), N - 1);
      counts[bin][s]++;
    });
    let lastStatus: "present" | "deleted" | "outside" = "outside";
    return counts.map(c => {
      let status: "present" | "deleted" | "outside";
      if (c.present === 0 && c.deleted === 0 && c.outside === 0) {
        status = lastStatus;
      } else if (c.present >= c.deleted && c.present >= c.outside) {
        status = "present";
      } else if (c.deleted >= c.present && c.deleted >= c.outside) {
        status = "deleted";
      } else {
        status = "outside";
      }
      lastStatus = status;
      return status;
    });
  }, [canonicalSeq, isoformSeq]);

  // ── Conditional returns (after all hooks) ─────────────────────────────────

  if (loading) {
    return <div className="text-sm text-muted-foreground animate-pulse p-4">Loading PTM sites…</div>;
  }
  if (sites.length === 0) {
    return (
      <div className="p-6 text-center text-muted-foreground">
        <Droplets className="w-8 h-8 mx-auto mb-2 opacity-40" />
        <p className="text-sm">No PTM sites to show for this isoform.</p>
        <p className="text-xs mt-1 opacity-60 max-w-md mx-auto">
          PTM sites are taken from the gene&apos;s reviewed UniProt entry, and we track only
          phosphorylation, ubiquitination and acetylation. Genes with no annotated site of those
          types — or with no reviewed accession — have nothing to display.
        </p>
      </div>
    );
  }

  // The bar always spans canonical sequence length — residue_position is canonical coords.
  const barLen = Math.max(...sites.map(s => s.residue_position), 1);

  // Summary counts
  const grouped = sites.reduce<Record<string, number>>((acc, s) => {
    acc[s.accessibility_change] = (acc[s.accessibility_change] || 0) + 1;
    return acc;
  }, {});

  const hasNovel = !isCanonical;

  return (
    <div className="space-y-4 p-1">

      {/* ── Summary badges ───────────────────────────────────────── */}
      <div className="flex flex-wrap gap-2">
        {Object.entries(grouped).map(([change, count]) => {
          const c = CHANGE_COLORS[change] || CHANGE_COLORS.unchanged;
          return (
            <span
              key={change}
              className={`${c.bg} ${c.text} px-3 py-1 rounded-full text-xs font-medium capitalize`}
            >
              {count} {change}
            </span>
          );
        })}
      </div>

      {/* ── Lollipop map ─────────────────────────────────────────── */}
      <div className="space-y-1">

        {/* Track header: canonical reference */}
        <div className="flex items-center justify-between gap-2 flex-wrap">
          <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-slate-500 dark:text-slate-400">
            {hasNovel ? "Canonical PTM Reference" : "PTM Site Map"}
            <span className="font-normal normal-case text-muted-foreground/50">
              {barLen} aa
            </span>
            <TooltipProvider>
              <Tooltip>
                <TooltipTrigger>
                  <Info className="w-3 h-3 opacity-40 cursor-help" />
                </TooltipTrigger>
                <TooltipContent className="bg-popover border-border text-xs max-w-xs shadow-xl">
                  <p className="text-foreground/90 font-semibold mb-1">PTM lollipop map</p>
                  <p className="text-foreground/80">
                    Each lollipop represents a known post-translational modification (PTM) site
                    on the canonical sequence. Hover over a site for position and SASA details.
                  </p>
                  {hasNovel && (
                    <p className="text-foreground/80 mt-1">
                      The upper track shows only PTM sites retained in this isoform, colored
                      by modification type. The lower track shows all sites — retained ones
                      colored by accessibility change, lost ones (spliced out) as ✕ marks.
                    </p>
                  )}
                </TooltipContent>
              </Tooltip>
            </TooltipProvider>
          </div>

          {/* Legend — PTM types */}
          <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
            <LegendCircle color="#3B82F6" label="Phosphorylation" />
            <LegendCircle color="#F59E0B" label="Ubiquitination" />
            <LegendCircle color="#8B5CF6" label="Acetylation" />
          </div>
        </div>

        {/* Canonical reference — all sites as circles colored by PTM type (no ✕ marks) */}
        <div className="bg-muted/20 rounded-lg border border-border px-3 pt-2 pb-1">
          <LollipopTrack
            sites={sites}
            barLen={barLen}
            colorMode="ptm-type"
            showLost={true}
            hoveredResidues={isCanonical ? hoveredResidues : undefined}
            selectedResidues={isCanonical ? selectedResidues : undefined}
            onHover={isCanonical ? onHover : undefined}
            onSelect={isCanonical ? onSelect : undefined}
          />
        </div>

        {/* ── Novel isoform PTM track ───────────────────────────── */}
        {hasNovel && (
          <>
            {/* Coverage ribbon (same as DomainsPanel) */}
            {coverageBins && (
              <div className="flex flex-col gap-0.5 mt-1 px-0.5">
                <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-muted-foreground/60">
                  Isoform Coverage
                  <TooltipProvider>
                    <Tooltip>
                      <TooltipTrigger><Info className="w-3 h-3 opacity-40 cursor-help" /></TooltipTrigger>
                      <TooltipContent className="bg-popover border-border text-xs max-w-xs shadow-xl">
                        <p className="text-foreground/90 font-semibold mb-1">Sequence alignment of isoform onto canonical</p>
                        <p className="text-foreground/80"><span className="text-emerald-400 font-bold">Green</span> — canonical region present in isoform</p>
                        <p className="text-foreground/80 mt-0.5"><span className="text-rose-400 font-bold">Red</span> — internal deletion (missing exon)</p>
                        <p className="text-foreground/80 mt-0.5"><span className="text-muted-foreground">Gray</span> — outside isoform's aligned extent</p>
                      </TooltipContent>
                    </Tooltip>
                  </TooltipProvider>
                </div>
                <div className="h-4 rounded overflow-hidden border border-border/40 flex w-full">
                  {(() => {
                    const segs: { status: string; count: number }[] = [];
                    for (const s of coverageBins) {
                      const last = segs[segs.length - 1];
                      if (last && last.status === s) last.count++;
                      else segs.push({ status: s, count: 1 });
                    }
                    return segs.map((seg, i) => (
                      <div key={i} style={{
                        width: `${(seg.count / coverageBins.length) * 100}%`,
                        backgroundColor:
                          seg.status === "present"  ? "rgba(52,211,153,0.45)"
                          : seg.status === "deleted" ? "rgba(251,113,133,0.50)"
                          : "rgba(100,116,139,0.12)",
                      }} />
                    ));
                  })()}
                </div>
              </div>
            )}

            {/* Divider */}
            <div className="flex items-center gap-1 my-0.5">
              <div className="flex-1 h-px bg-border/40" />
              <span className="text-[9px] text-muted-foreground/50 uppercase tracking-wider px-1">
                novel isoform
              </span>
              <div className="flex-1 h-px bg-border/40" />
            </div>

            {/* Novel track header */}
            <div className="flex items-center justify-between gap-2 flex-wrap">
              <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-sky-600 dark:text-[#65CBF3]">
                Isoform PTM Accessibility
                <span className="font-normal normal-case text-muted-foreground/50">
                  canonical positions
                </span>
              </div>

              {/* Legend — accessibility */}
              <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
                <LegendCircle color="#10B981" label="Exposed" />
                <LegendCircle color="#EF4444" label="Buried" />
                <LegendCircle color="#6B7280" label="Intact" />
                <LegendCircle color="#F59E0B" label="Undetermined" />
                <LegendX     color="#F97316" label="Lost" />
              </div>
            </div>

            {/* Novel lollipop track — all sites shown; lost = orange ✕, present = accessibility color */}
            <div className="bg-muted/30 rounded-lg border border-border px-3 pt-2 pb-1">
              <LollipopTrack
                sites={sites}
                barLen={barLen}
                colorMode="accessibility"
                showLost={true}
                hoveredResidues={hoveredResidues}
                selectedResidues={selectedResidues}
                onHover={onHover}
                onSelect={onSelect}
              />
            </div>
          </>
        )}
      </div>

      {/* ── Expandable table ─────────────────────────────────────── */}
      <button
        onClick={() => setExpanded(!expanded)}
        className="flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground/90 transition mt-2"
      >
        {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
        {expanded ? "Hide" : "Show"} PTM details ({sites.length} sites)
      </button>

      {expanded && (
        <div className="overflow-x-auto max-h-[400px] overflow-y-auto">
          <table className="w-full text-xs">
            <thead className="sticky top-0 bg-card">
              <tr className="border-b border-border text-muted-foreground">
                <th className="text-left py-2 px-2">Pos (canonical)</th>
                <th className="text-left py-2 px-2">Residue</th>
                <th className="text-left py-2 px-2">PTM Type</th>
                {/* On the canonical isoform there is nothing to compare against: it IS the
                    reference. Its own SASA is stored in novel_sasa (canonical_sasa is NULL),
                    so show a single correctly-labelled column and drop the comparison ones. */}
                {isCanonical ? (
                  <th className="text-right py-2 px-2">
                    <span className="inline-flex items-center justify-end gap-0.5">
                      SASA
                      <InfoTip text="Solvent-accessible surface area (Å²) in the canonical structure. Higher = more exposed." />
                    </span>
                  </th>
                ) : (
                  <>
                    <th className="text-right py-2 px-2">Novel Pos</th>
                    <th className="text-right py-2 px-2">
                      <span className="inline-flex items-center justify-end gap-0.5">
                        Canon. SASA
                        <InfoTip text="Solvent-accessible surface area (Å²) in the canonical structure. Higher = more exposed." />
                      </span>
                    </th>
                    <th className="text-right py-2 px-2">
                      <span className="inline-flex items-center justify-end gap-0.5">
                        Novel SASA
                        <InfoTip text="Solvent-accessible surface area (Å²) in this isoform. Compare with canonical." />
                      </span>
                    </th>
                    <th className="text-right py-2 px-2">
                      <span className="inline-flex items-center justify-end gap-0.5">
                        Δ SASA
                        <InfoTip text="Change in accessibility (isoform − canonical). + = more exposed, − = more buried." />
                      </span>
                    </th>
                    <th className="text-center py-2 px-2">
                      <span className="inline-flex items-center gap-0.5">
                        Change
                        <InfoTip text="Direction of change — gain (newly exposed) or loss (newly buried) can add or remove a PTM site." />
                      </span>
                    </th>
                  </>
                )}
              </tr>
            </thead>
            <tbody>
              {sites.map(s => {
                const c = CHANGE_COLORS[s.accessibility_change] || CHANGE_COLORS.unchanged;
                const targetRes = s.novel_residue_position ?? s.residue_position;
                const isLost = s.accessibility_change === "lost";
                const isSelected = selectedResidues && selectedResidues.length > 0 && selectedResidues[0] === targetRes;
                return (
                  <tr 
                    key={s.id} 
                    className={`border-b border-border/50 transition-colors cursor-pointer ${
                      isSelected ? "bg-muted/60" : "hover:bg-muted/30"
                    }`}
                    onMouseEnter={() => !isLost && onHover && onHover([targetRes])}
                    onMouseLeave={() => !isLost && onHover && onHover([])}
                    onClick={() => {
                      if (!isLost && onSelect) {
                        onSelect(isSelected ? [] : [targetRes]);
                      }
                    }}
                  >
                    <td className="py-2 px-2 font-mono">{s.residue_position}</td>
                    <td className="py-2 px-2 font-mono font-bold">{s.residue_type}</td>
                    <td className="py-2 px-2 capitalize">{s.ptm_type.replace(/_/g, " ")}</td>
                    {isCanonical ? (
                      // The canonical's own SASA lives in novel_sasa (see header comment).
                      <td className="py-2 px-2 text-right font-mono">
                        {(s.novel_sasa ?? s.canonical_sasa)?.toFixed(1) ?? "—"}
                      </td>
                    ) : (
                      <>
                        <td className="py-2 px-2 text-right font-mono">
                          {s.novel_residue_position ?? "—"}
                        </td>
                        <td className="py-2 px-2 text-right font-mono">
                          {s.canonical_sasa?.toFixed(1) ?? "—"}
                        </td>
                        <td className="py-2 px-2 text-right font-mono">
                          {s.novel_sasa?.toFixed(1) ?? "—"}
                        </td>
                        <td className="py-2 px-2 text-right font-mono">
                          {s.delta_sasa != null
                            ? (s.delta_sasa > 0 ? "+" : "") + s.delta_sasa.toFixed(1)
                            : "—"}
                        </td>
                        <td className="py-2 px-2 text-center">
                          <span
                            className={`${c.bg} ${c.text} px-2 py-0.5 rounded text-[10px] capitalize`}
                          >
                            {s.accessibility_change}
                          </span>
                        </td>
                      </>
                    )}
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
