"use client";
import { useState, useMemo } from "react";
import { Shield, ChevronDown, ChevronUp, Info } from "lucide-react";
import { spliceAwareAlign } from "@/utils/alignment";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";

export interface PfamDomain {
  id: number;
  pfam_accession: string;
  pfam_name: string;
  seq_start: number;
  seq_end: number;
  hmm_start: number;
  hmm_end: number;
  hmm_length: number;
  domain_coverage: number;
  e_value: number;
  domain_status: string;
  functional_class?: string | null;
  alteration_type?: string | null;
}

// Display labels for the curated functional classes (see
// backend/app/analysis/functional_classes.py).  Order is the legend order
// in the UI; "other" is rendered last so curated classes stand out.
const FUNCTIONAL_CLASS_LABELS: Record<string, string> = {
  catalytic_kinase:      "Kinase",
  catalytic_phosphatase: "Phosphatase",
  catalytic_protease:    "Protease",
  catalytic_other:       "Catalytic (other)",
  ub_ligase:             "Ubiquitin ligase / E2",
  dna_binding:           "DNA-binding",
  rna_binding:           "RNA-binding",
  signaling_adapter:     "Signaling adapter",
  oligomerization:       "Oligomerization / scaffold",
  other:                 "Other / unclassified",
};

const CATALYTIC_CLASSES = new Set([
  "catalytic_kinase", "catalytic_phosphatase",
  "catalytic_protease", "catalytic_other", "ub_ligase",
]);

// Human-readable alteration-geometry labels (replaces bare "truncated")
const ALTERATION_LABELS: Record<string, string> = {
  intact:              "Intact",
  truncated_n:         "N-terminal portion lost",
  truncated_c:         "C-terminal portion lost",
  truncated_both:      "Both termini truncated",
  internal_disruption: "Internal disruption",
  truncated_unknown:   "Truncated (geometry unknown)",
  lost:                "Lost",
  gained:              "Gained",
  unknown:             "Unknown",
};

export interface TedDomain {
  id: number;
  ted_id: string;
  domain_type: string;
  seq_start: number;
  seq_end: number;
  score?: number;
}

const PFAM_COLORS: Record<string, { bg: string; text: string; border: string }> = {
  present:   { bg: "bg-emerald-500/15", text: "text-emerald-700 dark:text-emerald-400", border: "border-emerald-500/30" },
  gained:    { bg: "bg-blue-500/15",    text: "text-blue-700 dark:text-blue-400",       border: "border-blue-500/30" },
  lost:      { bg: "bg-red-500/15",     text: "text-red-700 dark:text-red-400",         border: "border-red-500/30" },
  truncated: { bg: "bg-amber-500/15",   text: "text-amber-700 dark:text-amber-400",     border: "border-amber-500/30" },
};

const TED_COLORS = { bg: "bg-[#FF7D45]/15", text: "text-orange-700 dark:text-[#FF7D45]", border: "border-[#FF7D45]/40" };

export default function DomainsPanel({
  pfamDomains = [],
  tedDomains = [],
  seqLength,
  canonicalPfamDomains,
  canonicalSeqLength,
  isoformSeq,
  canonicalSeq,
  hoveredResidues,
  selectedResidues,
  onHover,
  onSelect,
}: {
  pfamDomains?: PfamDomain[];
  tedDomains?: TedDomain[];
  seqLength: number;
  canonicalPfamDomains?: PfamDomain[];
  canonicalSeqLength?: number;
  isoformSeq?: string;
  canonicalSeq?: string;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelect?: (residues: number[]) => void;
}) {
  const [expanded, setExpanded] = useState(false);

  // Alignment-based coverage in canonical coordinate space.
  // Returns:
  //   bins      – 200-bin array classifying each canonical region as
  //               'present' | 'deleted' | 'outside'
  //   isoToCan  – per-isoform-residue (0-based) → canonical 0-based position,
  //               used to re-project isoform domain boxes into canonical space
  const alignmentResult = useMemo<{
    bins: ("present" | "deleted" | "outside")[];
    isoToCan: number[];
  } | null>(() => {
    if (!canonicalSeq || !isoformSeq || !canonicalSeqLength) return null;

    // spliceAwareAlign chains exact k-mer anchors, so deletion blocks land
    // cleanly at exon boundaries rather than being scattered by NW.
    const { canonAln, queryAln } = spliceAwareAlign(canonicalSeq, isoformSeq);
    const alnLen = canonAln.length;

    // Build isoToCan: isoform 0-based residue index → canonical 0-based position.
    // Insertions (canonAln[i] === '-') map to the current canonical position.
    const isoToCan: number[] = [];
    let cPos = 0, iPos = 0;
    for (let i = 0; i < alnLen; i++) {
      if (canonAln[i] !== "-" && queryAln[i] !== "-") {
        // match / mismatch
        isoToCan[iPos] = cPos;
        cPos++; iPos++;
      } else if (canonAln[i] === "-") {
        // isoform insertion – no canonical residue; pin to adjacent canonical pos
        isoToCan[iPos] = cPos;
        iPos++;
      } else {
        // canonical residue, isoform gap – canonical advances only
        cPos++;
      }
    }

    // Find first and last alignment position where isoform is not a gap
    let firstIso = alnLen, lastIso = -1;
    for (let i = 0; i < alnLen; i++) {
      if (queryAln[i] !== "-") {
        if (firstIso === alnLen) firstIso = i;
        lastIso = i;
      }
    }

    // Build per-residue classification in canonical coordinate space
    const perResidue: ("present" | "deleted" | "outside")[] = [];
    for (let i = 0; i < alnLen; i++) {
      if (canonAln[i] === "-") continue; // isoform insertion – no canonical residue
      if (queryAln[i] !== "-") {
        perResidue.push("present");
      } else if (i < firstIso || i > lastIso) {
        perResidue.push("outside");
      } else {
        perResidue.push("deleted");
      }
    }

    // Bin into N display bins (dominant status wins).
    // Use canonicalSeq.length (= canonBarWidth's source of truth) as the bin
    // denominator so bin boundaries align 1:1 with percentage-positioned
    // domain boxes in the canonical track (both divide by the same length).
    const cLen = canonicalSeq.length; // same value canonBarWidth will use
    const N = Math.min(200, cLen || 1);
    const counts = Array.from({ length: N }, () => ({ present: 0, deleted: 0, outside: 0 }));
    perResidue.forEach((s, i) => {
      // FIX: use (i * N) / cLen before flooring to prevent empty bins due to floating point precision when N == cLen
      const bin = Math.min(Math.floor((i * N) / cLen), N - 1);
      counts[bin][s]++;
    });

    let lastStatus: "present" | "deleted" | "outside" = "outside";
    const bins = counts.map(c => {
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

    return { bins, isoToCan };
  }, [canonicalSeq, isoformSeq, canonicalSeqLength]);

  if (pfamDomains.length === 0 && tedDomains.length === 0 && (!canonicalPfamDomains || canonicalPfamDomains.length === 0)) {
    return (
      <div className="p-6 text-center text-muted-foreground">
        <Shield className="w-8 h-8 mx-auto mb-2 opacity-40" />
        <p className="text-sm">No Pfam or TED domain annotations available.</p>
        <p className="text-xs mt-1 opacity-60">Run the domain calculation pipelines mathematically to register domains.</p>
      </div>
    );
  }

  const coverageBins = alignmentResult?.bins ?? null;
  const isoToCan     = alignmentResult?.isoToCan ?? null;
  const showCoverage = coverageBins !== null && canonicalPfamDomains && canonicalPfamDomains.length > 0;

  // ── Coordinate spaces ────────────────────────────────────────────────────
  // canonBarWidth is the single source of truth for the canonical coordinate
  // axis. It must be consistent across:
  //   1. canonical domain box positions  (left = seq_start / canonBarWidth)
  //   2. coverage bin scale              (bin = floor(i / canonBarWidth * N))
  //   3. trapezoid connector fractions   (isoToCan[k] / canonBarWidth)
  //
  // We derive it from the actual canonicalSeq string when available, falling
  // back to the stored canonicalSeqLength integer, then barWidth. Using the
  // string length avoids drift when the stored integer and the fetched sequence
  // disagree (e.g. trailing whitespace stripped at ingest but not reflected in
  // the DB column).
  const canonBarWidth = canonicalSeq
    ? canonicalSeq.length
    : (canonicalSeqLength ?? Math.max(seqLength, 1));
  const isoBarWidth = seqLength > 0 ? seqLength : 1;
  const barWidth    = Math.max(isoBarWidth, canonBarWidth);

  const hasCanon = canonBarWidth > 0;
  const pct      = hasCanon ? Math.round((seqLength / canonBarWidth) * 100) : 100;

  // When blow-up layout is active (alignment + canonical track both present):
  //   • canonical track / coverage ribbon use canonBarWidth
  //   • isoform Pfam / TED tracks         use isoBarWidth
  //   • a trapezoid SVG connector bridges the two scales
  const useBlowup = showCoverage;

  // Fractional position [0,1] in canonical space where the isoform begins/ends.
  // Uses canonBarWidth so the fractions are on the same scale as domain boxes.
  const canonIsoStartFrac = (isoToCan && isoToCan.length > 0 && canonBarWidth > 0)
    ? isoToCan[0] / canonBarWidth
    : 0;
  const canonIsoEndFrac = (isoToCan && isoToCan.length > 0 && canonBarWidth > 0)
    ? Math.min((isoToCan[isoToCan.length - 1] + 1) / canonBarWidth, 1)
    : 1;

  return (
    <div className="space-y-4 p-1">

      <div className="flex gap-4 items-center justify-between">
         <div className="flex gap-3 text-xs flex-wrap">
           {canonicalPfamDomains && canonicalPfamDomains.length > 0 && (
             <span className="bg-slate-200 dark:bg-slate-500/20 text-slate-700 dark:text-slate-300 px-2 py-1 rounded border border-slate-400/70 dark:border-slate-400/40 capitalize">
               {canonicalPfamDomains.length} Canonical REF
             </span>
           )}
           {Object.entries(PFAM_COLORS).map(([status, c]) => {
             const count = pfamDomains.filter(d => d.domain_status === status).length;
             if (count === 0) return null;
             return (
               <span key={status} className={`${c.bg} ${c.text} px-2 py-1 rounded border ${c.border} capitalize`}>
                 {count} Pfam {status}
               </span>
             );
           })}
           {tedDomains.length > 0 && (() => {
              const chainsawCount = tedDomains.filter(d => d.domain_type === "CHAINSAW").length;
              const tedCount = tedDomains.length - chainsawCount;
              return (
                <span className={`${TED_COLORS.bg} ${TED_COLORS.text} px-2 py-1 rounded border ${TED_COLORS.border} capitalize`}>
                  {chainsawCount > 0 && tedCount === 0
                    ? `${chainsawCount} Chainsaw Domain${chainsawCount !== 1 ? "s" : ""}`
                    : chainsawCount === 0
                    ? `${tedCount} TED/TEDLH Region${tedCount !== 1 ? "s" : ""}`
                    : `${tedCount} TED + ${chainsawCount} Chainsaw`}
                </span>
              );
           })()}
           {hasCanon && pct !== 100 && (
             <span className="bg-muted text-muted-foreground px-2 py-1 rounded border border-border font-mono">
               {seqLength} aa · {pct}% of canonical
             </span>
           )}
         </div>
      </div>

      <div className="space-y-2">
         {/* Canonical Reference Track */}
         {canonicalPfamDomains && canonicalPfamDomains.length > 0 && (
           <div className="flex flex-col gap-1 mb-0">
             <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-slate-600 dark:text-slate-400">
               Canonical Reference (Pfam)
               <span className="font-normal normal-case text-muted-foreground/70">{canonBarWidth} aa</span>
             </div>
             <div className="relative h-10 bg-muted/60 rounded-lg overflow-hidden border border-border">
               {canonicalPfamDomains.map(d => {
                 const left = (d.seq_start / canonBarWidth) * 100;
                 const width = Math.max(((d.seq_end - d.seq_start) / canonBarWidth) * 100, 0.5);
                 return (
                   <div
                     key={`canon-pfam-${d.id}`}
                     className="absolute top-1 bottom-1 rounded bg-slate-400/50 dark:bg-slate-500/30 border border-slate-500/60 dark:border-slate-400/50 flex items-center justify-center overflow-hidden cursor-default"
                     style={{ left: `${left}%`, width: `${width}%`, minWidth: "4px" }}
                     title={`${d.pfam_name} (canonical aa ${d.seq_start}–${d.seq_end})`}
                   >
                     {width > 8 && (
                       <span className="text-[10px] font-medium truncate px-1 text-slate-800 dark:text-slate-300/80">
                         {d.pfam_name}
                       </span>
                     )}
                   </div>
                 );
               })}
             </div>

             {/* ── Alignment coverage ribbon ─────────────────────────────── */}
             {showCoverage && (
               <div className="flex flex-col gap-0.5 mt-0.5">
                 <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-muted-foreground">
                   Isoform Coverage
                   <TooltipProvider>
                     <Tooltip>
                       <TooltipTrigger><Info className="w-3 h-3 opacity-40 cursor-help" /></TooltipTrigger>
                       <TooltipContent className="bg-popover border-border text-xs max-w-xs shadow-xl">
                         <p className="text-foreground/90 font-semibold mb-1">Sequence alignment of isoform onto canonical</p>
                         <p className="text-foreground/80"><span className="text-emerald-400 font-bold">Green</span> — canonical region present in isoform (matched or substituted residues)</p>
                         <p className="text-foreground/80 mt-0.5"><span className="text-rose-400 font-bold">Red</span> — internal deletion: canonical region skipped by isoform (e.g. missing exon)</p>
                         <p className="text-foreground/80 mt-0.5"><span className="text-muted-foreground">Gray</span> — canonical region outside the isoform's aligned extent (N- or C-terminal difference)</p>
                       </TooltipContent>
                     </Tooltip>
                   </TooltipProvider>
                 </div>
                 <div
                   className="h-4 rounded overflow-hidden border border-border flex w-full"
                   title="Isoform vs canonical alignment coverage — hover segments for detail"
                 >
                   {(() => {
                     const segs: { status: string; count: number }[] = [];
                     for (const s of coverageBins!) {
                       const last = segs[segs.length - 1];
                       if (last && last.status === s) last.count++;
                       else segs.push({ status: s, count: 1 });
                     }
                     return segs.map((seg, i) => (
                       <div
                         key={i}
                         style={{
                           flexGrow: seg.count,
                           flexBasis: 0,
                           backgroundColor:
                             seg.status === "present"
                               ? "rgba(52,211,153,0.55)"   // emerald-400/55
                               : seg.status === "deleted"
                               ? "rgba(251,113,133,0.60)"  // rose-400/60
                               : "rgba(100,116,139,0.20)", // slate-500/20
                         }}
                       />
                     ));
                   })()}
                 </div>
                 <div className="flex justify-between text-[9px] text-muted-foreground/70 px-0.5">
                   <span>aa 1</span>
                   <span>aa {canonicalSeqLength}</span>
                 </div>
               </div>
             )}
           </div>
         )}

         {/* ── Blow-up connector OR fallback divider ──────────────────────── */}
         {canonicalPfamDomains && canonicalPfamDomains.length > 0 && (
           useBlowup ? (
             /* Trapezoid connector: shows how isoform coordinate space
                maps onto the canonical reference above */
             <div className="relative w-full" style={{ height: "24px" }}>
               <svg
                 className="absolute inset-0 w-full h-full"
                 preserveAspectRatio="none"
                 viewBox="0 0 1000 24"
               >
                 {/* Filled trapezoid */}
                 <polygon
                   points={`${canonIsoStartFrac * 1000},0 ${canonIsoEndFrac * 1000},0 1000,24 0,24`}
                   fill="rgba(100,116,139,0.10)"
                 />
                 {/* Left boundary line */}
                 <line
                   x1={canonIsoStartFrac * 1000} y1="0"
                   x2="0" y2="24"
                   stroke="rgba(71,85,105,0.55)"
                   strokeWidth="1"
                   strokeDasharray="4 3"
                 />
                 {/* Right boundary line */}
                 <line
                   x1={canonIsoEndFrac * 1000} y1="0"
                   x2="1000" y2="24"
                   stroke="rgba(71,85,105,0.55)"
                   strokeWidth="1"
                   strokeDasharray="4 3"
                 />
               </svg>
             </div>
           ) : (
             /* No alignment — simple visual divider */
             <div className="flex items-center gap-1 my-0.5">
               <div className="flex-1 h-px bg-border/70" />
               <span className="text-[9px] text-muted-foreground/80 uppercase tracking-wider px-1">isoform</span>
               <div className="flex-1 h-px bg-border/70" />
             </div>
           )
         )}

         {/* Pfam Track */}
         <div className="flex flex-col gap-1">
            <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-sky-600 dark:text-[#65CBF3]">
               Sequence Evolutionary Families (Pfam)
               {useBlowup && (
                 <span className="font-normal normal-case text-muted-foreground/70">{isoBarWidth} aa · isoform</span>
               )}
               <TooltipProvider>
                 <Tooltip>
                   <TooltipTrigger><Info className="w-3 h-3 opacity-50 cursor-help" /></TooltipTrigger>
                   <TooltipContent className="bg-popover border-border text-xs max-w-xs shadow-xl">
                      <p className="text-foreground/90">Biologically, sequence evolutionary domains (Pfam) can frequently overlap! This represents naturally nested domain architecture, identically repeating modular blocks, or parallel multi-clan HMM alignments over similar binding sites.</p>
                      <p className="text-muted-foreground mt-1">Lost domains are shown only in the canonical reference track above. Gaps in this track relative to the reference indicate domain loss.</p>
                   </TooltipContent>
                 </Tooltip>
               </TooltipProvider>
            </div>
            <div className="relative h-10 bg-muted/70 rounded-lg overflow-hidden border border-border">
               {pfamDomains.filter(d => d.domain_status !== "lost").map(d => {
                  const domW = useBlowup ? isoBarWidth : barWidth;
                  const left = (d.seq_start / domW) * 100;
                  const width = Math.max(((d.seq_end - d.seq_start) / domW) * 100, 0.5);
                  const colors = PFAM_COLORS[d.domain_status] || PFAM_COLORS.present;
                  const isTruncated = d.domain_status === "truncated";

                  // Determine which end(s) of the HMM model are missing.  Prefer
                  // the persisted alteration_type when available; fall back to
                  // computing from HMM coords for legacy rows that predate the
                  // c3e9a4b7d2f1 migration.
                  const altType = d.alteration_type || null;
                  const nTermCut = isTruncated && (
                    altType
                      ? altType === "truncated_n" || altType === "truncated_both"
                      : d.hmm_length > 0 && d.hmm_start > d.hmm_length * 0.10
                  );
                  const cTermCut = isTruncated && (
                    altType
                      ? altType === "truncated_c" || altType === "truncated_both"
                      : d.hmm_length > 0 && d.hmm_end < d.hmm_length * 0.90
                  );

                  // Build border-radius: flatten the truncated corner(s) to a square
                  // edge and keep the intact corner(s) rounded.
                  const borderRadius = isTruncated
                    ? `${nTermCut ? "0px" : "4px"} ${cTermCut ? "0px" : "4px"} ${cTermCut ? "0px" : "4px"} ${nTermCut ? "0px" : "4px"}`
                    : undefined;

                  // Gradient overlay for truncated ends (fade to transparent = "domain continues")
                  const gradientParts: string[] = [];
                  if (nTermCut) gradientParts.push("to right, rgba(0,0,0,0.25) 0%, transparent 25%");
                  if (cTermCut) gradientParts.push("to left, rgba(0,0,0,0.25) 0%, transparent 25%");
                  const gradientOverlay = gradientParts.length > 0
                    ? gradientParts.map(g => `linear-gradient(${g})`).join(", ")
                    : undefined;

                  const pctCov = Math.round(d.domain_coverage * 100);
                  const hmmInfo = d.hmm_length > 0
                    ? ` · HMM ${d.hmm_start}–${d.hmm_end}/${d.hmm_length} (${pctCov}% model coverage)`
                    : ` · ${pctCov}% model coverage`;
                  const altLabel = altType ? ALTERATION_LABELS[altType] ?? altType : d.domain_status;
                  const classLabel = d.functional_class
                    ? FUNCTIONAL_CLASS_LABELS[d.functional_class] ?? d.functional_class
                    : null;
                  const tooltip = `${d.pfam_name}${classLabel ? ` [${classLabel}]` : ""}`
                    + ` (isoform aa ${d.seq_start}–${d.seq_end}) — ${altLabel}${hmmInfo}`;

                  const range = Array.from({length: d.seq_end - d.seq_start + 1}, (_, i) => d.seq_start + i);
                  const isHovered = hoveredResidues?.some(r => range.includes(r));
                  const isSelected = selectedResidues?.some(r => range.includes(r));
                  const highlightClass = isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-background" : isHovered ? "ring-1 ring-primary/70" : "";

                  return (
                  <div
                     key={`pfam-${d.id}`}
                     className={`absolute top-1 bottom-1 ${colors.bg} border ${colors.border} flex items-center justify-center overflow-hidden cursor-pointer group transition-all ${highlightClass}`}
                     style={{
                       left: `${left}%`, width: `${width}%`, minWidth: "4px",
                       borderRadius,
                       // Dashed border on truncated side(s)
                       borderLeftStyle:  nTermCut ? "dashed" : undefined,
                       borderRightStyle: cTermCut ? "dashed" : undefined,
                     }}
                     title={tooltip}
                     onMouseEnter={() => onHover && onHover(range)}
                     onMouseLeave={() => onHover && onHover([])}
                     onClick={() => {
                       if (onSelect) {
                         const isCurrentlySelected = selectedResidues && selectedResidues.length > 0 && selectedResidues[0] === range[0];
                         onSelect(isCurrentlySelected ? [] : range);
                       }
                     }}
                  >
                     {/* Gradient fade overlay on truncated end(s) */}
                     {gradientOverlay && (
                       <div className="absolute inset-0 pointer-events-none" style={{ backgroundImage: gradientOverlay }} />
                     )}
                     {width > 8 && (
                       <span className={`relative text-[10px] font-medium truncate px-1 ${colors.text}`}>
                         {d.pfam_name}
                         {isTruncated && <span className="ml-0.5 opacity-60 text-[9px]">{pctCov}%</span>}
                       </span>
                     )}
                  </div>
                  );
               })}
               {pfamDomains.filter(d => d.domain_status !== "lost").length === 0 && (
                 <span className="absolute inset-0 flex items-center justify-center text-xs text-muted-foreground italic">
                   {pfamDomains.length > 0 ? "All detected domains lost — see canonical reference above" : "No Pfam Annotations"}
                 </span>
               )}
            </div>
         </div>

         {/* TED Track */}
         <div className="flex flex-col gap-1 mt-2">
            <div className="flex items-center gap-2 text-[10px] uppercase font-bold tracking-widest text-orange-600 dark:text-[#FF7D45]">
               Structure-Derived Domain Boundaries
               {useBlowup && tedDomains.length > 0 && (
                 <span className="font-normal normal-case text-muted-foreground/70">{isoBarWidth} aa · isoform</span>
               )}
               <TooltipProvider>
                 <Tooltip>
                   <TooltipTrigger><Info className="w-3 h-3 opacity-50 cursor-help" /></TooltipTrigger>
                   <TooltipContent className="bg-popover border-border text-xs max-w-xs shadow-xl">
                      <p className="text-foreground/90 font-semibold mb-1">Two annotation sources:</p>
                      <p className="text-foreground/80"><span className="font-mono text-[#FF7D45]">TED / TEDLH</span> — Canonical isoforms: consensus of 3 tools (ECOD, Chainsaw, Merizo) from the TED database. Includes CATH classification and pLDDT-validated boundaries.</p>
                      <p className="text-foreground/80 mt-1"><span className="font-mono text-[#FF7D45]">Chainsaw</span> — Novel long-read isoforms: single-tool deep-learning structural segmentation on AF3 predictions. Domain boundaries only, no structural taxonomy assigned.</p>
                   </TooltipContent>
                 </Tooltip>
               </TooltipProvider>
            </div>
            <div className="relative h-10 bg-muted/70 rounded-lg overflow-hidden border border-border">
               {tedDomains.map((d, i) => {
                  // Blow-up: isoform-native coordinate space.
                  const domW = useBlowup ? isoBarWidth : barWidth;
                  const left = (d.seq_start / domW) * 100;
                  const width = Math.max(((d.seq_end - d.seq_start) / domW) * 100, 0.5);
                  const isChainsaw = d.domain_type === "CHAINSAW";
                  const label = isChainsaw ? `Domain ${i + 1}` : d.domain_type;
                  const tooltip = isChainsaw
                    ? `Chainsaw domain ${i + 1} (isoform aa ${d.seq_start}–${d.seq_end}) — single-tool AF3 segmentation`
                    : `${d.ted_id} [${d.domain_type}] (isoform aa ${d.seq_start}–${d.seq_end})${d.score ? ` · score ${d.score.toFixed(3)}` : ""}`;
                  const range = Array.from({length: d.seq_end - d.seq_start + 1}, (_, j) => d.seq_start + j);
                  const isHovered = hoveredResidues?.some(r => range.includes(r));
                  const isSelected = selectedResidues?.some(r => range.includes(r));
                  const highlightClass = isSelected ? "ring-2 ring-orange-500 ring-offset-1 ring-offset-background" : isHovered ? "ring-1 ring-orange-400" : "";

                  return (
                  <div
                     key={`ted-${i}`}
                     className={`absolute top-1 bottom-1 rounded ${TED_COLORS.bg} border ${isChainsaw ? "border-dashed" : ""} ${TED_COLORS.border} flex items-center justify-center overflow-hidden cursor-pointer group backdrop-blur transition-all ${highlightClass}`}
                     style={{ left: `${left}%`, width: `${width}%`, minWidth: "4px" }}
                     title={tooltip}
                     onMouseEnter={() => onHover && onHover(range)}
                     onMouseLeave={() => onHover && onHover([])}
                     onClick={() => {
                       if (onSelect) {
                         const isCurrentlySelected = selectedResidues && selectedResidues.length > 0 && selectedResidues[0] === range[0];
                         onSelect(isCurrentlySelected ? [] : range);
                       }
                     }}
                  >
                     {width > 8 && <span className="text-[10px] font-medium truncate px-1 text-orange-900 dark:text-white">{label}</span>}
                  </div>
                  );
               })}
               {tedDomains.length === 0 && <span className="absolute inset-0 flex items-center justify-center text-xs text-muted-foreground italic">No Structure-Based Domain Boundaries</span>}
            </div>
         </div>
      </div>

      <button onClick={() => setExpanded(!expanded)} className="flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground transition mt-4 mb-2">
        {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
        {expanded ? "Collapse" : "Expand"} tabular domains details ({pfamDomains.length + tedDomains.length})
      </button>

      {expanded && (
         <div className="flex flex-col gap-6">
            {pfamDomains.length > 0 && (
               <div className="overflow-x-auto">
               <h4 className="text-xs font-bold text-foreground/80 mb-2 uppercase tracking-widest border-b border-border pb-1">Pfam Database Alignments</h4>
               <table className="w-full text-xs">
                  <thead>
                     <tr className="border-b border-border text-muted-foreground">
                     <th className="text-left py-2 px-2">Domain</th>
                     <th className="text-left py-2 px-2">Accession</th>
                     <th className="text-left py-2 px-2">Class</th>
                     <th className="text-right py-2 px-2">Range</th>
                     <th className="text-right py-2 px-2">Coverage</th>
                     <th className="text-center py-2 px-2">Alteration</th>
                     </tr>
                  </thead>
                  <tbody>
                     {pfamDomains.map(d => {
                     const c = PFAM_COLORS[d.domain_status] || PFAM_COLORS.present;
                     const isCatalytic = d.functional_class && CATALYTIC_CLASSES.has(d.functional_class);
                     const classLabel = d.functional_class
                       ? FUNCTIONAL_CLASS_LABELS[d.functional_class] ?? d.functional_class
                       : "—";
                     const altLabel = d.alteration_type
                       ? ALTERATION_LABELS[d.alteration_type] ?? d.alteration_type
                       : d.domain_status;
                     return (
                        <tr key={`t-pfam-${d.id}`} className="border-b border-border/70 hover:bg-muted/50 cursor-pointer" onClick={() => { if (onSelect) { const range = Array.from({length: d.seq_end - d.seq_start + 1}, (_, j) => d.seq_start + j); const isCurrentlySelected = selectedResidues && selectedResidues.length > 0 && selectedResidues[0] === range[0]; onSelect(isCurrentlySelected ? [] : range); } }}>
                           <td className="py-2 px-2 font-medium text-foreground">{d.pfam_name}</td>
                           <td className="py-2 px-2 font-mono text-foreground/80">{d.pfam_accession}</td>
                           <td className="py-2 px-2">
                             {d.functional_class && d.functional_class !== "other" ? (
                               <span className={`px-1.5 py-0.5 rounded text-[10px] border ${
                                 isCatalytic
                                   ? "bg-violet-500/10 text-violet-700 dark:text-violet-300 border-violet-500/30"
                                   : "bg-sky-500/10 text-sky-700 dark:text-sky-300 border-sky-500/30"
                               }`}>{classLabel}</span>
                             ) : (
                               <span className="text-muted-foreground/60 text-[10px]">{classLabel}</span>
                             )}
                           </td>
                           <td className="py-2 px-2 text-right font-mono">{d.seq_start}–{d.seq_end}</td>
                           <td className="py-2 px-2 text-right">{(d.domain_coverage * 100).toFixed(0)}%</td>
                           <td className="py-2 px-2 text-center">
                           <span className={`${c.bg} ${c.text} px-2 py-0.5 rounded text-[10px] border ${c.border}`}>{altLabel}</span>
                           </td>
                        </tr>
                     );
                     })}
                  </tbody>
               </table>
               </div>
            )}

            {tedDomains.length > 0 && (
               <div className="overflow-x-auto">
               <h4 className="text-xs font-bold text-[#FF7D45] mb-2 uppercase tracking-widest border-b border-border pb-1 flex items-center gap-2">
                  Structure-Derived Domain Boundaries
                  {tedDomains.some(d => d.domain_type === "CHAINSAW") && tedDomains.some(d => d.domain_type !== "CHAINSAW") ? (
                    <span className="text-[9px] text-muted-foreground normal-case bg-popover px-1 py-0.5 rounded">TED Consensus + Chainsaw</span>
                  ) : tedDomains.every(d => d.domain_type === "CHAINSAW") ? (
                    <span className="text-[9px] text-muted-foreground normal-case bg-popover px-1 py-0.5 rounded">Chainsaw · Single-tool AF3 Segmentation</span>
                  ) : (
                    <span className="text-[9px] text-muted-foreground normal-case bg-popover px-1 py-0.5 rounded">TED Consensus · 3-Tool Validated</span>
                  )}
               </h4>
               <table className="w-full text-xs">
                  <thead>
                     <tr className="border-b border-border text-muted-foreground">
                     <th className="text-left py-2 px-2">TED Topology ID</th>
                     <th className="text-left py-2 px-2">Classification Type</th>
                     <th className="text-right py-2 px-2">Structural Range</th>
                     <th className="text-right py-2 px-2">Confidence Score</th>
                     </tr>
                  </thead>
                  <tbody>
                     {tedDomains.map((d, i) => (
                        <tr key={`t-ted-${i}`} className="border-b border-border/70 hover:bg-muted/50 cursor-pointer" onClick={() => { if (onSelect) { const range = Array.from({length: d.seq_end - d.seq_start + 1}, (_, j) => d.seq_start + j); const isCurrentlySelected = selectedResidues && selectedResidues.length > 0 && selectedResidues[0] === range[0]; onSelect(isCurrentlySelected ? [] : range); } }}>
                           <td className="py-2 px-2 font-medium text-foreground/90 font-mono tracking-tight">{d.ted_id}</td>
                           <td className="py-2 px-2 text-muted-foreground font-bold tracking-widest">{d.domain_type}</td>
                           <td className="py-2 px-2 text-right font-mono">{d.seq_start}–{d.seq_end}</td>
                           <td className="py-2 px-2 text-right opacity-80">{d.score ? d.score.toFixed(3) : '-'}</td>
                        </tr>
                     ))}
                  </tbody>
               </table>
               </div>
            )}
         </div>
      )}
    </div>
  );
}
