"use client";

import { useMemo, useState, useEffect } from "react";
import { spliceAwareAlign } from "@/utils/alignment";
import { ChevronDown, ChevronUp } from "lucide-react";

interface Props {
  canonicalSequence: string;
  isoformSequence: string;
  canonicalLabel?: string;
  isoformLabel?: string;
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelectionChange?: (residues: number[]) => void;
  onClick?: (resNo: number) => void;
}

export default function IsoformCanonicalAlignment({
  canonicalSequence,
  isoformSequence,
  canonicalLabel = "Canonical",
  isoformLabel = "Isoform",
  hoveredResidues = [],
  selectedResidues = [],
  onHover,
  onSelectionChange,
  onClick,
}: Props) {
  const [expanded, setExpanded] = useState(true);
  const [selectionStart, setSelectionStart] = useState<number | null>(null);
  const [isDragging, setIsDragging] = useState(false);

  useEffect(() => {
    const handleGlobalMouseUp = () => setIsDragging(false);
    window.addEventListener("mouseup", handleGlobalMouseUp);
    return () => window.removeEventListener("mouseup", handleGlobalMouseUp);
  }, []);

  const handleResidueMouseDown = (e: React.MouseEvent, resNo: number) => {
    setIsDragging(true);
    setSelectionStart(resNo);

    if (e.shiftKey && selectionStart !== null) {
      const start = Math.min(selectionStart, resNo);
      const end = Math.max(selectionStart, resNo);
      const range = [];
      for (let i = start; i <= end; i++) range.push(i);
      
      if (e.metaKey || e.ctrlKey) {
         const combined = Array.from(new Set([...selectedResidues, ...range]));
         if (onSelectionChange) onSelectionChange(combined);
      } else {
         if (onSelectionChange) onSelectionChange(range);
      }
    } else if (e.metaKey || e.ctrlKey) {
      if (onSelectionChange) {
         if (selectedResidues.includes(resNo)) {
             onSelectionChange(selectedResidues.filter(r => r !== resNo));
         } else {
             onSelectionChange([...selectedResidues, resNo]);
         }
      }
    } else {
      if (onClick) onClick(resNo);
      if (onSelectionChange) onSelectionChange([resNo]);
    }
  };

  const handleResidueMouseEnter = (resNo: number) => {
    if (onHover) onHover([resNo]);
    
    if (isDragging && selectionStart !== null) {
      const start = Math.min(selectionStart, resNo);
      const end = Math.max(selectionStart, resNo);
      const range = [];
      for (let i = start; i <= end; i++) range.push(i);
      if (onSelectionChange) onSelectionChange(range);
    }
  };

  const { canonAln, queryAln, identity, isoPosArray } = useMemo(() => {
    // Splice-aware alignment: shared exon anchors ensure novel exon sequences
    // appear as green insertions rather than amber mismatches.
    const { canonAln, queryAln } = spliceAwareAlign(canonicalSequence, isoformSequence);

    // Identity = shared exon anchors matched / canonical length
    let matchCount = 0;
    const isoPosArray: (number | null)[] = [];
    let currIsoPos = 1;

    for (let i = 0; i < canonAln.length; i++) {
      if (canonAln[i] !== "-" && queryAln[i] !== "-" && canonAln[i] === queryAln[i]) {
        matchCount++;
      }
      if (queryAln[i] !== "-") {
        isoPosArray.push(currIsoPos++);
      } else {
        isoPosArray.push(null);
      }
    }
    const identity = canonicalSequence.length > 0
      ? (matchCount / canonicalSequence.length) * 100
      : 0;

    return { canonAln, queryAln, identity, isoPosArray };
  }, [canonicalSequence, isoformSequence]);

  return (
    <div className="bg-card/60 border border-primary/20 shadow-xl rounded-xl p-6 mt-4">
      {/* Header */}
      <div className="flex items-center justify-between mb-3">
        <div className="flex items-center gap-3">
          <h3 className="text-sm font-semibold text-foreground/90 uppercase tracking-widest">
            Sequence Alignment (Canonical vs. Isoform)
          </h3>
          <span className="text-xs bg-primary/10 text-primary border border-primary/20 px-2 py-0.5 rounded font-mono">
            {identity.toFixed(1)}% identity
          </span>
        </div>
        <button
          onClick={() => setExpanded(!expanded)}
          className="flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground transition"
        >
          {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
          {expanded ? "Collapse" : "Expand"}
        </button>
      </div>

      {expanded && (
        <>
          {/* Sticky legend — lives outside the scroll container so it's always visible */}
          <div className="sticky top-0 z-10 flex flex-wrap gap-4 mb-2 py-1.5 bg-card/95 backdrop-blur-sm border-b border-border/30 text-[10px] text-muted-foreground">
            <span className="flex items-center gap-1">
              <span className="text-foreground/40 font-mono font-bold">A</span> Shared exon (match)
            </span>
            <span className="flex items-center gap-1">
              <span className="text-amber-500 font-mono font-bold">A</span> SNP / point diff
            </span>
            <span className="flex items-center gap-1">
              <span className="text-rose-500 font-mono opacity-60">-</span> Skipped exon
            </span>
            <span className="flex items-center gap-1">
              <span className="text-emerald-500 font-mono font-bold">A</span> Novel / alternative exon
            </span>
          </div>

          <div className="overflow-x-auto custom-scrollbar pb-2">
            <div className="inline-block min-w-full">
              {/* Numbering row — must use text-xs (same as character rows) so slot widths match */}
              <div className="flex items-start gap-2 mb-1">
                <span className="shrink-0 w-28" />
                <div className="font-mono text-xs tracking-[0.2em] leading-relaxed whitespace-nowrap h-4 relative">
                  {queryAln.split("").map((_, i) => {
                    const num = isoPosArray[i];
                    const showLabel = num !== null && (num === 1 || num % 10 === 0);
                    const showTick  = num !== null && !showLabel && num % 5 === 0;
                    return (
                      <span key={`num-${i}`} className="inline-block w-[1.2em] text-center relative overflow-visible">
                        {showLabel ? (
                          <span className="absolute left-1/2 -translate-x-1/2 top-0 text-[9px] font-bold text-muted-foreground/70 leading-none">
                            {num}
                          </span>
                        ) : showTick ? (
                          <span className="text-[9px] text-muted-foreground/35 leading-none">·</span>
                        ) : null}
                      </span>
                    );
                  })}
                </div>
              </div>

              {/* Canonical row */}
              <div className="flex items-start gap-2 mb-1">
                <span className="shrink-0 text-[10px] font-bold bg-slate-700/50 text-slate-300 border border-slate-600/50 px-2 py-0.5 rounded uppercase tracking-wider w-28 text-center mt-0.5">
                  CANONICAL
                </span>
                <div className="font-mono text-xs tracking-[0.2em] leading-relaxed whitespace-nowrap">
                  {canonAln.split("").map((char, i) => (
                    <span
                      key={i}
                      className={`inline-block w-[1.2em] text-center ${
                        char === "-"
                          ? "text-muted-foreground/30"
                          : "text-foreground/70"
                      }`}
                    >
                      {char}
                    </span>
                  ))}
                </div>
              </div>

              {/* Isoform row */}
              <div className="flex items-start gap-2">
                <span
                  className="shrink-0 text-[10px] font-medium bg-primary/10 text-primary border border-primary/20 px-2 py-0.5 rounded w-28 text-center truncate mt-0.5"
                  title={isoformLabel}
                >
                  {isoformLabel}
                </span>
                <div className="font-mono text-xs tracking-[0.2em] leading-relaxed whitespace-nowrap">
                  {queryAln.split("").map((char, i) => {
                    const canonChar = canonAln[i];
                    const num = isoPosArray[i];
                    let cls = "";
                    if (char === "-") {
                      // Deletion relative to canonical
                      cls = "text-rose-500 opacity-60";
                    } else if (canonChar === "-") {
                      // Insertion not in canonical
                      cls = "text-emerald-500 font-bold";
                    } else if (char !== canonChar) {
                      // Mismatch
                      cls = "text-amber-500 font-bold";
                    } else {
                      // Match
                      cls = "text-foreground/40";
                    }
                    
                    const isHovered = num !== null && hoveredResidues.includes(num);
                    const isSelected = num !== null && selectedResidues.includes(num);
                    const highlightClass = isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-card bg-primary/20 rounded-sm z-10 font-bold" : isHovered ? "ring-1 ring-primary/80 bg-primary/10 rounded-sm z-10" : "";
                    
                    return (
                      <span 
                        key={i} 
                        className={`inline-block w-[1.2em] text-center cursor-pointer select-none relative ${cls} ${highlightClass}`}
                        onMouseDown={num !== null ? (e) => handleResidueMouseDown(e, num) : undefined}
                        onMouseEnter={num !== null ? () => handleResidueMouseEnter(num) : undefined}
                        onMouseLeave={() => {
                          if (onHover) onHover([]);
                        }}
                      >
                        {char}
                      </span>
                    );
                  })}
                </div>
              </div>
            </div>
          </div>
        </>
      )}
    </div>
  );
}
