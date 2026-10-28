"use client";

import React, { useMemo, useState, useEffect } from "react";
import { ListCollapse, Fingerprint, AlertTriangle } from "lucide-react";

interface SelectedIsoform {
  id: number;
  transcript_id: string;
  transcript_name?: string | null;
  sequence?: string;
}

interface Props {
  canonicalId?: string;
  canonicalSequence?: string;
  selectedIsoforms: SelectedIsoform[];
  hoveredResidues?: number[];
  selectedResidues?: number[];
  onHover?: (residues: number[]) => void;
  onSelectionChange?: (residues: number[]) => void;
  onClick?: (resNo: number) => void;
}

import { spliceAwareAlign, AlignedPair } from "@/utils/alignment";

/**
 * Detect whether a pairwise splice-aware alignment has any canonical overlap.
 *
 * spliceAwareAlign returns a specific pattern when no shared k-mers are found:
 *   canonAln = canonical + "-".repeat(novelLen)    (entire novel sequence as insertion)
 *   queryAln = "-".repeat(canonLen) + novelSeq     (entire canonical as deletion)
 *
 * We detect this by checking whether the first canonLen characters of queryAln
 * are all gaps — meaning the isoform covers none of the canonical sequence.
 */
function detectCanonOverlap(aln: AlignedPair, canonLen: number): boolean {
  if (canonLen === 0) return false;
  for (let i = 0; i < canonLen && i < aln.queryAln.length; i++) {
    if (aln.queryAln[i] !== "-") return true;
  }
  return false;
}

/**
 * Unifies multiple pairwise splice-aware alignments into a single MSA space.
 *
 * Each pairwise alignment is expressed relative to the canonical sequence.
 * Insertions (novel exon sequences) in any isoform widen the master canonical
 * row with gap columns, so all isoforms share the same column coordinate space.
 *
 * Also detects "contested" insertion columns — master gap slots where two or
 * more isoforms carry different residues.  Those columns are not biologically
 * co-aligned: each isoform's novel sequence just happens to occupy the same
 * canonical coordinate gap, but the residues have no positional correspondence.
 * Contested columns are flagged so the renderer can warn the user.
 */
function unifyMSA(
  canonical: string,
  queryAlignments: { id: string; aln: AlignedPair; hasCanonOverlap: boolean }[]
) {
  if (queryAlignments.length === 0) return { masterCanon: canonical, blocks: [], contestedCols: new Set<number>() };

  // 1. Map each alignment to a set of "insertions before canonical index i"
  // indexed 0 to canonical.length (canonical.length is for trailing insertions)
  const allInsertions = queryAlignments.map(({ id, aln, hasCanonOverlap }) => {
    const insertions: string[] = new Array(canonical.length + 1).fill("");
    let cIdx = 0;
    let currentInsertion = "";

    for (let i = 0; i < aln.canonAln.length; i++) {
      if (aln.canonAln[i] === "-") {
        currentInsertion += aln.queryAln[i];
      } else {
        insertions[cIdx] = currentInsertion;
        currentInsertion = "";
        cIdx++;
      }
    }
    insertions[cIdx] = currentInsertion; // Tail
    return { id, insertions, queryAln: aln.queryAln, canonAln: aln.canonAln, hasCanonOverlap };
  });

  // 2. Determine master gap widths for each canonical position
  const masterGapWidths = new Array(canonical.length + 1).fill(0);
  allInsertions.forEach(insObj => {
    insObj.insertions.forEach((str, i) => {
      masterGapWidths[i] = Math.max(masterGapWidths[i], str.length);
    });
  });

  // 3. Construct Master Canonical
  let masterCanon = "";
  for (let i = 0; i <= canonical.length; i++) {
    masterCanon += "-".repeat(masterGapWidths[i]);
    if (i < canonical.length) masterCanon += canonical[i];
  }

  // 4. Project each isoform into the master coordinate space
  const blocks = allInsertions.map(insObj => {
    let projectedQuery = "";
    const aln = insObj.canonAln;
    const qAln = insObj.queryAln;

    let cIdx = 0;
    let alnPos = 0;

    for (let i = 0; i <= canonical.length; i++) {
      // Add current insertion segment (padded to master gap width).
      //
      // An insertion before canonical position i is anchored on its RIGHT to
      // canonical[i] (the residue that immediately follows it), so isoforms that
      // share sequence adjacent to that anchor (e.g. two N-terminal extensions
      // whose shorter one is a suffix of the longer) co-align only when the
      // insertions are RIGHT-aligned. Left-padding scattered the shared residues
      // into different columns. Trailing insertions (i === canonical.length) have
      // no following canonical residue — they anchor on the LEFT, so keep padEnd.
      const gapWidth = masterGapWidths[i];
      const myIns = insObj.insertions[i];
      projectedQuery += i < canonical.length
        ? myIns.padStart(gapWidth, "-")
        : myIns.padEnd(gapWidth, "-");

      if (i < canonical.length) {
        // Advance past insertion columns in this isoform's pairwise alignment
        while (alnPos < aln.length && aln[alnPos] === "-") alnPos++;
        if (alnPos < aln.length) {
          projectedQuery += qAln[alnPos];
          alnPos++;
        } else {
          projectedQuery += "-";
        }
      }
    }
    return { id: insObj.id, queryAln: projectedQuery, hasCanonOverlap: insObj.hasCanonOverlap };
  });

  // 5. Find contested insertion columns.
  //
  //    A master gap column (masterCanon[col] === '-') is "contested" when two or
  //    more isoforms place *different* residues there.  Those residues are not
  //    co-aligned — they are simply different novel sequences that both happen to
  //    sit in the same canonical coordinate gap.  Reading them column-by-column as
  //    if they were aligned is misleading and should be flagged visually.
  const contestedCols = new Set<number>();
  if (blocks.length > 1) {
    for (let col = 0; col < masterCanon.length; col++) {
      if (masterCanon[col] !== "-") continue; // only insertion (gap) columns
      const seen = new Set<string>();
      for (const block of blocks) {
        const ch = block.queryAln[col];
        if (ch && ch !== "-") seen.add(ch);
      }
      if (seen.size > 1) contestedCols.add(col);
    }
  }

  return { masterCanon, blocks, contestedCols };
}

export default function GeneSequenceAlignment({
  canonicalId,
  canonicalSequence = "",
  selectedIsoforms,
  hoveredResidues = [],
  selectedResidues = [],
  onHover,
  onSelectionChange,
  onClick
}: Props) {
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

  const msa = useMemo(() => {
    if (!canonicalSequence) return { masterCanon: "", blocks: [], contestedCols: new Set<number>(), canonPosArray: [] };

    // Pairwise splice-aware alignment for each novel isoform vs canonical.
    // spliceAwareAlign uses k-mer anchors (shared exon segments) so novel exon
    // sequences appear as green insertions, not amber mutation columns.
    // We also detect whether each isoform shares ANY sequence with the canonical.
    const rawPairs = selectedIsoforms
      .filter(iso => iso.transcript_id !== canonicalId)
      .map(iso => {
        const aln = spliceAwareAlign(canonicalSequence, iso.sequence || "");
        return {
          id: iso.transcript_id,
          aln,
          hasCanonOverlap: detectCanonOverlap(aln, canonicalSequence.length),
        };
      });

    const unified = unifyMSA(canonicalSequence, rawPairs);

    // 5. Generate canonical positions mapping for the master sequence
    let cPos = 1;
    const canonPosArray = unified.masterCanon.split('').map(char => {
       if (char !== '-') return cPos++;
       return null;
    });

    return { ...unified, canonPosArray };
  }, [canonicalSequence, selectedIsoforms, canonicalId]);

  // transcript_id → human-readable transcript_name (falls back to the id itself).
  // block.id stays the stable transcript_id used for keys/matching; we only swap
  // the visible label.
  const nameById = useMemo(() => {
    const m = new Map<string, string>();
    selectedIsoforms.forEach(iso => m.set(iso.transcript_id, iso.transcript_name || iso.transcript_id));
    return m;
  }, [selectedIsoforms]);

  const hasAnyNoOverlap  = msa.blocks.some(b => !b.hasCanonOverlap);
  const hasContestedCols = msa.contestedCols.size > 0;

  // Legend item
  const Legend = ({ color, label }: { color: string; label: string }) => (
    <span className="flex items-center gap-1 text-[10px] text-muted-foreground">
      <span className={`inline-block w-3 h-3 rounded-sm ${color}`} />
      {label}
    </span>
  );

  return (
    <div className="bg-card/60 backdrop-blur-md border border-primary/20 shadow-xl rounded-xl p-6 w-full flex flex-col">
      <div className="flex items-center justify-between gap-3 mb-6 pb-4 border-b border-border/50">
        <div className="flex items-center gap-3">
          <div className="p-2 bg-primary/10 rounded-lg border border-primary/20">
            <ListCollapse className="w-5 h-5 text-primary" />
          </div>
          <div>
            <p className="text-sm text-muted-foreground">
              Splice-aware alignment — shared exon anchors guide novel isoforms against the canonical reference
            </p>
          </div>
        </div>
        {/* Legend */}
        <div className="flex items-center gap-3 flex-wrap shrink-0 justify-end">
          <Legend color="bg-foreground/20" label="shared exon" />
          <Legend color="bg-emerald-500/70" label="novel / alternative exon" />
          <span className="flex items-center gap-1 text-[10px] text-muted-foreground">
            <span className="font-mono text-rose-400 font-bold opacity-70 text-sm leading-none">–</span>
            skipped exon
          </span>
          <Legend color="bg-amber-500/70" label="SNP / point difference" />
          {hasContestedCols && (
            <Legend color="bg-orange-400/30 border border-orange-400/40" label="independent novel seqs" />
          )}
        </div>
      </div>

      {/* Warning callouts */}
      {hasAnyNoOverlap && (
        <div className="mb-4 flex items-start gap-2 bg-amber-500/5 border border-amber-500/20 rounded-lg px-4 py-3 text-xs text-amber-600 dark:text-amber-400">
          <AlertTriangle className="w-4 h-4 mt-0.5 flex-shrink-0" />
          <span>
            <strong>One or more isoforms share no sequence anchors with the canonical reference.</strong>{" "}
            Their novel exon regions are displayed relative to the canonical coordinate space but are{" "}
            <em>not</em> co-aligned to each other — residues in the same column do not correspond.
            The structures may still be similar if the isoforms share exons with each other.
          </span>
        </div>
      )}
      {!hasAnyNoOverlap && hasContestedCols && (
        <div className="mb-4 flex items-start gap-2 bg-orange-400/5 border border-orange-400/20 rounded-lg px-4 py-3 text-xs text-orange-600 dark:text-orange-400">
          <AlertTriangle className="w-4 h-4 mt-0.5 flex-shrink-0" />
          <span>
            <strong>Highlighted columns contain independent novel insertions.</strong>{" "}
            Multiple isoforms have different novel sequences at the same canonical gap position.
            These residues are not co-aligned — columns are shared by canonical coordinate only.
          </span>
        </div>
      )}

      {selectedIsoforms.length === 0 ? (
        <div className="flex-1 flex flex-col items-center justify-center text-center opacity-60 py-12">
          <Fingerprint className="w-12 h-12 text-slate-500 mb-4" />
          <p className="text-slate-400 max-w-sm">Select one or more transcripts from the table to perform a dynamic sequence alignment.</p>
        </div>
      ) : (
        <div className="flex max-w-full overflow-x-auto pb-4 custom-scrollbar bg-muted/50 rounded-lg p-4 border border-border">
          <div className="flex flex-col gap-6 min-w-max">

            {/* Master Canonical Reference Row */}
            <div className="flex items-center gap-4">
               <div className="w-32 flex-shrink-0 text-right">
                  <span className="text-xs font-bold text-foreground uppercase tracking-widest bg-muted px-2 py-1 rounded">Canonical</span>
               </div>
               <div className="font-mono text-sm tracking-[0.2em] bg-muted px-4 py-2 rounded-lg border border-border shadow-inner flex">
                  {msa.masterCanon.split('').map((char, charIdx) => {
                    const num = msa.canonPosArray[charIdx];
                    const isHovered = num !== null && hoveredResidues.includes(num);
                    const isSelected = num !== null && selectedResidues.includes(num);
                    const highlightClass = isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-card bg-primary/20 rounded-sm z-10 font-bold" : isHovered ? "ring-1 ring-primary/80 bg-primary/10 rounded-sm z-10" : "";
                    
                    return (
                      <span
                        key={charIdx}
                        className={`inline-block w-[1.2em] text-center cursor-pointer select-none relative ${char === '-' ? 'text-muted-foreground/30' : 'text-foreground'} ${msa.contestedCols.has(charIdx) ? 'bg-orange-400/20' : ''} ${highlightClass}`}
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

            {/* Isoform Rows */}
            {msa.blocks.map((block) => (
              <div key={block.id} className="flex items-center gap-4 group">
                  <div className="w-32 flex-shrink-0 text-right space-y-1">
                    <span className="text-[10px] font-mono font-medium text-primary bg-primary/10 border border-primary/20 px-2 py-1 rounded truncate block" title={block.id}>
                       {nameById.get(block.id) || block.id}
                    </span>
                    {(() => {
                        // Identity = shared exon anchors matched / canonical length
                        const matchesCount = block.queryAln.split('').filter((char, idx) => char !== '-' && char === msa.masterCanon[idx]).length;
                        const canonLength = msa.masterCanon.replace(/-/g, '').length;
                        const identity = canonLength > 0 ? (matchesCount / canonLength) * 100 : 0;
                        return (
                          <span className="text-[9px] font-bold text-muted-foreground uppercase tracking-tighter">
                            {identity.toFixed(1)}% ID
                          </span>
                        );
                    })()}
                    {!block.hasCanonOverlap && (
                      <span className="text-[9px] font-bold text-amber-500 uppercase tracking-tighter flex items-center gap-0.5">
                        <AlertTriangle className="w-2.5 h-2.5 inline" /> no shared anchors
                      </span>
                    )}
                  </div>
                 <div className="font-mono text-sm tracking-[0.2em] flex bg-muted px-4 py-2 rounded-lg border border-border shadow-inner relative">
                    {block.queryAln.split('').map((char, charIdx) => {
                       const canonChar = msa.masterCanon[charIdx];
                       const isInsertion  = canonChar === '-';   // novel exon in isoform
                       const isDeletion   = char === '-';         // skipped exon
                       const isMatch      = char === canonChar;
                       const isContested  = msa.contestedCols.has(charIdx);

                       let colorClass = "text-foreground/50"; // shared exon (dim)
                       if (isInsertion && char !== '-') {
                         colorClass = isContested
                           ? "text-emerald-500 font-bold bg-orange-400/20"   // contested novel
                           : "text-emerald-500 font-bold bg-emerald-500/10"; // uncontested novel
                       } else if (isDeletion) {
                         colorClass = "text-rose-500 font-bold opacity-40";
                       } else if (!isMatch && char !== '-') {
                         colorClass = "text-amber-500 font-bold bg-amber-500/10";
                       }

                       const num = msa.canonPosArray[charIdx];
                       const isHovered = num !== null && hoveredResidues.includes(num);
                       const isSelected = num !== null && selectedResidues.includes(num);
                       const highlightClass = isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-card bg-primary/20 rounded-sm z-10 font-bold" : isHovered ? "ring-1 ring-primary/80 bg-primary/10 rounded-sm z-10" : "";

                       return (
                         <span
                           key={charIdx}
                           className={`${colorClass} inline-block w-[1.2em] text-center rounded-sm transition-colors duration-200 cursor-pointer select-none relative ${highlightClass}`}
                           title={isContested ? "Independent novel sequence — not co-aligned with other isoforms" : undefined}
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
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
