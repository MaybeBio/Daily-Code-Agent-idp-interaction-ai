"use client";

import React, { useMemo } from "react";
import { GitCommit, AlignLeft, TriangleAlert } from "lucide-react";

interface Props {
  gene: any;
  selectedIsoformIds: number[];
  sortedIsoformIds?: number[];
  hideHeader?: boolean;
  hideCardBackground?: boolean;
}

const getIsCanonical = (iso: any, gene: any) => {
  return gene?.canonical_transcript_id ? iso.id === gene.canonical_transcript_id : !!iso.is_canonical;
};

export default function GeneExonArchitecture({ gene, selectedIsoformIds, sortedIsoformIds, hideHeader, hideCardBackground }: Props) {
  const { displayIsoforms, dynamicBounds, coordMismatch } = useMemo(() => {
    const isoById = new Map((gene.isoforms ?? []).map((i: any) => [i.id, i]));
    let isos: any[] = [];
    if (selectedIsoformIds.length > 0) {
      if (sortedIsoformIds && sortedIsoformIds.length > 0) {
        // Use the authoritative sort order from the table
        isos = sortedIsoformIds
          .filter(id => selectedIsoformIds.includes(id))
          .map(id => isoById.get(id))
          .filter(Boolean);
      } else {
        isos = (gene.isoforms?.filter((iso: any) => selectedIsoformIds.includes(iso.id)) || [])
          .sort((a: any, b: any) => {
            const aCanon = getIsCanonical(a, gene);
            const bCanon = getIsCanonical(b, gene);
            if (aCanon && !bCanon) return -1;
            if (!aCanon && bCanon) return 1;
            return 0;
          });
      }
    } else {
       const canon = gene.isoforms?.find((i: any) => getIsCanonical(i, gene));
       const others = (gene.isoforms?.filter((i: any) => !getIsCanonical(i, gene)) || []).slice(0, 2);
       isos = canon ? [canon, ...others] : others;
    }

    const sortedIsos = isos;

    // Parse all exons to dynamically determine visual bounds
    let minStart = Infinity;
    let maxEnd = -Infinity;
    
    const parsedIsoforms = sortedIsos.map((iso: any) => {
      let parsedExons: [number, number][] = [];
      try {
        parsedExons = iso.exons_json ? JSON.parse(iso.exons_json) : [];
      } catch (e) {
        console.error("Failed to parse exons for", iso.transcript_id);
      }
      
      let isFallback = false;
      if (parsedExons.length === 0) {
         if (iso.metadata_record?.genomic_start && iso.metadata_record?.genomic_end) {
             parsedExons = [[iso.metadata_record.genomic_start, iso.metadata_record.genomic_end]];
             isFallback = true;
         } else if (getIsCanonical(iso, gene) && gene.genomic_start && gene.genomic_end) {
             parsedExons = [[gene.genomic_start, gene.genomic_end]];
             isFallback = true;
         }
      }

      parsedExons.forEach(([s, e]) => {
        if (s < minStart) minStart = s;
        if (e > maxEnd) maxEnd = e;
      });

      return { ...iso, _parsedExons: parsedExons, _isFallback: isFallback };
    });

    // Detect coordinate mismatch: canonical and novel isoforms in non-overlapping ranges.
    // Threshold: if the gap between their ranges is > 50 kb we assume different assemblies.
    const canonIso = parsedIsoforms.find((i: any) => getIsCanonical(i, gene));
    const novelIsos = parsedIsoforms.filter((i: any) => !getIsCanonical(i, gene));
    let coordMismatch = false;
    if (canonIso && novelIsos.length > 0 && canonIso._parsedExons.length > 0) {
      const canonMin = Math.min(...canonIso._parsedExons.map(([s]: [number, number]) => s));
      const canonMax = Math.max(...canonIso._parsedExons.map(([, e]: [number, number]) => e));
      for (const nov of novelIsos) {
        if (nov._parsedExons.length === 0) continue;
        const novMin = Math.min(...nov._parsedExons.map(([s]: [number, number]) => s));
        const novMax = Math.max(...nov._parsedExons.map(([, e]: [number, number]) => e));
        const overlaps = novMin <= canonMax && novMax >= canonMin;
        const gap = overlaps ? 0 : Math.max(novMin - canonMax, canonMin - novMax);
        if (!overlaps && gap > 50_000) { coordMismatch = true; break; }
      }
    }

    // Add a 5% buffer on each side
    if (minStart !== Infinity && maxEnd !== -Infinity) {
      const buffer = Math.floor((maxEnd - minStart) * 0.05);
      return {
        displayIsoforms: parsedIsoforms,
        dynamicBounds: { start: minStart - buffer, end: maxEnd + buffer },
        coordMismatch,
      };
    }

    return {
      displayIsoforms: parsedIsoforms,
      dynamicBounds: { start: gene.genomic_start || 0, end: gene.genomic_end || 10000 },
      coordMismatch,
    };
  }, [gene, selectedIsoformIds, sortedIsoformIds]);

  if (!displayIsoforms.length) return null;

  const start = dynamicBounds.start;
  const end = dynamicBounds.end;
  const length = end - start;

  if (length <= 0) return null;

  const fmt = (n: number) => new Intl.NumberFormat("en-US").format(n);
  const chrStr = gene.chromosome.startsWith('chr') ? gene.chromosome : `chr${gene.chromosome}`;

  return (
    <div className={`${hideCardBackground ? "w-full flex flex-col" : "bg-card/60 backdrop-blur-md border border-primary/20 shadow-xl rounded-xl p-6 w-full flex flex-col"}`}>
      {/* Coordinate mismatch warning */}
      {coordMismatch && (
        <div className="flex items-start gap-2 mb-4 px-3 py-2.5 rounded-lg bg-amber-500/10 border border-amber-500/30 text-amber-400 text-xs">
          <TriangleAlert className="w-3.5 h-3.5 mt-0.5 flex-shrink-0" />
          <span>
            Genomic coordinates of novel isoform(s) do not overlap the canonical reference — the datasets may use different genome assemblies.
            Exon tracks are shown individually but are <strong>not directly comparable</strong> on this axis.
          </span>
        </div>
      )}

      {/* Legend - Always visible if there are multiple isoforms */}
      {displayIsoforms.length > 1 && (
         <div className="flex flex-wrap items-center gap-3 text-[10px] uppercase tracking-widest font-semibold px-3 py-1.5 bg-muted/30 rounded-md border border-border/50 w-fit mb-4">
            <span className="flex items-center gap-1.5"><div className="w-3 h-2 border border-emerald-400/50 bg-emerald-500/60"></div> Preserved Exon</span>
            <span className="flex items-center gap-1.5"><div className="w-3 h-2 border border-amber-500 bg-amber-400"></div> Novel Exon</span>
            <span className="flex items-center gap-1.5"><div className="w-3 h-2 border border-violet-600 bg-violet-500/80"></div> Intron Retention</span>
            <span className="flex items-center gap-1.5"><div className="w-3 h-2 border border-cyan-600 bg-cyan-500/80"></div> Alt Splice Site (&gt;5 bp)</span>
            <span className="flex items-center gap-1.5"><div className="w-3 h-2 border border-teal-500 bg-teal-400/80"></div> Micro-Alt Splice (1–5 bp)</span>
            <span className="flex items-center gap-1.5"><div className="w-3 h-2 border-2 border-dashed border-rose-500/60 bg-rose-500/10"></div> Skipped Exon</span>
         </div>
      )}
      
      <div className="flex flex-col gap-6 relative py-4">
         {/* Baseline Axis */}
         <div className="absolute top-0 bottom-0 left-[100px] border-l border-border/50" />
                  {displayIsoforms.map((iso: any, idx: number) => {
            const exons: [number, number][] = iso._parsedExons || [];
            const canonicalExons = displayIsoforms.find((i: any) => getIsCanonical(i, gene))?._parsedExons || [];

            return (
               <div key={iso.id} className="relative flex items-center gap-4 group">
                   <div className="w-24 flex-shrink-0 text-right pr-2 flex flex-col justify-center items-end gap-1">
                     <span className="text-[10px] font-mono font-medium text-muted-foreground group-hover:text-primary transition-colors block break-all">
                        {iso.transcript_name || iso.transcript_id}
                     </span>
                     {iso._isFallback && (
                        <span className="text-[8px] leading-tight text-amber-500 bg-amber-500/10 px-1 py-0.5 rounded border border-amber-500/20" title="Missing precise exon boundaries. Displaying full genomic span.">
                          Missing Exons
                        </span>
                     )}
                  </div>
                  
                  <div className="flex-1 relative h-6 bg-muted/30 rounded-md border border-border/50 overflow-hidden">
                     {/* Intron line background */}
                     <div className="absolute top-1/2 left-0 right-0 h-px bg-muted -translate-y-1/2" />
                     
                     {/* Ghost blocks: canonical exons absent in this novel isoform */}
                     {!getIsCanonical(iso, gene) && canonicalExons.map(([cStart, cEnd]: [number, number], cIdx: number) => {
                        const hasOverlap = exons.some(([ns, ne]) => ns <= cEnd && ne >= cStart);
                        if (hasOverlap) return null;
                        const ghostLeft  = ((cStart - start) / length) * 100;
                        const ghostWidth = ((cEnd - cStart) / length) * 100;
                        return (
                          <div
                            key={`ghost-${cIdx}`}
                            className="absolute top-1.5 bottom-1.5 border-2 border-dashed border-rose-500/60 bg-rose-500/10 z-[5]"
                            style={{
                              left:  `${Math.max(0, ghostLeft)}%`,
                              width: `${Math.max(0.5, ghostWidth)}%`,
                            }}
                            title={`Canonical exon ${cIdx + 1} (${cStart}–${cEnd}): absent in this isoform`}
                          />
                        );
                     })}

                     {/* Real Exon Blocks */}
                     {exons.map(([eStart, eEnd], eIdx) => {
                       const left = ((eStart - start) / length) * 100;
                       const width = ((eEnd - eStart) / length) * 100;

                       let highlightClass = "";
                       let statusText = `Exon ${eIdx + 1}: ${eStart}-${eEnd}`;

                       if (!getIsCanonical(iso, gene) && canonicalExons.length > 0) {
                          const overlaps = canonicalExons.filter(([cStart, cEnd]: [number, number]) => eStart <= cEnd && eEnd >= cStart);
                          if (overlaps.length === 0) {
                             highlightClass = "bg-amber-400 border-amber-500 z-10";
                             statusText += " (Novel Exon)";
                          } else if (overlaps.length > 1) {
                             highlightClass = "bg-violet-500/80 border-violet-600 z-10";
                             statusText += " (Intron Retention)";
                          } else {
                             const [cStart, cEnd] = overlaps[0];
                             const startDiff = Math.abs(eStart - cStart);
                             const endDiff   = Math.abs(eEnd - cEnd);
                             const maxDiff   = Math.max(startDiff, endDiff);
                             if (maxDiff > 5) {
                                // Major alternative splice site (>5 bp boundary shift)
                                highlightClass = "bg-cyan-500/80 border-cyan-600 z-10";
                                statusText += ` (Alternative Splice Site, Δ${maxDiff} bp)`;
                             } else if (maxDiff >= 1) {
                                // Micro alternative splice site (1–5 bp boundary shift)
                                highlightClass = "bg-teal-400/80 border-teal-500 z-10";
                                statusText += ` (Micro-Alt Splice, Δ${maxDiff} bp)`;
                             }
                          }
                       }

                       if (iso._isFallback) {
                           statusText += " (Missing Boundaries: Full Gene Span)";
                       }

                       return (
                         <div
                           key={eIdx}
                           className={`absolute top-1.5 bottom-1.5 shadow-sm border transition-all duration-300 ${
                             getIsCanonical(iso, gene)
                               ? (iso._isFallback ? 'bg-[repeating-linear-gradient(45deg,transparent,transparent_4px,rgba(99,102,241,0.2)_4px,rgba(99,102,241,0.2)_8px)] bg-indigo-500/30 border-indigo-500/50 border-dashed group-hover:bg-indigo-500/40' : 'bg-indigo-500/80 border-indigo-400 group-hover:bg-indigo-400')
                               : highlightClass ? highlightClass : 'bg-emerald-500/60 border-emerald-400/50 group-hover:bg-emerald-400'
                           }`}
                           style={{
                             left: `${Math.max(0, left)}%`,
                             width: `${Math.max(0.5, width)}%`
                           }}
                           title={statusText}
                         />
                       );
                     })}
                  </div>
               </div>
            );
         })}
      </div>
      
      <div className="flex items-center gap-4 py-2 px-2 mt-2 text-xs font-mono text-slate-500">
         <div className="w-24 shrink-0"></div>
         <div className="flex-1 flex justify-between mx-4">
            <span>{chrStr}:{fmt(start)}</span>
            <span>{fmt(end)} ({gene.strand})</span>
         </div>
      </div>
    </div>
  );
}
