import React from 'react';
import { Badge } from '@/components/ui/badge';
import GeneExonArchitecture from '@/app/gene/[id]/components/GeneExonArchitecture';

export interface JunctionEvent {
  id: number;
  junction_type: string;
  subcategory?: string;
  exon_count_novel?: number;
  exon_count_canonical?: number;
  transcript_pos?: number;
  is_novel: boolean;
}

interface SequenceMiniMapProps {
  sequence: string;
  junctionEvents?: JunctionEvent[];
  geneData?: any;
  currentIsoformId?: number;
}

export function SequenceMiniMap({ sequence, junctionEvents = [], geneData, currentIsoformId }: SequenceMiniMapProps) {
  if (!junctionEvents || junctionEvents.length === 0) {
    return null; // Return nothing if there are no junction events
  }

  // Usually there's only one SQANTI3 primary event per isoform
  const event = junctionEvents[0];

  const formatText = (text?: string) => {
    if (!text) return "Unknown";
    return text.replace(/_/g, " ");
  };

  const getExonCount = (iso: any) => {
    if (!iso || !iso.exons_json) return null;
    try {
      return JSON.parse(iso.exons_json).length;
    } catch {
      return null;
    }
  };

  const canonicalIso = geneData?.isoforms?.find((i: any) => i.is_canonical);
  const currentIso = geneData?.isoforms?.find((i: any) => i.id === currentIsoformId);
  const canonicalId = canonicalIso?.id;

  const selectedIds = [];
  if (canonicalId) selectedIds.push(canonicalId);
  if (currentIsoformId && currentIsoformId !== canonicalId) selectedIds.push(currentIsoformId);

  // Use actual parsed exon blocks if available, fallback to SQANTI3
  const actualNovelExonCount = (currentIso ? getExonCount(currentIso) : null) ?? event.exon_count_novel;
  const actualCanonExonCount = (canonicalIso ? getExonCount(canonicalIso) : null) ?? event.exon_count_canonical;

  const deltaExons = (actualNovelExonCount ?? 0) - (actualCanonExonCount ?? 0);
  const deltaColor = deltaExons > 0 ? "text-blue-400" : deltaExons < 0 ? "text-rose-400" : "text-emerald-400";
  const deltaSign = deltaExons > 0 ? "+" : "";

  return (
    <div className="w-full mt-6">
      <div className="flex items-center gap-2 mb-3">
        <h4 className="text-sm font-semibold text-foreground/90 uppercase tracking-widest flex items-center gap-2">
          Transcript Splicing Architecture
        </h4>
        {event.is_novel && (
          <Badge className="bg-primary/20 text-primary border-none shadow-sm text-[10px] h-5 py-0">Novel Splice Pattern</Badge>
        )}
      </div>

      <div className="bg-muted/40 border border-border/60 rounded-xl p-5 flex flex-col gap-6">
        
        {/* Top half: SQANTI3 Classification and Counts */}
        <div className="flex flex-col md:flex-row md:items-center justify-between gap-6">
        
        {/* Classification */}
        <div className="flex-1">
          <div className="flex flex-col gap-1">
            <span className="text-[10px] uppercase tracking-widest text-muted-foreground font-bold">Structural Classification</span>
            <span className="text-lg font-medium text-foreground capitalize">
              {formatText(event.junction_type)}
            </span>
            <span className="text-sm text-foreground/70 capitalize mt-1 border-l-2 border-border pl-3 py-0.5">
              {formatText(event.subcategory)}
            </span>
          </div>
        </div>

        {/* Exon Counts */}
        <div className="flex items-center gap-6 md:px-8 md:border-l border-border/50">
           <div className="flex flex-col items-center">
             <div className="text-3xl font-mono text-muted-foreground/60">{actualCanonExonCount ?? "?"}</div>
             <div className="text-[10px] uppercase font-semibold text-muted-foreground mt-1 tracking-wider">Canonical Exons</div>
           </div>
           
           <div className="h-8 w-px bg-border/50 hidden md:block"></div>
           
           <div className="flex flex-col items-center relative">
             <div className="text-3xl font-mono text-foreground">{actualNovelExonCount ?? "?"}</div>
             <div className="text-[10px] uppercase font-semibold text-foreground/80 mt-1 tracking-wider">Isoform Exons</div>
             {deltaExons !== 0 && actualNovelExonCount != null && actualCanonExonCount != null && (
               <div className={`absolute -top-3 -right-6 text-xs font-bold ${deltaColor}`}>
                 {deltaSign}{deltaExons}
               </div>
             )}
           </div>
        </div>

        </div>
      </div>
    </div>
  );
}
