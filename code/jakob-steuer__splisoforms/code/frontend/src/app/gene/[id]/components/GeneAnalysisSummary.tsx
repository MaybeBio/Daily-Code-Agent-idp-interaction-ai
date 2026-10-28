"use client";

import React, { useMemo } from "react";
import { Activity, ShieldAlert, Cpu, AlertTriangle, MinusCircle, PlusCircle, Link } from "lucide-react";

interface Props {
  gene: any;
  selectedIsoformIds: number[];
  sortedIsoformIds?: number[];
}

export default function GeneAnalysisSummary({ gene, selectedIsoformIds, sortedIsoformIds }: Props) {
  const displayIsoforms = useMemo(() => {
    if (selectedIsoformIds.length > 0) {
      const isoById = new Map((gene.isoforms ?? []).map((i: any) => [i.id, i]));
      if (sortedIsoformIds && sortedIsoformIds.length > 0) {
        // Follow the authoritative table sort order
        return sortedIsoformIds
          .filter(id => selectedIsoformIds.includes(id))
          .map(id => isoById.get(id))
          .filter(Boolean);
      }
      // Fallback: canonical first
      return (gene.isoforms?.filter((iso: any) => selectedIsoformIds.includes(iso.id)) || [])
        .sort((a: any, b: any) => {
          if (a.is_canonical !== b.is_canonical) return a.is_canonical ? -1 : 1;
          return 0;
        });
    }
    // Default: Show Canonical and top 3 by length
    const canon = gene.isoforms?.find((i: any) => i.is_canonical);
    const others = (gene.isoforms?.filter((i: any) => !i.is_canonical) || [])
                      .sort((a: any, b: any) => (b.sequence?.length || 0) - (a.sequence?.length || 0))
                      .slice(0, 3);
    return canon ? [canon, ...others] : others.slice(0, 4);
  }, [gene, selectedIsoformIds, sortedIsoformIds]);

  if (!displayIsoforms.length) return null;

  return (
     <div className="flex flex-col gap-6 w-full">
        {/* NMD & Pfam Domain Analysis block */}
        <div className="bg-card/60 backdrop-blur-md border border-primary/20 shadow-xl rounded-xl p-6">
           <div className="flex items-center gap-2 mb-4 pb-2 border-b border-border/50">
              <Activity className="w-4 h-4 text-emerald-400" />
              <p className="text-sm text-muted-foreground">
                 NMD & Pfam Domain Analysis
              </p>
           </div>
           
           <div className="space-y-4">
              {displayIsoforms.map((iso: any) => {
                 const isNMD = iso.nmd_prediction?.nmd_status === "nmd_target" || iso.nmd_prediction?.is_nmd === true;
                 const nmdLabel = iso.nmd_prediction?.nmd_status?.replace("_", " ").toUpperCase() || (isNMD ? "NMD TARGET" : "STABLE");
                 
                 const domains = iso.pfam_domains || [];
                 const total = domains.filter((d: any) => d.domain_status !== 'lost').length;
                 const lost = domains.filter((d: any) => d.domain_status === 'lost').length;
                 const gained = domains.filter((d: any) => d.domain_status === 'gained').length;
                 const truncated = domains.filter((d: any) => d.domain_status === 'truncated').length;

                 return (
                    <div key={iso.id} className="flex flex-col gap-3 p-3 rounded-lg bg-muted/50 border border-border hover:border-primary/30 transition-colors">
                       <div className="flex justify-between items-center">
                          <span className="font-mono text-sm text-primary flex items-center gap-1">
                             <Link className="px-1 w-4 h-4 text-primary/50" />
                             {iso.transcript_name || iso.transcript_id}
                          </span>
                          {iso.is_canonical && <span className="text-[10px] uppercase text-amber-500 font-bold tracking-wider px-2 py-0.5 bg-amber-500/10 border border-amber-500/20 rounded">Canonical Reference</span>}
                       </div>
                       
                       <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
                          {/* NMD Status */}
                          <div className={`flex items-center gap-1.5 text-[11px] px-2 py-1 rounded font-medium ${isNMD ? 'bg-rose-500/10 text-rose-400 border border-rose-500/20' : 'bg-emerald-500/10 text-emerald-400 border border-emerald-500/20'}`}>
                             {isNMD ? <ShieldAlert className="w-3 h-3" /> : <Activity className="w-3 h-3" />}
                             {nmdLabel}
                          </div>

                          {/* Domain Count */}
                          <div className="flex items-center gap-1.5 text-[11px] px-2 py-1 rounded bg-indigo-500/10 text-indigo-400 border border-indigo-500/20 font-medium">
                             <Cpu className="w-3 h-3" />
                             {total} Domains
                          </div>

                          {/* Domain Losses */}
                          {lost > 0 && (
                            <div className="flex items-center gap-1.5 text-[11px] px-2 py-1 rounded bg-rose-500/10 text-rose-400 border border-rose-500/30 font-bold">
                               <MinusCircle className="w-3 h-3" />
                               {lost} Lost
                            </div>
                          )}

                          {/* Domain Gains / Changes */}
                          {gained > 0 && (
                            <div className="flex items-center gap-1.5 text-[11px] px-2 py-1 rounded bg-cyan-500/10 text-cyan-400 border border-cyan-500/30 font-bold">
                               <PlusCircle className="w-3 h-3" />
                               {gained} Gained
                            </div>
                          )}

                          {/* Truncations */}
                          {truncated > 0 && (
                            <div className="flex items-center gap-1.5 text-[11px] px-2 py-1 rounded bg-amber-500/10 text-amber-400 border border-amber-500/30 font-bold col-span-2 md:col-span-1">
                               <AlertTriangle className="w-3 h-3" />
                               {truncated} Truncated
                            </div>
                          )}
                       </div>
                    </div>
                 );
              })}
           </div>
        </div>
     </div>
  );
}
