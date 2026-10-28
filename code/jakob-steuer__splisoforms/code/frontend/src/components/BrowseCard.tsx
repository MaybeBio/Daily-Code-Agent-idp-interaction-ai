"use client";

import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Progress } from "@/components/ui/progress";
import { Dna, Activity, Microscope, Zap, ShieldAlert, Binary, Layers } from "lucide-react";
import { useRouter } from "next/navigation";

interface BrowseIsoformSummary {
  id: number;
  transcript_id: string;
  gene_id: number;
  gene_symbol?: string;
  uniprot_id?: string;
  plddt_mean?: number;
  delta_length?: number;
  source_name?: string;
  has_neoantigen: boolean;
  num_domains_lost: number;
  num_domains_gained: number;
  isoform_count?: number;
  impact_score: number;
}

export function BrowseCard({ item }: { item: BrowseIsoformSummary }) {
  const router = useRouter();
  
  const plddt = item.plddt_mean || 0;
  let plddtColor = "text-primary";
  let plddtBg = "bg-primary/10";
  if (plddt > 90) {
    plddtColor = "text-blue-400";
    plddtBg = "bg-blue-400/10";
  } else if (plddt < 50) {
    plddtColor = "text-orange-400";
    plddtBg = "bg-orange-400/10";
  }

  return (
    <Card 
      className="bg-card/40 border-border/40 hover:border-primary/60 cursor-pointer transition-all hover:-translate-y-1 group relative overflow-hidden"
      onClick={() => router.push(`/isoform/${item.id}`)}
    >

      <CardHeader className="pb-2 space-y-1">
        <div className="flex justify-between items-start mb-1">
          <div className="flex gap-2 items-center">
            <Badge variant="outline" className="bg-muted/50 text-[10px] h-5 border-border/50 px-1.5 flex gap-1 items-center">
              {item.source_name || "Unknown"}
            </Badge>
            {item.isoform_count && item.isoform_count > 1 && (
              <Badge variant="outline" className="text-[10px] h-5 border-primary/20 bg-primary/5 text-primary/70 flex gap-1 items-center px-1.5">
                <Layers className="w-2.5 h-2.5" />
                {item.isoform_count} Isoforms
              </Badge>
            )}
          </div>
          <div className="flex gap-1">
             {item.has_neoantigen && (
               <Badge className="bg-magenta-500/20 text-magenta-400 border-magenta-500/30 text-[9px] h-5 px-1">
                 Neoantigen
               </Badge>
             )}
             {item.num_domains_lost > 0 && (
               <Badge className="bg-orange-500/20 text-orange-400 border-orange-500/30 text-[9px] h-5 px-1">
                 -{item.num_domains_lost} Domains
               </Badge>
             )}
          </div>
        </div>
        
        <CardTitle className="text-md flex items-center gap-2 text-foreground group-hover:text-primary transition-colors">
          <Activity className="w-3.5 h-3.5 text-primary" />
          <span className="truncate max-w-[120px]">{item.gene_symbol || item.transcript_id}</span>
          <span className="text-[10px] font-mono text-muted-foreground ml-auto">{item.uniprot_id}</span>
        </CardTitle>
        <CardDescription className="text-[10px] font-mono truncate">
          {item.transcript_id}
        </CardDescription>
      </CardHeader>

      <CardContent className="pt-2">
        <div className="grid grid-cols-2 gap-4 mb-4">
          <div className={`${plddtBg} rounded-md p-2 border border-border/50`}>
            <div className="text-[9px] text-muted-foreground uppercase tracking-tight mb-0.5">pLDDT</div>
            <div className={`text-lg font-bold ${plddtColor}`}>{plddt.toFixed(1)}</div>
          </div>
          <div className="bg-muted/10 rounded-md p-2 border border-border/50">
            <div className="text-[9px] text-muted-foreground uppercase tracking-tight mb-0.5">Seq Δ</div>
            <div className={`text-lg font-bold ${item.delta_length && item.delta_length < 0 ? 'text-destructive' : 'text-primary'}`}>
              {item.delta_length && item.delta_length > 0 ? '+' : ''}{item.delta_length || 0}
            </div>
          </div>
        </div>

        <div className="space-y-2">
           <div className="flex justify-between items-center text-[10px]">
             <span className="text-muted-foreground flex items-center gap-1">
               <Zap className="w-3 h-3 text-amber-500" /> Impact Score
             </span>
             <span className="font-bold text-amber-500">{item.impact_score?.toFixed(1) ?? "0.0"}<span className="text-muted-foreground font-normal">/100</span></span>
           </div>
           {/* Structural Impact Score (SIS v2) is a normalised 0–100 composite */}
           <Progress value={Math.min(100, Math.max(0, item.impact_score ?? 0))} className="h-1 bg-muted" />
        </div>
      </CardContent>
    </Card>
  );
}
