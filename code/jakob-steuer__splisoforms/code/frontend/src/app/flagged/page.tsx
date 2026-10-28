"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { 
  Box, 
  Flag, 
  ShieldAlert, 
  Database, 
  ArrowRight,
  ExternalLink,
  ChevronRight,
  Filter,
  Search,
  Zap
} from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { getApiUrl } from "@/lib/api";

interface FlaggedIsoform {
  id: number;
  gene_id: number;
  transcript_id: string;
  transcript_name?: string | null;
  sequence?: string;
  delta_length?: number;
  gene?: {
    symbol?: string;
    uniprot_id: string;
  };
  metadata_record?: {
    orf_type?: string;
    subcategory?: string;
    uniprot_id?: string;
  };
}

export default function FlaggedPage() {
  const [flagged, setFlagged] = useState<FlaggedIsoform[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [search, setSearch] = useState("");

  useEffect(() => {
    const fetchFlagged = async () => {
      setLoading(true);
      const apiUrl = getApiUrl();
      try {
        const res = await fetch(`${apiUrl}/api/isoforms/flagged`);
        if (!res.ok) throw new Error("Failed to fetch flagged isoforms");
        const data = await res.json();
        setFlagged(data);
      } catch (err: any) {
        setError(err.message || "An error occurred fetching curated isoforms.");
      } finally {
        setLoading(false);
      }
    };
    fetchFlagged();
  }, []);

  const filtered = flagged.filter(iso => 
    iso.gene?.symbol?.toLowerCase().includes(search.toLowerCase()) ||
    iso.transcript_id.toLowerCase().includes(search.toLowerCase()) ||
    iso.transcript_name?.toLowerCase().includes(search.toLowerCase()) ||
    iso.metadata_record?.subcategory?.toLowerCase().includes(search.toLowerCase())
  );

  return (
    <div className="min-h-screen bg-background text-foreground pb-20">
      <header className="border-b border-primary/20 bg-card/50 backdrop-blur sticky top-0 z-50">
        <div className="w-full max-w-[1400px] mx-auto px-6 py-6">
           <div className="flex items-center justify-between gap-6 flex-wrap">
              <div className="flex items-center gap-4">
                 <div className="w-12 h-12 bg-primary/10 rounded-xl flex items-center justify-center border border-primary/20">
                    <Flag className="w-6 h-6 text-primary fill-primary/20" />
                 </div>
                 <div>
                    <h1 className="text-2xl font-bold tracking-tight bg-gradient-to-r from-primary to-secondary bg-clip-text text-transparent">
                       Curated Biological Interest
                    </h1>
                    <p className="text-muted-foreground text-sm flex items-center gap-2">
                       <Database className="w-3.5 h-3.5" /> 
                       Global Shelf of Flagged Isoforms
                    </p>
                 </div>
              </div>

              <div className="flex-1 max-w-md relative">
                 <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
                 <Input 
                    placeholder="Search curated items..." 
                    className="bg-card/50 border-border pl-10 focus:border-primary/50"
                    value={search}
                    onChange={(e) => setSearch(e.target.value)}
                 />
              </div>
           </div>
        </div>
      </header>

      <main className="w-full max-w-[1400px] mx-auto px-6 py-10">
        
        {loading ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 animate-pulse">
            {[1, 2, 3, 4, 5, 6].map(i => (
              <div key={i} className="h-48 bg-card/40 rounded-xl border border-border/60" />
            ))}
          </div>
        ) : filtered.length === 0 ? (
          <div className="flex flex-col items-center justify-center min-h-[40vh] text-center p-12 bg-muted/20 rounded-2xl border border-border/60 border-dashed">
            <Zap className="w-12 h-12 text-muted-foreground/30 mb-4" />
            <h3 className="text-xl font-semibold text-foreground">No curators have flagged items matching your criteria.</h3>
            <p className="text-muted-foreground max-w-sm mt-2">Curated isoforms appear here once they are flagged from their respective Gene dashboards.</p>
          </div>
        ) : (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {filtered.map((iso) => (
              <Card key={iso.id} className="bg-card/40 border-border/60 backdrop-blur group hover:border-primary/50 transition-all duration-300 shadow-xl overflow-hidden flex flex-col">
                <CardHeader className="pb-4">
                   <div className="flex justify-between items-start">
                      <div>
                        <CardTitle className="text-xl font-bold flex items-center gap-2 text-primary group-hover:text-foreground transition-colors">
                           {iso.gene?.symbol}
                        </CardTitle>
                        <CardDescription className="font-mono text-xs text-muted-foreground mt-0.5">
                           {iso.transcript_name || iso.transcript_id}
                        </CardDescription>
                      </div>
                      <Badge variant="outline" className="bg-background/50 text-[10px] text-muted-foreground capitalize">
                         {iso.metadata_record?.subcategory?.replace(/_/g, " ") || "Isoform"}
                      </Badge>
                   </div>
                </CardHeader>
                
                <CardContent className="flex-1 flex flex-col">
                   <div className="space-y-4 flex-1">
                      <div className="flex gap-6">
                         <div className="space-y-1">
                            <span className="text-[10px] font-bold uppercase text-muted-foreground block tracking-widest">Length</span>
                            <span className="text-sm font-medium text-foreground">{iso.sequence?.length || 0} AA</span>
                         </div>
                         <div className="space-y-1">
                            <span className="text-[10px] font-bold uppercase text-muted-foreground block tracking-widest">Delta</span>
                            <Badge variant={iso.delta_length && iso.delta_length > 0 ? "default" : "secondary"} className="h-5 text-[10px]">
                               {iso.delta_length && iso.delta_length > 0 ? `+${iso.delta_length}` : iso.delta_length}
                            </Badge>
                         </div>
                      </div>

                      <div className="pt-2">
                         <div className="text-[10px] font-bold uppercase text-muted-foreground mb-2 tracking-widest">Key Structural Feature</div>
                         <div className="h-1.5 w-full bg-muted rounded-full overflow-hidden">
                            <div className="h-full bg-primary/40 w-2/3" />
                         </div>
                      </div>
                   </div>

                   <div className="mt-8 flex gap-2">
                      <Link href={`/isoform/${iso.id}`} className="flex-1">
                        <Button className="w-full h-9 bg-primary/10 hover:bg-primary text-primary hover:text-foreground border-primary/20 transition-all gap-2 text-xs">
                           Analyze Structure <ChevronRight className="w-4 h-4" />
                        </Button>
                      </Link>
                      <Link href={`/gene/${iso.gene_id}`} className="shrink-0">
                        <Button variant="outline" size="icon" title="Go to Gene Dashboard" className="h-9 w-9 border-border bg-card hover:bg-muted">
                           <ArrowRight className="w-4 h-4" />
                        </Button>
                      </Link>
                   </div>
                </CardContent>
              </Card>
            ))}
          </div>
        )}
      </main>
    </div>
  );
}
