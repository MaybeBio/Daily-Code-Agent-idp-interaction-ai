"use client";

import { useEffect, useState } from "react";
import { useRouter } from "next/navigation";
import { Search, GitCompare, Dna, Check, ArrowRight, X, Loader2 } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Badge } from "@/components/ui/badge";
import { getApiUrl } from "@/lib/api";

const MAX_COMPARE = 5;

interface GeneHit { id: number; gene_symbol: string; uniprot_id?: string | null }
interface IsoRow { id: number; transcript_id: string; transcript_name?: string | null; is_canonical?: boolean; delta_length?: number | null }

export default function CompareBuilder() {
  const router = useRouter();
  const apiUrl = getApiUrl();

  const [geneQuery, setGeneQuery] = useState("");
  const [geneResults, setGeneResults] = useState<GeneHit[]>([]);
  const [gene, setGene] = useState<{ id: number; symbol: string } | null>(null);
  const [isos, setIsos] = useState<IsoRow[]>([]);
  const [picked, setPicked] = useState<number[]>([]);
  const [loadingIsos, setLoadingIsos] = useState(false);

  // Debounced gene search (until a gene is chosen)
  useEffect(() => {
    if (gene) return;
    const q = geneQuery.trim();
    if (q.length < 2) { setGeneResults([]); return; }
    const t = setTimeout(async () => {
      try {
        const res = await fetch(`${apiUrl}/api/genes/search?query=${encodeURIComponent(q)}`).then((r) => r.ok ? r.json() : { genes: [] });
        setGeneResults(res.genes || []);
      } catch { setGeneResults([]); }
    }, 250);
    return () => clearTimeout(t);
  }, [geneQuery, gene, apiUrl]);

  const selectGene = async (g: GeneHit) => {
    setGene({ id: g.id, symbol: g.gene_symbol });
    setGeneResults([]);
    setLoadingIsos(true);
    try {
      const list = await fetch(`${apiUrl}/api/genes/${g.id}/isoforms`).then((r) => r.ok ? r.json() : []);
      setIsos(Array.isArray(list) ? list : []);
    } catch { setIsos([]); }
    finally { setLoadingIsos(false); }
  };

  const reset = () => { setGene(null); setIsos([]); setPicked([]); setGeneQuery(""); };
  const togglePick = (id: number) =>
    setPicked((p) => p.includes(id) ? p.filter((x) => x !== id) : (p.length < MAX_COMPARE ? [...p, id] : p));

  return (
    <div className="max-w-2xl mx-auto py-10">
      <div className="flex items-center gap-3 mb-2">
        <GitCompare className="w-6 h-6 text-primary" />
        <h2 className="text-2xl font-bold text-foreground">Build a comparison</h2>
      </div>
      <p className="text-muted-foreground mb-8 text-sm">
        Pick a gene, then choose 2–{MAX_COMPARE} of its isoforms to compare side-by-side.
      </p>

      {!gene ? (
        /* ── Gene search ────────────────────────────────────────────── */
        <div className="relative">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
          <Input
            autoFocus
            placeholder="Search a gene (symbol or UniProt ID)…"
            className="pl-10 h-12 bg-card/50 border-border/60"
            value={geneQuery}
            onChange={(e) => setGeneQuery(e.target.value)}
          />
          {geneResults.length > 0 && (
            <div className="absolute z-20 mt-2 w-full rounded-xl border border-border/60 bg-background/95 backdrop-blur-md shadow-xl overflow-hidden">
              {geneResults.map((g) => (
                <button
                  key={g.id}
                  onClick={() => selectGene(g)}
                  className="w-full flex items-center justify-between px-4 py-3 hover:bg-muted/50 text-left transition-colors"
                >
                  <span className="flex items-center gap-2 font-medium text-foreground"><Dna className="w-3.5 h-3.5 text-primary" /> {g.gene_symbol}</span>
                  {g.uniprot_id && <span className="text-xs font-mono text-muted-foreground">{g.uniprot_id}</span>}
                </button>
              ))}
            </div>
          )}
        </div>
      ) : (
        /* ── Isoform picker for the chosen gene ─────────────────────── */
        <div className="space-y-4">
          <div className="flex items-center justify-between">
            <span className="flex items-center gap-2 text-sm">
              <Badge className="bg-primary/15 text-primary border-primary/30">{gene.symbol}</Badge>
              <span className="text-muted-foreground">{picked.length} selected</span>
            </span>
            <Button variant="ghost" size="sm" onClick={reset} className="h-7 gap-1 text-xs text-muted-foreground">
              <X className="w-3 h-3" /> Change gene
            </Button>
          </div>

          {loadingIsos ? (
            <div className="flex items-center gap-2 text-muted-foreground py-10 justify-center"><Loader2 className="w-5 h-5 animate-spin text-primary" /> Loading isoforms…</div>
          ) : (
            <div className="max-h-[50vh] overflow-y-auto rounded-xl border border-border/40 divide-y divide-border/30">
              {isos.map((iso) => {
                const sel = picked.includes(iso.id);
                const capped = !sel && picked.length >= MAX_COMPARE;
                return (
                  <button
                    key={iso.id}
                    disabled={capped}
                    onClick={() => togglePick(iso.id)}
                    className={`w-full flex items-center gap-3 px-4 py-3 text-left transition-colors disabled:opacity-40 disabled:cursor-not-allowed ${sel ? "bg-primary/10" : "hover:bg-muted/40"}`}
                  >
                    <span className={`w-5 h-5 rounded border flex items-center justify-center shrink-0 ${sel ? "bg-primary border-primary text-primary-foreground" : "border-border"}`}>
                      {sel && <Check className="w-3.5 h-3.5" />}
                    </span>
                    <span className="font-mono text-sm text-foreground">{iso.transcript_name || iso.transcript_id}</span>
                    {iso.is_canonical && <Badge variant="outline" className="text-[9px] bg-amber-500/10 text-amber-500 border-amber-500/20">canonical</Badge>}
                    {iso.delta_length != null && iso.delta_length !== 0 && (
                      <span className={`text-xs ml-auto ${iso.delta_length > 0 ? "text-emerald-400" : "text-rose-400"}`}>
                        {iso.delta_length > 0 ? "+" : ""}{iso.delta_length} aa
                      </span>
                    )}
                  </button>
                );
              })}
              {isos.length === 0 && <div className="px-4 py-6 text-sm text-muted-foreground text-center">No isoforms found for this gene.</div>}
            </div>
          )}

          <Button
            disabled={picked.length < 2}
            onClick={() => router.push(`/compare?ids=${picked.join(",")}`)}
            className="w-full h-11 gap-2"
          >
            Compare {picked.length >= 2 ? `${picked.length} isoforms` : "(select at least 2)"} <ArrowRight className="w-4 h-4" />
          </Button>
        </div>
      )}
    </div>
  );
}
