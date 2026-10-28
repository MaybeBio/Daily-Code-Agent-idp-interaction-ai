"use client";

import { useState, useEffect, useRef } from "react";
import { useRouter } from "next/navigation";
import { Search, Pill } from "lucide-react";
import {
  Command,
  CommandEmpty,
  CommandGroup,
  CommandItem,
  CommandList,
  CommandSeparator,
} from "@/components/ui/command";
import { Input } from "@/components/ui/input";
import { useDebounce } from "@/hooks/use-debounce";
import { getApiUrl } from "@/lib/api";

interface SearchResultItem {
  match_type: "gene" | "transcript";
  id: number;              // for transcripts this is the GENE id
  isoform_id?: number | null;   // numeric isoform id → #isoforms=<id> selection
  gene_symbol?: string | null;
  uniprot_id?: string | null;
  transcript_id?: string | null;
  transcript_name?: string | null;
  source_name?: string | null;
  description?: string | null;
}

/** A target–drug pair. `id` is the GENE id, so selecting it opens /drugs/{id}. */
interface DrugSearchItem {
  match_type: "drug";
  id: number;
  gene_symbol?: string | null;
  drug_name: string;
  pubchem_cid?: string | null;
  moa?: string | null;
  n_complexes: number;
}

interface SearchResponse {
  genes: SearchResultItem[];
  transcripts: SearchResultItem[];
  drugs?: DrugSearchItem[];
}

export function SearchCatalog({ compact = false }: { compact?: boolean }) {
  const [open, setOpen] = useState(false);
  const [query, setQuery] = useState("");
  const debouncedQuery = useDebounce(query, 300);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<SearchResponse>({ genes: [], transcripts: [], drugs: [] });
  const [sources, setSources] = useState<{id: number, name: string}[]>([]);
  const [selectedSource, setSelectedSource] = useState<string>("all");
  const router = useRouter();
  const wrapperRef = useRef<HTMLDivElement>(null);

  // Close the dropdown on outside click + Escape.
  useEffect(() => {
    if (!open) return;
    const handlePointer = (e: MouseEvent | TouchEvent) => {
      if (wrapperRef.current && !wrapperRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    const handleKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") setOpen(false);
    };
    document.addEventListener("mousedown", handlePointer);
    document.addEventListener("touchstart", handlePointer);
    document.addEventListener("keydown", handleKey);
    return () => {
      document.removeEventListener("mousedown", handlePointer);
      document.removeEventListener("touchstart", handlePointer);
      document.removeEventListener("keydown", handleKey);
    };
  }, [open]);

  useEffect(() => {
    const fetchSources = async () => {
      try {
        const apiUrl = getApiUrl();
        const res = await fetch(`${apiUrl}/api/sources`);
        if (res.ok) setSources(await res.json());
      } catch (e) {
        console.error("Failed to fetch sources", e);
      }
    };
    fetchSources();
  }, []);

  useEffect(() => {
    if (debouncedQuery.length < 2) {
      setResults({ genes: [], transcripts: [], drugs: [] });
      return;
    }

    const fetchGenes = async () => {
      setLoading(true);
      try {
        const apiUrl = getApiUrl();
        let url = `${apiUrl}/api/genes/search?query=${encodeURIComponent(debouncedQuery)}`;
        if (selectedSource !== "all") {
          url += `&source_id=${selectedSource}`;
        }
        const res = await fetch(url);
        if (res.ok) {
          const data = await res.json();
          setResults(data);
          setOpen(true);
        } else {
          setResults({ genes: [], transcripts: [], drugs: [] });
        }
      } catch (e) {
        console.error("Failed to search", e);
      } finally {
        setLoading(false);
      }
    };

    fetchGenes();
  }, [debouncedQuery, selectedSource]);

  const handleSelectGene = (geneId: number) => {
    setOpen(false);
    setQuery("");
    router.push(`/gene/${geneId}`);
  };

  const handleSelectTranscript = (geneId: number, isoformId?: number | null) => {
    setOpen(false);
    setQuery("");
    // The gene dashboard selects isoforms via the URL hash (#isoforms=<numeric id>),
    // not a ?highlight query param. Fall back to the plain gene page if we lack the id.
    router.push(isoformId != null ? `/gene/${geneId}#isoforms=${isoformId}` : `/gene/${geneId}`);
  };

  // A drug hit is a target–drug PAIR, so it resolves to that target's drug-pocket page.
  const handleSelectDrug = (geneId: number) => {
    setOpen(false);
    setQuery("");
    router.push(`/drugs/${geneId}`);
  };

  const handleKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === "Enter" && query.trim() !== "") {
      router.push(`/search?query=${encodeURIComponent(query)}`);
      setOpen(false);
    }
  };

  return (
    <div ref={wrapperRef} className="w-full relative flex-1">
      <div className="relative w-full text-left" onClick={() => setOpen(true)}>
        <Search className={`absolute left-4 top-1/2 -translate-y-1/2 text-muted-foreground pointer-events-none ${compact ? 'w-4 h-4' : 'w-5 h-5'}`} />
            <Input
              type="text"
              placeholder="Search by gene, drug, UniProt ID, or transcript..."
              className={`w-full bg-background/80 border-primary/40 focus-visible:ring-primary backdrop-blur-sm shadow-xl ${
                 compact ? 'h-10 text-xs pl-10' : 'h-14 text-base pl-12 pr-32'
              }`}
              value={query}
              onChange={(e) => {
                setQuery(e.target.value);
                setOpen(true);
              }}
              onKeyDown={handleKeyDown}
            />
            {!compact && (
              <div className="absolute right-2 top-1/2 -translate-y-1/2 flex gap-1">
                <select 
                  className="bg-primary/10 border-primary/20 text-[10px] uppercase font-bold text-primary rounded px-2 py-1 focus:outline-none cursor-pointer hover:bg-primary/20 transition-colors"
                  value={selectedSource}
                  onChange={(e) => setSelectedSource(e.target.value)}
                  onClick={(e) => e.stopPropagation()}
                >
                  <option value="all" className="bg-slate-900 text-foreground">All Data</option>
                  {sources.map(s => (
                    <option key={s.id} value={s.id} className="bg-slate-900 text-foreground">{s.name}</option>
                  ))}
                </select>
              </div>
            )}
          </div>
      
      {open && (
        <div className="absolute top-full left-0 w-full mt-2 z-50 bg-background/95 backdrop-blur-md border border-primary/20 rounded-md shadow-lg overflow-hidden">
          <Command>
            <CommandList>
              {loading && <div className="p-4 text-sm text-center text-muted-foreground">Searching...</div>}
              {!loading && query.length >= 2 && results.genes.length === 0 && results.transcripts.length === 0 && (results.drugs?.length ?? 0) === 0 && <CommandEmpty>No results found.</CommandEmpty>}
              
              {!loading && query.length < 2 && (
                <CommandGroup heading="Example targets">
                  {[
                    { gene: "MAOB",   hint: "drugs · homo- & heterodimers" },
                    { gene: "TP53",   hint: "neoepitopes · drugs · network" },
                    { gene: "PLCB4",  hint: "heterodimer interfaces" },
                    { gene: "F2",     hint: "drug-pocket fate (argatroban)" },
                    { gene: "ZNF268", hint: "domain loss (zinc fingers)" },
                  ].map(({ gene, hint }) => (
                    <CommandItem
                      key={gene}
                      value={gene}
                      onSelect={() => {
                        setOpen(false);
                        setQuery("");
                        router.push(`/search?query=${gene}`);
                      }}
                      className="flex justify-between items-center py-3 cursor-pointer gap-3"
                    >
                      <span className="font-medium text-foreground">{gene}</span>
                      <span className="text-xs text-muted-foreground bg-muted px-2 py-1 rounded-full whitespace-nowrap">
                        {hint}
                      </span>
                    </CommandItem>
                  ))}
                </CommandGroup>
              )}

              {!loading && query.length >= 2 && results.genes.length > 0 && (
                <CommandGroup heading="Genes">
                  {results.genes.map((gene) => (
                    <CommandItem
                      key={`gene-${gene.id}`}
                      value={`${gene.gene_symbol} ${gene.uniprot_id}`}
                      onSelect={() => handleSelectGene(gene.id)}
                      className="flex justify-between items-center py-3 cursor-pointer"
                    >
                      <span className="font-medium text-foreground">{gene.gene_symbol}</span>
                      <span className="text-sm text-muted-foreground font-mono">{gene.uniprot_id}</span>
                    </CommandItem>
                  ))}
                </CommandGroup>
              )}

              {!loading && query.length >= 2 && results.transcripts.length > 0 && (
                <CommandGroup heading="Transcripts">
                  {results.transcripts.map((t) => (
                    <CommandItem
                      key={`transcript-${t.transcript_id}`}
                      value={`${t.transcript_name ?? ""} ${t.transcript_id} ${t.gene_symbol}`}
                      onSelect={() => handleSelectTranscript(t.id, t.isoform_id)}
                      className="flex justify-between items-center py-3 cursor-pointer"
                    >
                      <span className="font-medium text-foreground font-mono">{t.transcript_name || t.transcript_id}</span>
                      <span className="text-xs text-muted-foreground bg-muted px-2 py-1 rounded-full">
                        {t.source_name}
                      </span>
                    </CommandItem>
                  ))}
                </CommandGroup>
              )}

              {!loading && query.length >= 2 && (results.drugs?.length ?? 0) > 0 && (
                <CommandGroup heading="Drug – Target Pairs">
                  {results.drugs!.map((d) => (
                    <CommandItem
                      key={`drug-${d.id}-${d.pubchem_cid ?? d.drug_name}`}
                      value={`${d.drug_name} ${d.gene_symbol} ${d.pubchem_cid ?? ""}`}
                      onSelect={() => handleSelectDrug(d.id)}
                      className="flex justify-between items-center gap-3 py-3 cursor-pointer"
                    >
                      <span className="flex items-center gap-2 min-w-0">
                        <Pill className="w-3.5 h-3.5 text-emerald-500 shrink-0" />
                        <span className="font-medium text-foreground truncate">{d.drug_name}</span>
                        <span className="text-muted-foreground shrink-0">·</span>
                        <span className="font-medium text-foreground/80 shrink-0">{d.gene_symbol}</span>
                      </span>
                      <span className="text-xs text-muted-foreground shrink-0">
                        {d.n_complexes} complex{d.n_complexes === 1 ? "" : "es"}
                      </span>
                    </CommandItem>
                  ))}
                </CommandGroup>
              )}
            </CommandList>
          </Command>
        </div>
      )}
    </div>
  );
}
