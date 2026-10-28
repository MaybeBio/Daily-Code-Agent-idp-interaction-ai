"use client";

import { Suspense, useEffect, useState } from "react";
import { useSearchParams, useRouter } from "next/navigation";
import Link from "next/link";
import { Search, ArrowLeft, Dna, Layers, BarChart3, Syringe, Box, Pill, FlaskConical, Fingerprint } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Textarea } from "@/components/ui/textarea";
import { SearchCatalog } from "@/components/SearchCatalog";
import { getApiUrl } from "@/lib/api";

interface GeneHit {
  match_type: "gene";
  id: number;
  gene_symbol: string;
  uniprot_id?: string | null;
  description?: string | null;
}
interface TranscriptHit {
  match_type: "transcript";
  id: number; // gene_id
  transcript_id: string;
  transcript_name?: string | null;
  gene_symbol?: string | null;
  source_name?: string | null;
}
/** A target–drug pair; `id` is the GENE id, so the card links to /drugs/{id}. */
interface DrugHit {
  match_type: "drug";
  id: number;
  gene_symbol?: string | null;
  drug_name: string;
  pubchem_cid?: string | null;
  moa?: string | null;
  n_complexes: number;
}
interface SearchResponse {
  genes: GeneHit[];
  transcripts: TranscriptHit[];
  drugs?: DrugHit[];
  total_genes: number;
  total_transcripts: number;
  total_drugs?: number;
}
interface SequenceHit {
  match_type: string;
  isoform_id: number;
  gene_id?: number | null;
  transcript_id: string;
  transcript_name?: string | null;
  gene_symbol?: string | null;
  source_name?: string | null;
  seq_length: number;
  identity?: number | null;
}
interface SequenceSearchResponse {
  query_length: number;
  query_md5: string;
  exact: SequenceHit[];
  message?: string | null;
}
interface ResolveResponse {
  query: string;
  normalized_id?: string | null;
  id_type?: string | null;
  status: string; // in_db | resolved_gene | gene_not_in_db | not_found | not_ensembl | error
  gene_id?: number | null;
  gene_symbol?: string | null;
  resolved_ensg?: string | null;
  isoform_count?: number | null;
  message?: string | null;
  closest_isoform_id?: number | null;
  closest_transcript_id?: string | null;
  closest_transcript_name?: string | null;
  closest_gene_id?: number | null;
  closest_gene_symbol?: string | null;
  closest_identity?: number | null;
  closest_exact?: boolean | null;
}

const ENSEMBL_RE = /^ENS[A-Z]*[GTP]\d{6,}(\.\d+)?$/i;

function SearchInner() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const query = searchParams.get("query") || "";

  const [data, setData] = useState<SearchResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(25);

  // Sequence-search mode (Phase 1: exact protein match by md5).
  const [mode, setMode] = useState<"id" | "seq">("id");
  const [seqInput, setSeqInput] = useState("");
  const [seqResult, setSeqResult] = useState<SequenceSearchResponse | null>(null);
  const [seqLoading, setSeqLoading] = useState(false);
  const [seqError, setSeqError] = useState<string | null>(null);

  // Ensembl-ID resolution: when an identifier search comes back empty and the query
  // looks like an Ensembl ID, ask the backend to resolve it to its gene (Phase 3).
  const [resolve, setResolve] = useState<ResolveResponse | null>(null);
  const [resolveLoading, setResolveLoading] = useState(false);

  const handleSeqSearch = async () => {
    if (seqInput.trim().length < 10) return;
    setSeqLoading(true);
    setSeqError(null);
    setSeqResult(null);
    try {
      const res = await fetch(`${getApiUrl()}/api/search/sequence`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence: seqInput }),
      });
      const body = await res.json();
      if (!res.ok) {
        setSeqError(body?.detail || "Search failed.");
      } else {
        setSeqResult(body);
      }
    } catch {
      setSeqError("Could not reach the search service.");
    } finally {
      setSeqLoading(false);
    }
  };

  useEffect(() => {
    setPage(1);
  }, [query]);

  useEffect(() => {
    if (!query.trim()) {
      setData(null);
      return;
    }
    let cancelled = false;
    setLoading(true);
    const apiUrl = getApiUrl();
    fetch(`${apiUrl}/api/genes/search?query=${encodeURIComponent(query)}&limit=${pageSize}&offset=${(page - 1) * pageSize}`)
      .then((r) => (r.ok ? r.json() : { genes: [], transcripts: [], drugs: [], total_genes: 0, total_transcripts: 0, total_drugs: 0 }))
      .then((d) => { if (!cancelled) setData(d); })
      .catch(() => { if (!cancelled) setData({ genes: [], transcripts: [], drugs: [], total_genes: 0, total_transcripts: 0, total_drugs: 0 }); })
      .finally(() => { if (!cancelled) setLoading(false); });
    return () => { cancelled = true; };
  }, [query, page, pageSize]);

  // When the identifier search returns nothing for an Ensembl-shaped query, try to
  // resolve it to a gene via the backend (cached Ensembl REST).
  useEffect(() => {
    const q = query.trim();
    const empty =
      !!data &&
      (data.total_genes ?? 0) === 0 &&
      (data.total_transcripts ?? 0) === 0 &&
      (data.total_drugs ?? 0) === 0;
    if (loading || !q || !empty || !ENSEMBL_RE.test(q)) {
      setResolve(null);
      return;
    }
    let cancelled = false;
    setResolveLoading(true);
    setResolve(null);
    fetch(`${getApiUrl()}/api/search/resolve?id=${encodeURIComponent(q)}`)
      .then((r) => (r.ok ? r.json() : null))
      .then((d) => { if (!cancelled) setResolve(d); })
      .catch(() => { if (!cancelled) setResolve(null); })
      .finally(() => { if (!cancelled) setResolveLoading(false); });
    return () => { cancelled = true; };
  }, [data, query, loading]);

  const geneHits = data?.genes ?? [];
  const txHits = data?.transcripts ?? [];
  const totalGenes = data?.total_genes ?? 0;
  const totalTx = data?.total_transcripts ?? 0;
  const drugHits = data?.drugs ?? [];
  const totalDrugs = data?.total_drugs ?? 0;
  const total = geneHits.length + txHits.length + drugHits.length;
  const totalPages = Math.max(Math.ceil(totalGenes / pageSize), Math.ceil(totalTx / pageSize));
  const hasMore = page < totalPages;

  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-4xl mx-auto p-6 md:p-10">
        <header className="flex items-center gap-4 mb-8">
          <Button variant="ghost" size="icon" onClick={() => router.push("/")} className="hover:bg-primary/10">
            <ArrowLeft className="w-5 h-5" />
          </Button>
          <div className="flex items-center gap-3">
            <Search className="w-6 h-6 text-primary" />
            <div>
              <h1 className="text-3xl font-bold tracking-tight">Search</h1>
              <p className="text-sm text-muted-foreground mt-0.5">
                {query
                  ? <>Results for <span className="text-primary font-medium">&ldquo;{query}&rdquo;</span></>
                  : "Search a gene, UniProt ID, or transcript."}
              </p>
            </div>
          </div>
        </header>

        {/* Mode toggle: identifier search (default) vs protein-sequence search. */}
        <div className="inline-flex items-center gap-1 mb-6 rounded-lg border border-border/50 bg-muted/20 p-1">
          <button
            onClick={() => setMode("id")}
            className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
              mode === "id" ? "bg-primary/15 text-primary" : "text-muted-foreground hover:text-foreground"
            }`}
          >
            <Fingerprint className="w-3.5 h-3.5" /> Identifier
          </button>
          <button
            onClick={() => setMode("seq")}
            className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-colors ${
              mode === "seq" ? "bg-primary/15 text-primary" : "text-muted-foreground hover:text-foreground"
            }`}
          >
            <FlaskConical className="w-3.5 h-3.5" /> Sequence
          </button>
        </div>

        {mode === "id" && (
          <div className="mb-10">
            <SearchCatalog />
          </div>
        )}

        {/* ── Protein-sequence search (exact match) ── */}
        {mode === "seq" && (
          <div className="space-y-6 mb-10">
            <div className="rounded-2xl border border-border/50 bg-card/40 p-5">
              <label className="text-xs uppercase font-bold tracking-widest text-muted-foreground">
                Protein sequence
              </label>
              <Textarea
                value={seqInput}
                onChange={(e) => setSeqInput(e.target.value)}
                placeholder={">optional_header\nMKVLA... (paste a protein sequence or FASTA)"}
                className="mt-2 font-mono text-xs min-h-[140px] resize-y"
                spellCheck={false}
              />
              <div className="flex flex-col sm:flex-row sm:items-center justify-between mt-3 gap-3">
                <p className="text-xs text-muted-foreground">
                  Exact match by protein identity. Similarity search (nearest hits) is coming soon.
                </p>
                <Button
                  onClick={handleSeqSearch}
                  disabled={seqLoading || seqInput.trim().length < 10}
                  className="shrink-0"
                >
                  {seqLoading ? "Searching…" : "Search by sequence"}
                </Button>
              </div>
              {seqError && <p className="text-xs text-rose-400 mt-2">{seqError}</p>}
            </div>

            {seqResult && (
              <section className="space-y-3">
                <h2 className="text-xs uppercase font-bold tracking-widest text-muted-foreground flex items-center gap-2">
                  <Fingerprint className="w-3.5 h-3.5" /> Exact matches ({seqResult.exact.length})
                </h2>
                {seqResult.exact.length === 0 ? (
                  <div className="flex flex-col items-center justify-center py-16 text-center rounded-2xl border border-dashed border-border/60 bg-muted/20">
                    <Box className="w-10 h-10 text-muted-foreground/40 mb-3" />
                    <p className="text-muted-foreground">{seqResult.message || "No exact match."}</p>
                  </div>
                ) : (
                  seqResult.exact.map((hit) => (
                    <Link
                      key={`seq-${hit.isoform_id}`}
                      href={hit.gene_id != null ? `/gene/${hit.gene_id}#isoforms=${hit.isoform_id}` : "#"}
                      className="group flex items-center justify-between gap-3 rounded-xl border border-border/40 bg-card/40 p-4 hover:bg-card/70 hover:border-primary/40 transition-colors"
                    >
                      <div className="min-w-0 flex items-center gap-2 flex-wrap">
                        <span className="font-mono text-sm text-foreground group-hover:text-primary transition-colors">
                          {hit.transcript_name || hit.transcript_id}
                        </span>
                        {hit.gene_symbol && <span className="text-xs text-muted-foreground">{hit.gene_symbol}</span>}
                        <Badge variant="outline" className="text-[10px] border-emerald-500/40 text-emerald-500">Exact match</Badge>
                      </div>
                      <div className="flex items-center gap-2 shrink-0">
                        <Badge variant="outline" className="text-[10px]">{hit.seq_length} aa</Badge>
                        {hit.source_name && <Badge variant="outline" className="text-[10px]">{hit.source_name}</Badge>}
                      </div>
                    </Link>
                  ))
                )}
                <p className="text-[11px] text-muted-foreground/70 font-mono pt-1">
                  query: {seqResult.query_length} aa · md5 {seqResult.query_md5.slice(0, 12)}…
                </p>
              </section>
            )}
          </div>
        )}

        {mode === "id" && loading && (
          <div className="text-center py-16 text-muted-foreground animate-pulse">Searching…</div>
        )}

        {mode === "id" && !loading && query && total === 0 && (
          <>
            {/* Ensembl-ID resolution: "not in database — did you mean…" */}
            {resolveLoading && (
              <div className="rounded-2xl border border-border/60 bg-muted/20 p-5 mb-4 text-sm text-muted-foreground animate-pulse">
                Looking up <span className="font-mono">{query}</span> in Ensembl…
              </div>
            )}

            {!resolveLoading && resolve?.status === "resolved_gene" && resolve.gene_id != null && (
              <div className="rounded-2xl border border-primary/40 bg-primary/5 p-5 mb-4">
                <p className="text-xs uppercase font-bold tracking-widest text-primary mb-2">Not in the database — did you mean…</p>
                <p className="text-sm text-foreground">
                  <span className="font-mono">{resolve.normalized_id || query}</span> isn&rsquo;t stored directly, but it belongs to{" "}
                  <span className="font-bold">{resolve.gene_symbol}</span>
                  {resolve.resolved_ensg && <span className="text-muted-foreground font-mono text-xs ml-1">({resolve.resolved_ensg})</span>}.
                </p>

                {/* Most-similar isoform, computed live against the current DB. */}
                {resolve.closest_isoform_id != null && resolve.closest_identity != null && (
                  <p className="text-sm text-muted-foreground mt-2">
                    {resolve.closest_exact ? (
                      <>Its protein is <span className="text-emerald-500 font-medium">identical</span> to </>
                    ) : (
                      <>Most similar isoform: </>
                    )}
                    <Link
                      href={`/gene/${resolve.closest_gene_id ?? resolve.gene_id}#isoforms=${resolve.closest_isoform_id}`}
                      className="font-mono text-foreground underline decoration-dotted underline-offset-2 hover:text-primary"
                    >
                      {resolve.closest_transcript_name || resolve.closest_transcript_id}
                    </Link>
                    {resolve.closest_exact ? (
                      <> (100% identity).</>
                    ) : (
                      <> — ~{Math.round(resolve.closest_identity * 100)}% identical.</>
                    )}
                  </p>
                )}

                <Link
                  href={`/gene/${resolve.gene_id}`}
                  className="inline-flex items-center gap-1.5 mt-3 px-3 py-1.5 rounded-lg border border-primary/40 bg-primary/10 hover:bg-primary/20 text-xs font-medium transition-colors"
                >
                  <Layers className="w-3.5 h-3.5 text-primary" />
                  View {resolve.gene_symbol}
                  {resolve.isoform_count != null && <> · {resolve.isoform_count} isoform{resolve.isoform_count === 1 ? "" : "s"}</>}
                </Link>
              </div>
            )}

            {!resolveLoading && resolve?.status === "gene_not_in_db" && (
              <div className="rounded-2xl border border-amber-500/40 bg-amber-500/5 p-5 mb-4">
                <p className="text-sm text-foreground">
                  <span className="font-mono">{resolve.normalized_id || query}</span> maps to gene{" "}
                  <span className="font-bold">{resolve.gene_symbol || resolve.resolved_ensg}</span>, which isn&rsquo;t in this database yet.
                </p>
              </div>
            )}

            <div className="flex flex-col items-center justify-center py-20 text-center rounded-2xl border border-dashed border-border/60 bg-muted/20">
              <Box className="w-10 h-10 text-muted-foreground/40 mb-3" />
              <p className="text-muted-foreground">No direct matches for &ldquo;{query}&rdquo;.</p>
              <p className="text-xs text-muted-foreground/70 mt-1">Try a gene symbol (e.g. ASCC2), a drug name (e.g. imatinib), a UniProt ID, or an ENST/ENSP.</p>
            </div>
          </>
        )}

        {mode === "id" && !loading && geneHits.length > 0 && (
          <section className="space-y-3 mb-10">
            <h2 className="text-xs uppercase font-bold tracking-widest text-muted-foreground flex items-center gap-2">
              <Dna className="w-3.5 h-3.5" /> Genes ({totalGenes})
            </h2>
            {geneHits.map((g) => (
              <div key={`g-${g.id}`} className="rounded-xl border border-border/40 bg-card/40 p-4 flex flex-col sm:flex-row sm:items-center gap-3">
                <div className="min-w-0 flex-1">
                  <div className="flex items-center gap-2">
                    <span className="font-bold text-foreground">{g.gene_symbol}</span>
                    {g.uniprot_id && <Badge variant="outline" className="text-[10px] font-mono">{g.uniprot_id}</Badge>}
                  </div>
                  {g.description && <p className="text-xs text-muted-foreground truncate mt-0.5">{g.description}</p>}
                </div>
                {/* Lens shortcuts */}
                <div className="flex flex-wrap gap-2 shrink-0">
                  <Link href={`/gene/${g.id}`} className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg border border-border/60 bg-muted/30 hover:bg-muted/60 hover:border-primary/40 text-xs font-medium transition-colors">
                    <Layers className="w-3.5 h-3.5 text-primary" /> Structure
                  </Link>
                  <Link href={`/explore/impact?query=${encodeURIComponent(g.gene_symbol)}`} className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg border border-border/60 bg-muted/30 hover:bg-muted/60 hover:border-secondary/40 text-xs font-medium transition-colors">
                    <BarChart3 className="w-3.5 h-3.5 text-secondary" /> Impact
                  </Link>
                  <Link href={`/explore/neoepitopes?gene=${encodeURIComponent(g.gene_symbol)}`} className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg border border-border/60 bg-muted/30 hover:bg-muted/60 hover:border-rose-400/40 text-xs font-medium transition-colors">
                    <Syringe className="w-3.5 h-3.5 text-rose-400" /> Neoepitopes
                  </Link>
                </div>
              </div>
            ))}
          </section>
        )}

        {/* Drug hits are target–drug PAIRS: a compound only means something here in the
            context of the protein it was co-folded against. Each card opens that target's
            drug-pocket page. */}
        {mode === "id" && !loading && drugHits.length > 0 && (
          <section className="space-y-3 mb-10">
            <h2 className="text-xs uppercase font-bold tracking-widest text-muted-foreground flex items-center gap-2">
              <Pill className="w-3.5 h-3.5" /> Drug – Target Pairs ({totalDrugs})
            </h2>
            {drugHits.map((d) => (
              <Link
                key={`d-${d.id}-${d.pubchem_cid ?? d.drug_name}`}
                href={`/drugs/${d.id}`}
                className="group flex items-center justify-between gap-3 rounded-xl border border-border/40 bg-card/40 p-4 hover:bg-card/70 hover:border-emerald-500/40 transition-colors"
              >
                <div className="min-w-0">
                  <div className="flex items-center gap-2 flex-wrap">
                    <Pill className="w-3.5 h-3.5 text-emerald-500 shrink-0" />
                    <span className="font-medium text-foreground group-hover:text-emerald-500 transition-colors">
                      {d.drug_name}
                    </span>
                    <span className="text-muted-foreground">·</span>
                    <span className="font-bold text-foreground">{d.gene_symbol}</span>
                    {d.pubchem_cid && (
                      <Badge variant="outline" className="text-[10px] font-mono">CID {d.pubchem_cid}</Badge>
                    )}
                  </div>
                  {d.moa && <p className="text-xs text-muted-foreground truncate mt-0.5">{d.moa}</p>}
                </div>
                <span className="text-xs text-muted-foreground shrink-0">
                  {d.n_complexes} complex{d.n_complexes === 1 ? "" : "es"}
                </span>
              </Link>
            ))}
          </section>
        )}

        {mode === "id" && !loading && txHits.length > 0 && (
          <section className="space-y-3">
            <h2 className="text-xs uppercase font-bold tracking-widest text-muted-foreground flex items-center gap-2">
              <Box className="w-3.5 h-3.5" /> Transcripts ({totalTx})
            </h2>
            {txHits.map((t) => (
              <Link
                key={`t-${t.transcript_id}`}
                href={`/gene/${t.id}?highlight=${encodeURIComponent(t.transcript_id)}`}
                className="group flex items-center justify-between gap-3 rounded-xl border border-border/40 bg-card/40 p-4 hover:bg-card/70 hover:border-primary/40 transition-colors"
              >
                <div className="min-w-0">
                  <span className="font-mono text-sm text-foreground group-hover:text-primary transition-colors">{t.transcript_name || t.transcript_id}</span>
                  {t.gene_symbol && <span className="text-xs text-muted-foreground ml-2">{t.gene_symbol}</span>}
                </div>
                {t.source_name && <Badge variant="outline" className="text-[10px] shrink-0">{t.source_name}</Badge>}
              </Link>
            ))}
          </section>
        )}

        {mode === "id" && !loading && (totalGenes > 0 || totalTx > 0) && (
          <div className="flex flex-col sm:flex-row items-center justify-between gap-4 mt-8 pt-8 border-t border-border/40 pb-8">
            <div className="flex items-center gap-2 text-sm text-muted-foreground">
              <span>Show</span>
              <select
                value={pageSize}
                onChange={(e) => { setPageSize(Number(e.target.value)); setPage(1); }}
                className="bg-transparent border border-border/60 rounded px-2 py-1 outline-none focus:border-primary/50 text-foreground"
              >
                <option value={25}>25</option>
                <option value={50}>50</option>
                <option value={100}>100</option>
              </select>
              <span>per page</span>
            </div>
            
            <div className="flex items-center gap-4 text-sm">
              <span className="text-muted-foreground">Page {page} of {totalPages || 1}</span>
              <div className="flex items-center gap-2">
                <Button
                  variant="outline"
                  size="sm"
                  disabled={page === 1}
                  onClick={() => setPage(p => Math.max(1, p - 1))}
                >
                  Previous
                </Button>
                <Button
                  variant="outline"
                  size="sm"
                  disabled={!hasMore}
                  onClick={() => setPage(p => p + 1)}
                >
                  Next
                </Button>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

export default function SearchPage() {
  return (
    <Suspense fallback={<div className="text-center py-24 animate-pulse text-primary tracking-widest uppercase text-xs">Loading search…</div>}>
      <SearchInner />
    </Suspense>
  );
}
