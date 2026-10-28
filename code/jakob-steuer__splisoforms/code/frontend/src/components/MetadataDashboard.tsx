"use client";

import { Badge } from "@/components/ui/badge";
import { Database, ExternalLink, Info } from "lucide-react";

export interface MetadataRecord {
  metadata_kind?: "sqanti3" | "gencode" | "afdb_reference" | "custom" | null;
  chrom?: string;
  strand?: string;
  length?: number;
  biotype?: string;
  orf_type?: string;
  protein_length_aa?: number;
  delta_length_aa?: number;
  sample_source?: string;
  associated_gene?: string;
  canonical_gene_name?: string;
  canonical_protein_length?: number;
  folded_orf_mismatch?: boolean;
  folded_best_ref?: string;
  folded_best_ref_identity?: number;
  structural_category?: string;
  subcategory?: string;
  all_canonical?: string;
  appris_tag?: string;
  tsl?: number;
  gencode_tags?: string[];
  uniprot_id?: string;
  ensp_id?: string;
  enst_id?: string;
  refseq_id?: string;
  hgnc_id?: string;
  is_mane_select?: boolean;
  pdb_ids?: string[];
}

export function MetadataDashboard({ metadata, fallbackId, isoform }: { metadata?: MetadataRecord | null, fallbackId?: string, isoform?: any }) {
  // Protein name + function come from the DB (genes.protein_name / genes.function_summary,
  // backfilled from UniProt), NOT a live UniProt call. The panel used to fetch
  // /api/external/uniprot/{metadata.uniprot_id} on every page view. That was wrong twice
  // over:
  //   * Correctness — metadata.uniprot_id is ISOFORM-level and diverges from the gene's
  //     accession on 12,806 isoforms, where it is typically a UniProt *best-hit homolog*.
  //     A card titled "Canonical Biological Function" was therefore describing a different
  //     protein than the gene the isoform belongs to. (Same flawed best-hit signal that
  //     caused the 2026-07 gene mis-attribution.)
  //   * Reliability — an uncached outbound call to rest.uniprot.org on every isoform view
  //     of a public site, 5s timeout, degrading on any UniProt outage or rate-limit.
  // Both fields are already returned by /api/isoform/{id} under `gene`, so this costs
  // nothing extra. Refresh them by re-running backfill_gene_descriptions.py.
  const gene = isoform?.gene;

  // Basic Identity. transcript_name is the readable handle (MET-201 /
  // MET-206-NOV-a6d00a); transcript_id remains the stable source identifier. Novel isoforms
  // with no close GENCODE match have no name.
  const tId = isoform?.transcript_id || "Unknown";
  const tName: string | undefined = isoform?.transcript_name || undefined;
  const geneSymbol = gene?.symbol || metadata?.associated_gene || "N/A";
  const proteinName = gene?.protein_name || "Unknown Protein";
  const biotype = metadata?.biotype?.replace(/_/g, " ");
  
  // Exons & Genomic Coordinates
  let exonCount: number | null = null;
  let genomicCoords: string | null = null;
  
  if (isoform?.exons_json) {
    try {
      const exons = JSON.parse(isoform.exons_json);
      exonCount = exons.length;
      if (exonCount !== null && exonCount > 0) {
        const starts = exons.map((e: any) => e[0]);
        const ends = exons.map((e: any) => e[1]);
        const min = Math.min(...starts, ...ends);
        const max = Math.max(...starts, ...ends);
        genomicCoords = `${min.toLocaleString()}-${max.toLocaleString()}`;
      }
    } catch(e) {}
  }
  
  // Locus
  const locus = (metadata?.chrom && metadata?.strand && metadata.chrom !== "unknown")
    ? `${metadata.chrom}${genomicCoords ? ':' + genomicCoords : ''} (${metadata.strand})`
    : null;

  // SQANTI3 classification
  const sqanti3Categories = new Set([
    "FSM", "ISM", "NIC", "NNC",
    "full-splice_match", "incomplete-splice_match",
    "novel_in_catalog", "novel_not_in_catalog",
    "antisense", "fusion", "genic", "genic_intron", "intergenic",
  ]);
  const inferredKind = metadata?.metadata_kind ??
    (metadata?.structural_category && sqanti3Categories.has(metadata.structural_category) ? "sqanti3" : null);

  // Delta Bar Setup
  const pLen = metadata?.protein_length_aa ?? isoform?.sequence?.length ?? 0;
  const delta: number | null = typeof isoform?.delta_length === "number" ? isoform.delta_length : null;
  const cLen = delta !== null ? pLen - delta : (isoform?.gene?.canonical_length ?? pLen);

  const maxLen = Math.max(pLen, cLen, 100);
  const pWidth = `${Math.min(100, (pLen / (maxLen || 1)) * 100)}%`;
  const cWidth = `${Math.min(100, (cLen / (maxLen || 1)) * 100)}%`;

  const isTruncated = delta !== null && delta < -0.01;
  const isMatch     = delta !== null && Math.abs(delta) < 0.01;
  const isTrueCanonical = isoform?.is_canonical === true;
  const isUnannotated = !metadata && isoform?.transcript_id && !isoform.transcript_id.startsWith("PB.") && !isoform.transcript_id.startsWith("ENST");

  // Function summary. Some accessions (e.g. Q92802 / N4BP2L2) genuinely carry no FUNCTION
  // or SIMILARITY comment in UniProt at all, so there is nothing to show — render an
  // explicit placeholder rather than leaving half the Overview grid blank.
  const funcText: string | undefined = gene?.function_summary || undefined;
  const geneProteinName: string | undefined = gene?.protein_name || undefined;
  const geneUniprotId: string | undefined =
    gene?.uniprot_id || metadata?.uniprot_id?.replace(/^AFDB_/, "") || fallbackId || undefined;
  const uniprotUrl = geneUniprotId
    ? `https://www.uniprot.org/uniprotkb/${geneUniprotId}/entry`
    : undefined;

  return (
    <div className="bg-card/50 border border-border rounded-xl p-6 min-h-full flex flex-col gap-6">
      <div className="flex justify-start">
        {isTrueCanonical ? (
          <Badge className="bg-secondary/10 text-secondary border-secondary/20 shadow-sm shadow-secondary/10">
            Canonical Reference
          </Badge>
        ) : inferredKind === "sqanti3" ? (
          <Badge className="bg-primary/10 text-primary border-primary/20 shadow-sm">
            SQANTI3 Annotated
          </Badge>
        ) : isUnannotated ? (
          <Badge className="bg-muted text-foreground/90 border-none shadow-sm">
            Unannotated
          </Badge>
        ) : (
          <Badge className="bg-amber-500/10 text-amber-400 border-amber-400/20 shadow-sm">
            Alternative Isoform
          </Badge>
        )}
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 items-stretch">
        {/* Row 1, Col 1: Identity & Genomics */}
        <div className="bg-muted/30 p-4 rounded-lg border border-border/40 flex flex-col gap-4">
          <div className="grid grid-cols-2 gap-4">
            <div className="flex flex-col gap-1 overflow-hidden">
              <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Transcript Name</span>
              <span className="text-sm font-medium text-foreground mt-1 truncate font-mono" title={tName || "—"}>
                {tName || <span className="text-muted-foreground/60 font-sans">—</span>}
              </span>
            </div>
            <div className="flex flex-col gap-1 overflow-hidden">
              <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Transcript ID</span>
              <span className="text-sm font-medium text-foreground mt-1 truncate font-mono" title={tId}>{tId}</span>
            </div>
            <div className="flex flex-col gap-1">
              <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Gene Symbol</span>
              <span className="text-sm font-medium text-foreground mt-1">{geneSymbol}</span>
            </div>
            <div className="flex flex-col gap-1">
              <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">UniProt</span>
              <span className="text-sm font-medium text-foreground mt-1">
                {geneUniprotId ? (
                  <a
                    href={`https://www.uniprot.org/uniprotkb/${geneUniprotId}/entry`}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="font-mono text-primary/90 hover:text-primary hover:underline inline-flex items-center gap-1"
                  >
                    {geneUniprotId}
                    <ExternalLink className="w-2.5 h-2.5 opacity-70" />
                  </a>
                ) : (
                  <span className="text-muted-foreground/60">—</span>
                )}
              </span>
            </div>
            <div className="flex flex-col gap-1 col-span-2">
              <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Protein Name</span>
              <span className="text-sm font-medium text-foreground mt-1">{proteinName}</span>
            </div>
          </div>
          
          <div className="h-px w-full bg-border/50" />
          
          <div className="grid grid-cols-2 gap-4">
            {biotype && (
              <div className="flex flex-col gap-1">
                <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Biotype</span>
                <span className="text-sm font-medium text-foreground mt-1 capitalize">{biotype}</span>
              </div>
            )}
            {exonCount !== null && (
              <div className="flex flex-col gap-1">
                <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Number of Exons</span>
                <span className="text-sm font-medium font-mono text-foreground mt-1">{exonCount}</span>
              </div>
            )}
            {locus && (
              <div className="flex flex-col gap-1 col-span-2">
                <span className="text-[10px] text-muted-foreground uppercase tracking-widest font-semibold">Genomic Locus</span>
                <span className="text-sm font-medium font-mono text-foreground mt-1">{locus}</span>
              </div>
            )}
          </div>
        </div>

        {/* Row 1, Col 2: Biological Function (Canonical).
            The card is absolutely positioned inside its grid cell at lg+, so its (possibly
            long) text does NOT drive the row height — the row is sized by the identity card
            on the left, and this card fills exactly that height with the summary scrolling
            inside it. Below lg the grid is single-column, so it flows normally instead. */}
        <div className="relative">
          {funcText ? (
            <div className="bg-muted/30 p-4 rounded-lg border border-border/40 flex flex-col gap-2 lg:absolute lg:inset-0">
              <span className="text-[10px] uppercase tracking-widest text-muted-foreground font-semibold flex items-center gap-1 shrink-0">
                <Database className="w-3 h-3" /> Canonical Biological Function
              </span>
              <p className="text-sm text-foreground/80 leading-relaxed flex-1 min-h-0 overflow-y-auto pr-2 lg:max-h-none max-h-56">
                {funcText}
              </p>
              {uniprotUrl && (
                <a
                  href={uniprotUrl}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="shrink-0 inline-flex items-center gap-1 text-[10px] text-muted-foreground hover:text-primary transition-colors w-fit"
                  title={`UniProt entry ${geneUniprotId}`}
                >
                  Source: UniProt <span className="font-mono">{geneUniprotId}</span>
                  <ExternalLink className="w-2.5 h-2.5 opacity-70" />
                </a>
              )}
            </div>
          ) : (
            <div className="bg-muted/20 p-4 rounded-lg border border-dashed border-border/50 flex flex-col gap-2 lg:absolute lg:inset-0">
              <span className="text-[10px] uppercase tracking-widest text-muted-foreground font-semibold flex items-center gap-1 shrink-0">
                <Database className="w-3 h-3" /> Canonical Biological Function
              </span>
              <div className="flex-1 min-h-0 overflow-y-auto pr-2">
                {geneProteinName && (
                  <p className="text-sm text-foreground/70 leading-relaxed">{geneProteinName}</p>
                )}
                <p className="text-xs text-muted-foreground/70 leading-relaxed mt-1">
                  No functional description available
                  {geneUniprotId ? <> for <span className="font-mono">{geneUniprotId}</span></> : null}
                  {" "}— UniProt has no function annotation for this protein.
                </p>
              </div>
              {uniprotUrl && (
                <a
                  href={uniprotUrl}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="shrink-0 inline-flex items-center gap-1 text-[10px] text-muted-foreground hover:text-primary transition-colors w-fit"
                  title={`UniProt entry ${geneUniprotId}`}
                >
                  View on UniProt <span className="font-mono">{geneUniprotId}</span>
                  <ExternalLink className="w-2.5 h-2.5 opacity-70" />
                </a>
              )}
            </div>
          )}
        </div>

      </div>

      {/* Full Width Bottom: Sequence Length Comparison */}
      <div className="flex flex-col gap-3 mt-6 pt-4 border-t border-border/50">
        <div className="flex justify-between items-end mb-1">
          <span className="text-sm font-semibold text-foreground">Sequence Length Comparison</span>
          <span className="text-xs font-mono text-foreground/80">
             {isTrueCanonical ? "Match" : delta === null ? "Δ ?" : `Δ ${delta > 0 ? '+' : ''}${delta} aa`}
          </span>
        </div>

        <div className="flex flex-col gap-3 bg-muted/20 p-4 rounded-lg border border-border/50">
          <div className="flex items-center gap-4">
             <div className="w-24 text-xs text-right text-muted-foreground">{isTrueCanonical ? "Sequence" : "Canonical"} ({cLen}aa)</div>
             <div className="flex-1 h-3 bg-muted rounded-full overflow-hidden">
               <div className="h-full bg-slate-500/50 rounded-full transition-all duration-1000" style={{ width: cWidth }} />
             </div>
          </div>
          {!isTrueCanonical && (
            <div className="flex items-center gap-4">
               <div className="w-24 text-xs text-right font-medium text-foreground">Isoform ({pLen}aa)</div>
               <div className="flex-1 h-3 bg-muted rounded-full overflow-hidden flex">
                 <div className={`h-full rounded-full transition-all duration-1000 ${isMatch ? 'bg-secondary/80' : isTruncated ? 'bg-destructive/80' : 'bg-primary/80'}`} style={{ width: pWidth }} />
               </div>
            </div>
          )}
        </div>
      </div>

    </div>
  );
}
