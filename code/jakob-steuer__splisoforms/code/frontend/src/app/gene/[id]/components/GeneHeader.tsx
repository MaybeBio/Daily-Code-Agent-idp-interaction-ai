"use client";

import { Badge } from "@/components/ui/badge";
import { Database, Dna, ExternalLink, FileText } from "lucide-react";

// Auto-generated stub descriptions from the BC ingestion fallback — hide them
// rather than show "Auto-generated for PacBio isoform PB.x.y" to the user.
const AUTO_DESC_RE = /^Auto-generated for PacBio isoform/i;

function ExternalChip({
  href,
  icon,
  label,
  className,
}: {
  href: string;
  icon: React.ReactNode;
  label: string;
  className: string;
}) {
  return (
    <a
      href={href}
      target="_blank"
      rel="noreferrer noopener"
      className={`flex items-center gap-1.5 px-3 py-1 rounded-full font-mono text-xs transition-colors ${className}`}
    >
      {icon}
      {label}
      <ExternalLink className="w-3 h-3 opacity-50" />
    </a>
  );
}

export default function GeneHeader({ gene }: { gene: any }) {
  const isoformCount = gene.isoforms?.length || 0;

  const description: string | undefined =
    gene.function_summary ||
    (gene.description && !AUTO_DESC_RE.test(gene.description) ? gene.description : undefined);
  const proteinName: string | undefined = gene.protein_name || undefined;

  return (
    <div className="bg-card/60 backdrop-blur-md border border-primary/30 shadow-2xl rounded-2xl p-6 md:p-8 flex flex-col gap-4">
      {/* Top row: title + external links */}
      <div className="flex flex-col md:flex-row md:items-start md:justify-between gap-4">
        <div className="space-y-1 min-w-0">
          <div className="flex items-center gap-3 flex-wrap">
            <h1 className="text-4xl md:text-5xl font-black tracking-tight text-foreground drop-shadow-md">
              {gene.symbol}
            </h1>
            {gene.hgnc_symbol && gene.hgnc_symbol !== gene.symbol && (
              <Badge variant="secondary" className="text-sm px-3 py-1 font-mono">
                HGNC: {gene.hgnc_symbol}
              </Badge>
            )}
          </div>
          {proteinName && (
            <p className="text-sm md:text-base text-foreground/70 mt-1 leading-snug">
              {proteinName}
            </p>
          )}
        </div>

        <div className="flex flex-wrap items-center gap-2 md:justify-end">
          {gene.uniprot_id && (
            <ExternalChip
              href={`https://www.uniprot.org/uniprotkb/${gene.uniprot_id}/entry`}
              icon={<Database className="w-3.5 h-3.5" />}
              label={gene.uniprot_id}
              className="bg-primary/10 border border-primary/20 text-primary hover:bg-primary/20"
            />
          )}
          {gene.ensembl_gene_id && (
            <ExternalChip
              href={`https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${gene.ensembl_gene_id}`}
              icon={<Dna className="w-3.5 h-3.5" />}
              label={gene.ensembl_gene_id}
              className="bg-muted text-foreground/80 hover:bg-muted/80 border border-border/40"
            />
          )}
          {gene.uniprot_id && (
            <ExternalChip
              href={`https://alphafold.ebi.ac.uk/entry/${gene.uniprot_id}`}
              icon={<Database className="w-3.5 h-3.5" />}
              label="AFDB"
              className="bg-blue-500/10 text-blue-400 border border-blue-500/20 hover:bg-blue-500/20"
            />
          )}
        </div>
      </div>

      {/* Description / function summary */}
      {description && (
        <p className="text-sm text-foreground/80 leading-relaxed max-w-4xl">
          {description}
        </p>
      )}

      {/* Bottom row: chromosome, isoform count, AFDB badges, connectivity */}
      <div className="flex flex-wrap items-center gap-3 text-sm text-foreground/80 pt-2 border-t border-border/40">
        <div className="flex items-center gap-1.5 bg-muted px-3 py-1 rounded-full font-mono text-xs">
          <Dna className="w-3.5 h-3.5" />
          {gene.chromosome}:{gene.genomic_start?.toLocaleString("en-US")}-
          {gene.genomic_end?.toLocaleString("en-US")} ({gene.strand})
        </div>
        {gene.canonical_length != null && (
          <div className="flex items-center gap-1.5 bg-muted px-3 py-1 rounded-full font-mono text-xs">
            <FileText className="w-3.5 h-3.5" />
            {gene.canonical_length} aa canonical
          </div>
        )}
        <div className="flex items-center gap-1.5 bg-emerald-500/10 text-emerald-400 border border-emerald-500/20 px-3 py-1 rounded-full font-medium text-xs">
          <FileText className="w-3.5 h-3.5" />
          {isoformCount} Isoform{isoformCount !== 1 ? "s" : ""} Mapped
        </div>
        {gene.afdb_model_version && (
          <Badge variant="outline" className="border-blue-500/30 text-blue-400 text-xs">
            AFDB {gene.afdb_model_version}
          </Badge>
        )}
        {gene.afdb_monomer_cached && (
          <Badge variant="outline" className="border-indigo-500/30 text-indigo-400 text-xs">
            Canonical Monomer Cached
          </Badge>
        )}
      </div>
    </div>
  );
}
