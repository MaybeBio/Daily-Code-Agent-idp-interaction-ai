"use client";

import { useEffect, useState } from "react";
import { getApiUrl } from "@/lib/api";
import { ChevronDown, ChevronUp, Link as LinkIcon } from "lucide-react";
import { StringNetworkGraph } from "@/components/analysis/InteractionsPanel";

interface Partner {
  partner_gene_symbol: string;
  string_score: number;
}

// Gene-level STRING interactions: the network graph (STRING's own node-edge view,
// proxied same-origin) plus a collapsible, deduped partner list. STRING evidence
// is gene-level — shared across all isoforms — so it lives naturally on the gene page.
export default function GeneStringInteractions({ geneSymbol }: { geneSymbol?: string }) {
  const [partners, setPartners] = useState<Partner[]>([]);
  const [loading, setLoading] = useState(true);
  const [showPartners, setShowPartners] = useState(false);

  useEffect(() => {
    if (!geneSymbol) return;
    setLoading(true);
    fetch(`${getApiUrl()}/api/genes/${encodeURIComponent(geneSymbol)}/string_interactions`)
      .then((r) => (r.ok ? r.json() : { partners: [] }))
      .then((d) => setPartners(d.partners ?? []))
      .catch(() => setPartners([]))
      .finally(() => setLoading(false));
  }, [geneSymbol]);

  if (!geneSymbol) return null;

  return (
    <div className="space-y-4">
      <p className="text-xs text-muted-foreground leading-relaxed">
        STRING functional-interaction network for{" "}
        <strong className="text-foreground">{geneSymbol}</strong>. Gene-level evidence — shared
        across all isoforms; not isoform-specific.
      </p>

      <StringNetworkGraph geneSymbol={geneSymbol} />

      {!loading && partners.length > 0 && (
        <>
          <button
            type="button"
            onClick={() => setShowPartners((v) => !v)}
            className="flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground/90 transition"
          >
            {showPartners ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
            {showPartners ? "Hide" : "Show"} interaction partners ({partners.length})
          </button>

          {showPartners && (
            <div className="rounded-md border border-border divide-y divide-border/40">
              {partners.map((p) => (
                <div key={p.partner_gene_symbol} className="flex items-center gap-3 px-4 py-2 text-sm">
                  <span className="font-semibold text-foreground/90 w-40 shrink-0">
                    {p.partner_gene_symbol}
                  </span>
                  <div className="flex-1 flex items-center gap-2">
                    <div className="w-24 bg-slate-700 h-1.5 rounded-full overflow-hidden">
                      <div
                        className="h-full bg-amber-500"
                        style={{ width: `${(p.string_score || 0) / 10}%` }}
                      />
                    </div>
                    <span className="text-xs font-mono text-foreground/80">{p.string_score}</span>
                  </div>
                  <a
                    href={`https://string-db.org/cgi/network?identifiers=${encodeURIComponent(
                      [geneSymbol, p.partner_gene_symbol].join("\r")
                    )}&species=9606`}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-xs text-primary hover:underline inline-flex items-center gap-1 shrink-0"
                  >
                    <LinkIcon className="w-3 h-3" /> View
                  </a>
                </div>
              ))}
            </div>
          )}
        </>
      )}

      {!loading && partners.length === 0 && (
        <p className="text-xs text-muted-foreground/70 italic">
          No STRING interaction partners stored for this gene.
        </p>
      )}

      <p className="text-[10px] text-muted-foreground/70 leading-snug">
        STRING combined scores draw on seven evidence channels (experiments, databases,
        co-expression, text-mining, phylogenetic profile, gene fusion, neighborhood); a high score
        doesn&apos;t guarantee direct physical binding.
      </p>
    </div>
  );
}
