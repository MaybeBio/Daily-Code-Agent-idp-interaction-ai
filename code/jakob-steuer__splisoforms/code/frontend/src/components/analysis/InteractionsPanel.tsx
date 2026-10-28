"use client";

import React, { useMemo, useState } from "react";
import { Badge } from "@/components/ui/badge";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { AlertCircle, ChevronDown, ChevronUp, Link as LinkIcon, Scan, Shield, Zap, HelpCircle } from "lucide-react";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";

function InfoTip({ text }: { text: string }) {
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger>
          <HelpCircle className="w-2.5 h-2.5 text-muted-foreground/50 hover:text-muted-foreground transition-colors cursor-help ml-0.5" />
        </TooltipTrigger>
        <TooltipContent className="bg-popover border-border text-foreground text-xs max-w-xs shadow-xl font-normal normal-case tracking-normal">
          {text}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}
import { Progress } from "@/components/ui/progress";
import { Button, buttonVariants } from "@/components/ui/button";
import { cn } from "@/lib/utils";
import Link from "next/link";
import { getAlignedIndices } from "@/utils/alignment";
import { getApiUrl } from "@/lib/api";

// STRING's own node-edge network, served through our backend proxy (same-origin,
// so it clears the site CSP) and centred on the gene. Hides itself if the image
// can't be loaded (e.g. STRING unavailable or no graph for the gene).
export function StringNetworkGraph({ geneSymbol }: { geneSymbol: string }) {
  const [failed, setFailed] = useState(false);
  if (failed) return null;
  const src = `${getApiUrl()}/api/genes/${encodeURIComponent(geneSymbol)}/string_network`;
  const interactiveUrl = `https://string-db.org/cgi/network?identifiers=${encodeURIComponent(geneSymbol)}&species=9606`;
  return (
    <div className="rounded-md border border-border bg-popover/30 p-3 space-y-2">
      <div className="flex items-center justify-between gap-2 flex-wrap">
        <p className="text-xs text-muted-foreground">
          High-confidence functional neighbourhood (STRING, score ≥ 0.7). Gene-level — not isoform-specific.
        </p>
        <a
          href={interactiveUrl}
          target="_blank"
          rel="noopener noreferrer"
          className="text-xs font-semibold text-primary hover:underline inline-flex items-center gap-1 shrink-0"
        >
          Open interactive on STRING <LinkIcon className="w-3 h-3" />
        </a>
      </div>
      {/* White backing so STRING's dark node labels stay legible in dark mode. */}
      <div className="bg-white rounded overflow-hidden">
        <img
          src={src}
          alt={`STRING interaction network for ${geneSymbol}`}
          onError={() => setFailed(true)}
          className="w-full h-auto max-h-[440px] object-contain"
        />
      </div>
    </div>
  );
}


interface HeterodimerAnalysis {
  id: number;
  partner_uniprot_id?: string;
  partner_name?: string;
  afdb_accession?: string;
  total_interface_residues?: number;
  count_absent?: number;
  count_sasa_changed?: number;
  disruption_score?: number;
  interface_residues_json?: string;
}

interface HomodimerAnalysis {
  id: number;
  afdb_accession?: string;
  total_interface_residues?: number;
  count_absent?: number;
  count_sasa_changed?: number;
  disruption_score?: number;
  interface_residues_json?: string;
}

interface StringInteraction {
  id: number;
  partner_gene_symbol?: string;
  partner_ensp_id?: string;       // renamed from partner_uniprot_id (audit-5)
  string_score?: number;
}

interface InteractionsPanelProps {
  homodimer?: HomodimerAnalysis;
  heterodimers?: HeterodimerAnalysis[];
  stringInteractions?: StringInteraction[];
  isoformId: string | number;
  /** This isoform's gene symbol — used to build the STRING network link for the pair. */
  geneSymbol?: string;
  isCanonical?: boolean;
  /** Canonical protein sequence — used to map interface positions to isoform coords */
  canonicalSequence?: string;
  /** This isoform's sequence — used to map interface positions to isoform coords */
  isoformSequence?: string;
  /** Called with isoform-space residue positions to highlight in the Mol* viewer */
  onSelectResidues?: (positions: number[]) => void;
}


interface ComplexSectionProps {
  complex: HomodimerAnalysis | HeterodimerAnalysis;
  isCanonical: boolean;
  canonicalSequence?: string;
  isoformSequence?: string;
  onSelectResidues?: (positions: number[]) => void;
  title: string;
  type: "homodimer" | "heterodimer";
}

const ComplexSection: React.FC<ComplexSectionProps> = ({
  complex, isCanonical, canonicalSequence, isoformSequence, onSelectResidues, title, type
}) => {
  const [showResidueTable, setShowResidueTable] = useState(false);
  const [showAllResidues, setShowAllResidues] = useState(false);
  const [interfaceActive, setInterfaceActive] = useState(false);

  const interfaceResidues: any[] = useMemo(
    () => (complex.interface_residues_json ? JSON.parse(complex.interface_residues_json) : []),
    [complex.interface_residues_json]
  );

  const disruptionScore = complex.disruption_score ?? 0;
  const totalAffected = (complex.count_absent ?? 0) + (complex.count_sasa_changed ?? 0);
  const isLikelyDisrupted = disruptionScore > 0.5;

  const canonToIsoform = useMemo<Map<number, number> | null>(() => {
    if (isCanonical || !canonicalSequence || !isoformSequence) return null;
    const { idxCanon, idxQuery } = getAlignedIndices(canonicalSequence, isoformSequence);
    const map = new Map<number, number>();
    idxCanon.forEach((c, i) => map.set(c, idxQuery[i]));
    return map;
  }, [isCanonical, canonicalSequence, isoformSequence]);

  const isoformInterfacePositions = useMemo<number[]>(() => {
    return interfaceResidues
      .filter(r => r.status !== "absent")
      .map(r => {
        if (isCanonical || !canonToIsoform) return r.pos as number;
        return canonToIsoform.get(r.pos) ?? null;
      })
      .filter((p): p is number => p !== null);
  }, [interfaceResidues, isCanonical, canonToIsoform]);

  const handleVisualizeInterface = () => {
    if (!onSelectResidues || isoformInterfacePositions.length === 0) return;
    if (interfaceActive) {
      onSelectResidues([]);
      setInterfaceActive(false);
    } else {
      onSelectResidues(isoformInterfacePositions);
      setInterfaceActive(true);
    }
  };

  const barColor =
    isCanonical              ? "bg-emerald-500"
    : disruptionScore >= 0.6 ? "bg-red-500"
    : disruptionScore >= 0.3 ? "bg-amber-500"
    :                          "bg-emerald-500";

  const visibleResidues = showAllResidues ? interfaceResidues : interfaceResidues.slice(0, 10);

  return (
    <div className="space-y-4 pt-4 first:pt-0">
      <div className="flex items-center justify-between">
        <h4 className="text-md font-semibold text-foreground flex items-center gap-2">
          <Shield className="w-4 h-4 text-blue-400" /> {title}
        </h4>
        {isCanonical ? (
          <Badge variant="outline" className="text-emerald-400 border-emerald-500/40 bg-emerald-500/10 flex items-center gap-1">
            <Shield className="w-3 h-3" /> Reference structure
          </Badge>
        ) : isLikelyDisrupted ? (
          <Badge variant="destructive" className="bg-red-500/20 text-red-400 border-red-500/50 flex items-center gap-1">
            <AlertCircle className="w-3 h-3" /> Interaction likely disrupted
          </Badge>
        ) : totalAffected > 0 ? (
          <Badge variant="outline" className="text-amber-400 border-amber-400/40 bg-amber-500/10 flex items-center gap-1">
            <AlertCircle className="w-3 h-3" /> Partial interface loss
          </Badge>
        ) : (
          <Badge variant="outline" className="text-emerald-400 border-emerald-400/40 bg-emerald-500/10 flex items-center gap-1">
            <Shield className="w-3 h-3" /> Interface preserved
          </Badge>
        )}
      </div>

      <div className="bg-muted/30 rounded-lg p-4 border border-border/50 space-y-3">
        {isCanonical ? (
          <p className="text-sm text-foreground/80">
            This is the <span className="font-semibold text-emerald-400">reference {type} structure</span>.
            All <span className="font-mono font-semibold">{complex.total_interface_residues}</span> interface residues are intact. Isoform-specific disruption scores are computed relative to this interface.
          </p>
        ) : (
          <>
            <div className="flex justify-between items-center">
              <span className="text-sm text-foreground/80">
                Interface disruption: <span className="font-mono font-semibold">{totalAffected}/{complex.total_interface_residues}</span> residues affected
              </span>
              <span className={`text-sm font-mono font-bold ${
                disruptionScore >= 0.6 ? "text-red-400"
                : disruptionScore >= 0.3 ? "text-amber-400"
                : "text-emerald-400"
              }`}>
                {(disruptionScore * 100).toFixed(0)}%
              </span>
            </div>
            {type === "heterodimer" && "partner_name" in complex && complex.partner_name && (
              <p className="text-xs text-muted-foreground/80 font-medium">Complex of {complex.partner_name}</p>
            )}
            <div className="h-2 w-full bg-slate-200 dark:bg-slate-700 rounded-full overflow-hidden flex">
              <div className={`h-full transition-all duration-500 ${barColor}`} style={{ width: `${disruptionScore * 100}%` }} />
            </div>
            <div className="flex gap-4 text-xs text-muted-foreground pt-1">
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-red-500/70 inline-block" />{complex.count_absent} absent</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-amber-500/70 inline-block" />{complex.count_sasa_changed} perturbed</span>
              <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-emerald-500/70 inline-block" />{(complex.total_interface_residues ?? 0) - totalAffected} present</span>
            </div>
          </>
        )}
        <p className="text-xs text-muted-foreground italic">
          Complex: <a href={`https://alphafold.ebi.ac.uk/entry/${complex.afdb_accession}`} target="_blank" rel="noopener noreferrer" className="font-mono hover:text-primary underline underline-offset-2">{complex.afdb_accession}</a> · AlphaFold-Multimer (NVIDIA)
          {type === "heterodimer" && "partner_uniprot_id" in complex && (
            <span> · Interacting with: <a href={`https://www.uniprot.org/uniprotkb/${complex.partner_uniprot_id}/entry`} target="_blank" rel="noopener noreferrer" className="font-mono hover:text-primary underline underline-offset-2">{complex.partner_uniprot_id}</a></span>
          )}
        </p>
      </div>

      <div className="flex items-center gap-2">
        {onSelectResidues && isoformInterfacePositions.length > 0 && (
          <button onClick={handleVisualizeInterface} className={cn("flex items-center gap-1.5 text-xs font-medium px-3 py-1.5 rounded-md border transition-colors", interfaceActive ? "border-blue-400/70 bg-blue-500/20 text-blue-300 hover:bg-blue-500/10" : "border-blue-500/40 bg-blue-500/10 text-blue-400 hover:bg-blue-500/20 hover:text-blue-300")}>
            <Scan className="w-3.5 h-3.5" />
            {interfaceActive ? "Deselect interface" : "Visualize interface"}
            {!isCanonical && <span className="opacity-60 font-normal">({isoformInterfacePositions.length} residues)</span>}
          </button>
        )}
        <button onClick={() => { setShowResidueTable(v => !v); setShowAllResidues(false); }} className="flex items-center gap-1.5 text-xs text-muted-foreground hover:text-foreground transition-colors ml-auto">
          {showResidueTable ? <><ChevronUp className="w-3 h-3" /> Hide residue table</> : <><ChevronDown className="w-3 h-3" /> View {interfaceResidues.length} interface residues</>}
        </button>
      </div>

      {showResidueTable && (
        <div className="rounded-md border border-border">
          <Table>
            <TableHeader className="bg-popover/50">
              <TableRow>
                <TableHead className="w-20">Position</TableHead>
                <TableHead>{isCanonical ? "Interface role" : "Status in isoform"}</TableHead>
                <TableHead className="text-right">Can. SASA</TableHead>
                <TableHead className="text-right">Nov. SASA</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {visibleResidues.map((res: any, idx: number) => {
                const isClickable = res.status !== "absent" && onSelectResidues;
                return (
                  <TableRow 
                    key={idx} 
                    className={cn("border-border/50 hover:bg-muted/30 transition-colors", isClickable && "cursor-pointer hover:bg-muted/50")}
                    onClick={() => {
                      if (isClickable) {
                        const isoPos = isCanonical ? res.pos : (canonToIsoform?.get(res.pos) ?? null);
                        if (isoPos !== null) onSelectResidues([isoPos]);
                      }
                    }}
                  >
                    <TableCell className="font-mono text-foreground/90">{res.pos}</TableCell>
                    <TableCell>
                      {isCanonical ? (
                        <Badge variant="outline" className="text-blue-400 border-blue-400/30 bg-blue-400/5">Interface residue</Badge>
                      ) : res.status === "absent" ? (
                        <Badge variant="outline" className="text-red-400 border-red-400/30 bg-red-400/5">Absent</Badge>
                      ) : res.status === "perturbed" ? (
                        <Badge variant="outline" className="text-amber-400 border-amber-400/30 bg-amber-400/5">Perturbed</Badge>
                      ) : (
                        <Badge variant="outline" className="text-emerald-400 border-emerald-400/30 bg-emerald-400/5">Present</Badge>
                      )}
                    </TableCell>
                    <TableCell className="text-right font-mono text-foreground/80">{res.sasa_canonical?.toFixed(1) ?? "--"}</TableCell>
                    <TableCell className="text-right font-mono text-foreground/80">{res.sasa_novel?.toFixed(1) ?? "--"}</TableCell>
                  </TableRow>
                );
              })}
            </TableBody>
          </Table>

          {interfaceResidues.length > 10 && (
            <button onClick={() => setShowAllResidues(v => !v)} className="w-full p-2 flex items-center justify-center gap-1.5 text-xs text-muted-foreground hover:text-foreground hover:bg-muted/40 bg-popover/20 border-t border-border transition-colors">
              {showAllResidues ? <><ChevronUp className="w-3 h-3" /> Show fewer</> : <><ChevronDown className="w-3 h-3" /> Show {interfaceResidues.length - 10} more</>}
            </button>
          )}
        </div>
      )}
    </div>
  );
};

const InteractionsPanel: React.FC<InteractionsPanelProps> = ({
  homodimer,
  heterodimers = [],
  stringInteractions = [],
  isoformId,
  geneSymbol,
  isCanonical = false,
  canonicalSequence,
  isoformSequence,
  onSelectResidues,
}) => {
    const hasComplexes = homodimer || (heterodimers && heterodimers.length > 0);
    // STRING partner table is collapsed by default so the network graph is the focus.
    const [showPartners, setShowPartners] = useState(false);

  return (
    <div className="space-y-8">
      {/* ── AFDB Complexes ───────────────────────────────────────────────── */}
      {hasComplexes ? (
        <div className="space-y-6 divide-y divide-border/50">
          {homodimer && (
            <ComplexSection
              complex={homodimer}
              isCanonical={isCanonical}
              canonicalSequence={canonicalSequence}
              isoformSequence={isoformSequence}
              onSelectResidues={onSelectResidues}
              title="AFDB Homodimer Interface"
              type="homodimer"
            />
          )}
          {heterodimers && heterodimers.map((het, idx) => (
            <ComplexSection
              key={het.id || idx}
              complex={het}
              isCanonical={isCanonical}
              canonicalSequence={canonicalSequence}
              isoformSequence={isoformSequence}
              onSelectResidues={onSelectResidues}
              title={`AFDB Heterodimer Interface (Partner: ${het.partner_uniprot_id})`}
              type="heterodimer"
            />
          ))}
        </div>
      ) : (
        <div className="bg-popover/40 border border-dashed border-border p-6 rounded-xl text-center">
          <AlertCircle className="w-8 h-8 text-slate-600 mx-auto mb-2" />
          <p className="text-sm text-muted-foreground font-medium">No complex structure available for this protein.</p>
          <p className="text-xs text-slate-600 mt-1">No AlphaFold-Multimer complex entry found in AFDB for this gene.</p>
        </div>
      )}

      {/* ── STRING Interactions ──────────────────────────────────────────────── */}
      <div className="space-y-4 pt-4 border-t border-border">
        <h4 className="text-md font-semibold text-foreground flex items-center gap-2">
          <Zap className="w-4 h-4 text-amber-400" /> Functional Network Partners (STRING)
        </h4>

        {geneSymbol && stringInteractions.length > 0 && (
          <StringNetworkGraph geneSymbol={geneSymbol} />
        )}

        {stringInteractions.length > 0 && (
          <button
            type="button"
            onClick={() => setShowPartners((v) => !v)}
            className="flex items-center gap-1 text-xs text-muted-foreground hover:text-foreground/90 transition"
          >
            {showPartners ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
            {showPartners ? "Hide" : "Show"} interaction partners ({stringInteractions.length})
          </button>
        )}

        {(showPartners || stringInteractions.length === 0) && (
        <div className="rounded-md border border-border">
          <Table>
            <TableHeader className="bg-popover/50">
              <TableRow>
                <TableHead>Partner Gene</TableHead>
                <TableHead><span className="inline-flex items-center">STRING score<InfoTip text="STRING combined score (0–1) across 7 evidence channels. Gene-level — not isoform-specific." /></span></TableHead>
                <TableHead className="text-right">Action</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {stringInteractions.length > 0 ? (
                [...stringInteractions]
                  .sort((a, b) => (b.string_score || 0) - (a.string_score || 0))
                  .map((partner, idx) => {
                    return (
                      <TableRow
                        key={idx}
                        className="border-border/50 hover:bg-muted/30"
                      >
                        <TableCell className="font-semibold text-foreground/90">
                          {partner.partner_gene_symbol}
                        </TableCell>
                        <TableCell>
                          <div className="flex items-center gap-2">
                            <div className="w-16 bg-slate-700 h-1.5 rounded-full overflow-hidden">
                              <div
                                className="h-full bg-amber-500"
                                style={{ width: `${(partner.string_score || 0) / 10}%` }}
                              />
                            </div>
                            <span className="text-xs font-mono text-foreground/80">{partner.string_score}</span>
                          </div>
                        </TableCell>
                        <TableCell className="text-right">
                          {partner.partner_gene_symbol ? (
                            <Link
                              // partner_ensp_id is not populated, so link by gene symbol via
                              // STRING's identifier query (resolves species-scoped). Include this
                              // isoform's gene too, so the network shows the interacting pair.
                              href={`https://string-db.org/cgi/network?identifiers=${encodeURIComponent(
                                [geneSymbol, partner.partner_gene_symbol].filter(Boolean).join("\r")
                              )}&species=9606`}
                              target="_blank"
                              rel="noopener noreferrer"
                              className={cn(buttonVariants({ variant: "ghost", size: "sm" }), "h-7 gap-1 text-foreground/80 hover:text-primary")}
                            >
                              <LinkIcon className="w-3 h-3" /> View
                            </Link>
                          ) : (
                            <span className="text-muted-foreground/40 text-xs">—</span>
                          )}
                        </TableCell>
                      </TableRow>
                    );
                  })
              ) : (
                <TableRow>
                  <TableCell colSpan={3} className="text-center py-8 text-muted-foreground italic">
                    No STRING interactions stored for this isoform.
                    <div className="text-[10px] mt-1 not-italic max-w-md mx-auto">
                      STRING partners are looked up via the gene&apos;s UniProt accession. Isoforms on
                      genes with no reviewed accession — or whose accession STRING has no interaction
                      data for — have no partners to show.
                    </div>
                  </TableCell>
                </TableRow>
              )}
            </TableBody>
          </Table>
        </div>
        )}
        <p className="text-[10px] text-muted-foreground/70 leading-snug">
          STRING combined scores are derived from 7 evidence channels (experiments, database,
          co-expression, text-mining, phylogenetic profile, gene fusion, neighborhood). High
          combined score doesn't guarantee physical binding — text-mining and co-expression can
          drive a 700+ score on their own.  See the Interactions &amp; Expression doc page for details.
        </p>
      </div>
    </div>
  );
};

export default InteractionsPanel;
