"use client";

import { useState, useEffect, useMemo } from "react";
import dynamic from "next/dynamic";
import IsoformTable from "./IsoformTable";
import GeneSequenceAlignment from "./GeneSequenceAlignment";
import GeneAnalysisSummary from "./GeneAnalysisSummary";
import GeneExpressionHeatmap from "./GeneExpressionHeatmap";
import GeneExonArchitecture from "./GeneExonArchitecture";
import GeneDrugPanel from "./GeneDrugPanel";
import GeneStringInteractions from "./GeneStringInteractions";
import { getApiUrl } from "@/lib/api";
import {
  Accordion,
  AccordionContent,
  AccordionItem,
  AccordionTrigger,
} from "@/components/ui/accordion";
import { AlertCircle } from "lucide-react";

// Mol* (molstar) loads browser-only React contexts at module-eval time.
// Importing GeneStructureViewerClient directly causes SSR to crash with
//   TypeError: Cannot read properties of undefined (reading 'Provider')
// at module evaluation, because Mol*'s plugin-UI context module assumes a
// browser environment. Wrap in next/dynamic with ssr:false so the entire
// transitive Mol* graph is excluded from the server bundle — the same
// pattern StructureViewer.tsx uses for the isoform page.
const GeneStructureViewerClient = dynamic(
  () => import("./GeneStructureViewerClient"),
  {
    ssr: false,
    loading: () => (
      <div className="relative w-full h-full min-h-[400px] bg-slate-900/50 rounded-lg overflow-hidden border border-slate-800 flex items-center justify-center">
        <span className="text-muted-foreground text-sm">Loading 3D viewer…</span>
      </div>
    ),
  }
);

// "id" sorts the raw source identifier (ENST…/PB.…); "name" sorts the readable
// transcript_name (MET-201). They are separate columns and sort differently.
export type SortKey = "id" | "name" | "length" | "source";
export type SortOrder = "asc" | "desc";

export default function GeneDashboardClient({ gene }: { gene: any }) {
  const [selectedIsoformIds, setSelectedIsoformIds] = useState<number[]>([]);

  // Sort state lifted here so all panels share the same order as the table
  const [sortKey, setSortKey] = useState<SortKey>("name");
  const [sortOrder, setSortOrder] = useState<"asc" | "desc">("asc");
  const [hoveredResidues, setHoveredResidues] = useState<number[]>([]);
  const [selectedResidues, setSelectedResidues] = useState<number[]>([]);
  const [drugData, setDrugData] = useState<any | null>(null);

  // Drug complexes for this gene (panel only shown when the gene is drugged).
  useEffect(() => {
    if (!gene?.id) return;
    let alive = true;
    fetch(`${getApiUrl()}/api/gene/${gene.id}/drugs`)
      .then(r => (r.ok ? r.json() : null))
      .then(d => { if (alive) setDrugData(d); })
      .catch(() => {});
    return () => { alive = false; };
  }, [gene?.id]);

  // Restore selection from URL hash on mount
  useEffect(() => {
    if (typeof window === "undefined") return;
    const hash = window.location.hash;
    const match = hash.match(/isoforms=([0-9,]+)/);
    if (match) {
      const ids = match[1].split(",").map(Number).filter(Boolean);
      setSelectedIsoformIds(ids);
    }
  }, []);

  // Persist selection in URL hash whenever it changes
  useEffect(() => {
    if (typeof window === "undefined") return;
    window.history.replaceState(null, "",
      selectedIsoformIds.length > 0
        ? `${window.location.pathname}#isoforms=${selectedIsoformIds.join(",")}`
        : window.location.pathname
    );
  }, [selectedIsoformIds]);

  const toggleIsoformSelection = (id: number) => {
    setSelectedIsoformIds(prev =>
      prev.includes(id) ? prev.filter(i => i !== id) : [...prev, id]
    );
  };

  // Canonical isoform reference
  const canonicalIso = gene.isoforms?.find((iso: any) => iso.id === gene.canonical_transcript_id)
                       || gene.isoforms?.find((iso: any) => iso.is_canonical);

  // Authoritative sort order — mirrors IsoformTable's sort logic (no text/source filter).
  // Every panel that shows a list of isoforms consumes this instead of deriving its own order.
  const sortedIsoformIds: number[] = useMemo(() => {
    const all: any[] = gene.isoforms ?? [];
    return [...all]
      .sort((a, b) => {
        // 1. Canonical always first
        if (a.is_canonical !== b.is_canonical) return a.is_canonical ? -1 : 1;
        // 2. User sort key
        let cmp = 0;
        if (sortKey === "id") {
          cmp = a.transcript_id.localeCompare(b.transcript_id);
        } else if (sortKey === "length") {
          cmp = (a.sequence?.length || 0) - (b.sequence?.length || 0);
        } else if (sortKey === "source") {
          cmp = (a.source?.name || "Unknown").localeCompare(b.source?.name || "Unknown");
        }
        if (cmp !== 0) return sortOrder === "asc" ? cmp : -cmp;
        // 3. Stable fallback: flagged first, then transcript_id
        if (a.is_flagged !== b.is_flagged) return a.is_flagged ? -1 : 1;
        return a.transcript_id.localeCompare(b.transcript_id);
      })
      .map((iso: any) => iso.id);
  }, [gene.isoforms, sortKey, sortOrder]);

  // Selected isoforms in table order (used by panels that need objects, not just IDs)
  const selectedIsoforms = useMemo(() => {
    const selectedSet = new Set(selectedIsoformIds);
    const isoById = new Map((gene.isoforms ?? []).map((i: any) => [i.id, i]));
    return sortedIsoformIds
      .filter(id => selectedSet.has(id))
      .map(id => isoById.get(id))
      .filter(Boolean);
  }, [sortedIsoformIds, selectedIsoformIds, gene.isoforms]);

  const mappedSelections = useMemo(() => {
    return selectedIsoforms.map((iso: any) => ({
      id: iso.id,
      transcript_id: iso.transcript_id,
      transcript_name: iso.transcript_name,   // readable label for viewer + alignment
      sequence: iso.sequence,
      cifUrl: iso.structures?.[0]?.af3_cif_url,
    }));
  }, [selectedIsoforms]);

  // Structural reference for the 3D viewer: the canonical when it has a structure,
  // otherwise fall back to the first isoform (table order) that does — so isoforms
  // still render when the canonical is too long for AF3 (e.g. >2,700 aa) and has no
  // model. The reference is excluded from the overlay list to avoid double-loading.
  const viewerRef = useMemo(() => {
    if (canonicalIso?.structures?.[0]?.af3_cif_url) return canonicalIso;
    const byId = new Map((gene.isoforms ?? []).map((i: any) => [i.id, i]));
    for (const id of sortedIsoformIds) {
      const i: any = byId.get(id);
      if (i?.structures?.[0]?.af3_cif_url) return i;
    }
    return undefined;
  }, [gene.isoforms, canonicalIso, sortedIsoformIds]);

  const refCifUrl: string | undefined = viewerRef?.structures?.[0]?.af3_cif_url;
  const refSequence: string | undefined =
    viewerRef?.sequence || (viewerRef?.id === canonicalIso?.id ? gene.canonical_sequence : undefined);
  const usingFallbackRef = !!viewerRef && viewerRef.id !== canonicalIso?.id;
  const viewerSelections = useMemo(
    () => mappedSelections.filter((s: any) => s.id !== viewerRef?.id),
    [mappedSelections, viewerRef]);

  return (
    <div className="flex flex-col gap-8 w-full mt-8 animate-in slide-in-from-bottom-4 duration-700">
       <div className="grid grid-cols-1 lg:grid-cols-12 gap-8 items-start">
         
         {/* Left Pane: 5 Columns, Sticky 3D Viewer */}
         <div className="lg:col-span-5 flex flex-col aspect-square min-h-[400px] lg:sticky lg:top-[6.5rem] z-10">
            {!canonicalIso?.structures?.[0]?.af3_cif_url && (
              <div className="mb-4 flex items-start gap-3 rounded-lg border border-amber-500/30 bg-amber-500/5 px-4 py-3 text-sm text-amber-400">
                 <AlertCircle className="w-5 h-5 mt-0.5 shrink-0" />
                 <div>
                   {(canonicalIso?.sequence || gene.canonical_sequence)?.length > 2700 ? (
                     <>
                       <div className="font-semibold mb-1">Canonical Sequence Too Long for AF3</div>
                       <p className="text-amber-500/90 leading-relaxed text-xs">
                         The canonical sequence is {(canonicalIso?.sequence || gene.canonical_sequence)?.length} amino acids long, exceeding the maximum length limit (2,700 AAs) for our AlphaFold 3 prediction pipeline.
                       </p>
                     </>
                   ) : (
                     <div className="font-semibold mb-1">Canonical structure unavailable</div>
                   )}
                   {usingFallbackRef && (
                     <p className="text-amber-500/90 leading-relaxed text-xs mt-1">
                       Showing <span className="font-mono">{viewerRef?.transcript_name || viewerRef?.transcript_id}</span> as the structural reference instead.
                     </p>
                   )}
                 </div>
              </div>
            )}
            <GeneStructureViewerClient
               canonicalCifUrl={refCifUrl}
               canonicalSequence={refSequence}
               selectedIsoforms={viewerSelections}
               hoveredResidues={hoveredResidues}
               selectedResidues={selectedResidues}
            />
         </div>

         {/* Right Pane: 7 Columns, Vertical Dashboard */}
         <div className="lg:col-span-7 flex flex-col gap-8">
            {/* Table stays outside Accordion to act as a sticky control center */}
            <IsoformTable
               isoforms={gene.isoforms}
               selectedIds={selectedIsoformIds}
               onToggleSelect={toggleIsoformSelection}
               onSetSelection={setSelectedIsoformIds}
               sortKey={sortKey}
               sortOrder={sortOrder}
               onSortChange={(key, order) => { setSortKey(key); setSortOrder(order); }}
            />

            <Accordion 
                multiple 
                defaultValue={["exon", "msa"]} 
                className="w-full space-y-4"
            >
              {/* 1. Exon Architecture */}
              <AccordionItem value="exon" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                 <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                   Exon Architecture
                 </AccordionTrigger>
                 <AccordionContent className="pb-4">
                    <GeneExonArchitecture
                       gene={gene}
                       selectedIsoformIds={selectedIsoformIds}
                       sortedIsoformIds={sortedIsoformIds}
                    />
                 </AccordionContent>
              </AccordionItem>

              {/* 2. Sequence Alignment (MSA) */}
              <AccordionItem value="msa" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                 <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                   Sequence Alignment (MSA)
                 </AccordionTrigger>
                 <AccordionContent className="pb-4">
                    <GeneSequenceAlignment
                       canonicalId={canonicalIso?.transcript_id}
                       canonicalSequence={canonicalIso?.sequence || gene.canonical_sequence}
                       selectedIsoforms={mappedSelections}
                       hoveredResidues={hoveredResidues}
                       selectedResidues={selectedResidues}
                       onHover={setHoveredResidues}
                       onSelectionChange={setSelectedResidues}
                    />
                 </AccordionContent>
              </AccordionItem>

              {/* 3. Functional Impact */}
              <AccordionItem value="functional" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                 <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                   Functional Impact
                 </AccordionTrigger>
                 <AccordionContent className="pb-4">
                    <GeneAnalysisSummary
                       gene={gene}
                       selectedIsoformIds={selectedIsoformIds}
                       sortedIsoformIds={sortedIsoformIds}
                    />
                 </AccordionContent>
              </AccordionItem>

              {/* 4. Drug Interactions — only for drugged genes */}
              {drugData?.complexes?.length > 0 && (
                <AccordionItem value="drugs" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                     Drug Interactions
                   </AccordionTrigger>
                   <AccordionContent className="pb-4">
                      <GeneDrugPanel data={drugData} gene={gene} />
                   </AccordionContent>
                </AccordionItem>
              )}

              {/* 5. Protein Interactions — gene-level STRING network */}
              {gene.symbol && (
                <AccordionItem value="interactions" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                     Protein Interactions
                   </AccordionTrigger>
                   <AccordionContent className="pb-4">
                      <GeneStringInteractions geneSymbol={gene.symbol} />
                   </AccordionContent>
                </AccordionItem>
              )}

              {/* 6. Tissue Detection */}
              <AccordionItem value="tissue" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                 <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                   Tissue Detection
                 </AccordionTrigger>
                 <AccordionContent className="pb-4">
                    <GeneExpressionHeatmap
                       gene={gene}
                       selectedIsoformIds={selectedIsoformIds}
                       sortedIsoformIds={sortedIsoformIds}
                    />
                 </AccordionContent>
              </AccordionItem>
            </Accordion>
         </div>
       </div>
    </div>
  );
}
