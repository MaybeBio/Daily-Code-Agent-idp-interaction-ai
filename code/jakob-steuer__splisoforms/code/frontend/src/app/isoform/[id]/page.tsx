"use client";

import { useEffect, useState, useMemo } from "react";
import { useParams, useRouter } from "next/navigation";
import { ArrowLeft, Check, X, Link as LinkIcon, Database, AlertCircle, Download, Box } from "lucide-react";
import { Button } from "@/components/ui/button";
import StructureViewer from "@/components/StructureViewer";
import PAEMatrix from "@/components/PAEMatrix";
import { SequenceMiniMap } from "@/components/SequenceMiniMap";
import { SequenceBar } from "@/components/SequenceBar";
import { MetadataDashboard, MetadataRecord } from "@/components/MetadataDashboard";
import { Badge } from "@/components/ui/badge";
import { Accordion, AccordionItem, AccordionTrigger, AccordionContent } from "@/components/ui/accordion";
import DomainsPanel from "@/components/analysis/DomainsPanel";
import IsoformCanonicalAlignment from "@/components/IsoformCanonicalAlignment";
import PtmAccessibilityPanel from "@/components/analysis/PtmAccessibilityPanel";
import NmdPredictionPanel from "@/components/analysis/NmdPredictionPanel";
import StructuralLandscapePanel from "@/components/analysis/StructuralLandscapePanel";
import CellTypePanel from "@/components/analysis/ExpressionInfoPanel";
import NeoantigenPanel from "@/components/analysis/NeoantigenPanel";
import InteractionsPanel from "@/components/analysis/InteractionsPanel";
import IsoformDrugPanel from "@/components/analysis/IsoformDrugPanel";
import TissueVisualizer from "@/components/TissueVisualizer";
import GeneExonArchitecture from "@/app/gene/[id]/components/GeneExonArchitecture";
import { getApiUrl } from "@/lib/api";

interface Structure {
  id: number;
  storage_uri: string;
  af3_cif_url?: string;
  af3_pae_url?: string;
  confidence_score?: number;
  avg_plddt?: number;
}

interface PPINetwork {
  id: number;
  target_uniprot_id: string;
  interaction_score: number;
  interaction_retained: boolean;
}

interface HomodimerAnalysis {
  id: number;
  isoform_id: number;
  afdb_accession?: string;
  total_interface_residues?: number;
  count_absent?: number;
  count_sasa_changed?: number;
  disruption_score?: number;
  interface_residues_json?: string;
}

interface HeterodimerAnalysis {
  id: number;
  isoform_id: number;
  partner_uniprot_id?: string;
  afdb_accession?: string;
  total_interface_residues?: number;
  count_absent?: number;
  count_sasa_changed?: number;
  disruption_score?: number;
  interface_residues_json?: string;
}


interface StringInteraction {
  id: number;
  isoform_id: number;
  partner_gene_symbol?: string;
  partner_uniprot_id?: string;
  string_score?: number;
  evidence_json?: string;
  domain_loss_possible?: boolean;
}

interface JunctionEvent {
  id: number;
  junction_type: string;
  subcategory?: string;
  exon_count_novel?: number;
  exon_count_canonical?: number;
  transcript_pos?: number;
  is_novel: boolean;
}

interface IsoformDetail {
  id: number;
  gene_id: number;
  transcript_id: string;
  transcript_name?: string | null;
  sequence?: string;
  fasta_sequence?: string;
  delta_length?: number;
  junction_coords?: string;
  structures?: Structure[];
  ppi_networks?: PPINetwork[];
  homodimer_analysis?: HomodimerAnalysis;
  heterodimer_analyses?: HeterodimerAnalysis[];
  string_interactions?: StringInteraction[];
  metadata_record?: MetadataRecord;
  gene?: {
    symbol?: string;
    uniprot_id: string;
    canonical_sequence?: string;
    canonical_length?: number;
  };
  is_canonical?: boolean;
  in_normal?: boolean | null;
  in_primary_tumor?: boolean | null;
  in_metastasis?: boolean | null;
  source?: { name?: string };
  dataset_presence?: Array<{
    source_id: number;
    source_name: string;
    in_normal: boolean;
    in_primary_tumor: boolean;
    in_metastasis: boolean;
  }>;
  pfam_domains?: any[];
  canonical_pfam_domains?: Array<{
    id: number; pfam_accession: string; pfam_name: string;
    seq_start: number; seq_end: number;
    hmm_start: number; hmm_end: number; hmm_length: number;
    domain_coverage: number; e_value: number; domain_status: string;
  }>;
  ted_domains?: any[];
  junction_events?: JunctionEvent[];
  expression_data?: Array<{
    id: number;
    isoform_id: number;
    sample: string;
    condition: string;
    fl_count: number;
    cpm?: number | null;
    isoform_frac?: number | null;
    n_cells?: number | null;
    cell_frac?: number | null;
  }>;
  structural_metrics?: {
    pae_json_path?: string;
    plddt_series?: number[];
    pae_blocks?: string;
  };
}

export default function IsoformPage() {
  const { id } = useParams();
  const router = useRouter();
  const [data, setData] = useState<IsoformDetail | null>(null);
  const [hoveredResidues, setHoveredResidues] = useState<number[]>([]);
  const [selectedResidues, setSelectedResidues] = useState<number[]>([]);
  const [calcPlddt, setCalcPlddt] = useState<number | null>(null);
  const [geneData, setGeneData] = useState<any>(null);

  useEffect(() => {
    if (data?.gene_id) {
       const apiUrl = getApiUrl();
       fetch(`${apiUrl}/api/genes/${data.gene_id}`)
         .then(res => res.json())
         .then(setGeneData)
         .catch(console.error);
    }
  }, [data?.gene_id]);

  useEffect(() => {
    if (id) {
       const apiUrl = getApiUrl();
       fetch(`${apiUrl}/api/isoform/${id}`)
         .then(res => res.json())
         .then(setData)
         .catch(console.error);
    }
  }, [id]);

  const [drugFate, setDrugFate] = useState<any[]>([]);

  // Drug-pocket fate for this isoform (panel only shown when there is data).
  useEffect(() => {
    if (!id) return;
    let alive = true;
    fetch(`${getApiUrl()}/api/isoform/${id}/drug-fate`)
      .then(r => (r.ok ? r.json() : []))
      .then(d => { if (alive) setDrugFate(Array.isArray(d) ? d : []); })
      .catch(() => {});
    return () => { alive = false; };
  }, [id]);

  const [cellTypeData, setCellTypeData] = useState<any[]>([]);

  useEffect(() => {
    if (data?.gene_id) {
      const apiUrl = getApiUrl();
      fetch(`${apiUrl}/api/gene/${data.gene_id}/cell-type`)
        .then(r => r.ok ? r.json() : [])
        .then(setCellTypeData)
        .catch(() => setCellTypeData([]));
    }
  }, [data?.gene_id]);

  const { cohortTissues, hpaTissues } = useMemo(() => {
    const cohort = new Set<string>();
    const hpa = new Set<string>();
    
    if (!data) return { cohortTissues: [], hpaTissues: [] };

    // 1. From cohort dataset presence
    if (data.dataset_presence) {
      data.dataset_presence.forEach(dp => {
        const name = dp.source_name.toLowerCase();
        if (name.includes('ccrcc') || name.includes('kidney')) cohort.add('kidney');
        if (name.includes('brca') || name.includes('bc') || name.includes('breast')) cohort.add('breast');
      });
    }

    // 2. From HPA expression data
    if (data.expression_data) {
      data.expression_data.forEach(expr => {
        if (expr.sample) hpa.add(expr.sample);
      });
    }
    
    // 3. Add all cell types from the HPA API endpoint
    cellTypeData.forEach(entry => {
      if (entry.cell_type) hpa.add(entry.cell_type);
    });

    return { cohortTissues: Array.from(cohort), hpaTissues: Array.from(hpa) };
  }, [data, cellTypeData]);

  if (!data) {
    return (
      <div className="flex items-center justify-center min-h-screen bg-background text-muted-foreground">
         <div className="flex flex-col items-center gap-4 animate-pulse">
           <Database className="w-12 h-12 text-primary" />
           <p className="text-lg text-foreground">Loading Structural Data...</p>
         </div>
      </div>
    );
  }

  const struct = data.structures?.[0];
  const cifUrl = struct?.af3_cif_url || "";
  const paeUrl = data.structural_metrics?.pae_json_path || struct?.af3_pae_url || "";
  const finalUniprotId = data.metadata_record?.uniprot_id || data.gene?.uniprot_id;
  const isPacBioIsoform = data.transcript_id?.startsWith("PB.");
  const isCanonical = data.is_canonical === true;

  const canonicalCifUrl: string | undefined = (() => {
    if (isCanonical) return undefined; 

    if (geneData?.isoforms) {
      const canIso = geneData.isoforms.find((i: any) => i.is_canonical);
      const url = canIso?.structures?.[0]?.af3_cif_url;
      if (url) return url;
    }

    const uniprotId = data.gene?.uniprot_id || data.metadata_record?.uniprot_id;
    if (uniprotId) return `/static/data/${uniprotId}.cif`;
    return undefined;
  })();

  const sequenceSimilarityToCanonical: number | null = (() => {
    const seq1 = data.gene?.canonical_sequence || "";
    const seq2 = data.sequence || "";
    if (!seq1 || !seq2 || isCanonical) return null;
    const K = 5;
    if (seq1.length < K || seq2.length < K) return 0;
    const set1 = new Set<string>();
    for (let i = 0; i <= seq1.length - K; i++) set1.add(seq1.substring(i, i + K));
    let shared = 0;
    const seen = new Set<string>();
    for (let i = 0; i <= seq2.length - K; i++) {
      const kmer = seq2.substring(i, i + K);
      if (set1.has(kmer) && !seen.has(kmer)) { shared++; seen.add(kmer); }
    }
    const total = set1.size;
    return total > 0 ? shared / total : 0;
  })();
  const canShowOverlay = sequenceSimilarityToCanonical === null || sequenceSimilarityToCanonical >= 0.05;

  const seq = data.sequence || "";

  // Why an isoform has no AF3 structure. These mirror the filters in
  // build_predict_queue.py — an isoform failing any of them is never queued, so telling the
  // user it is "queued for prediction" would be false. Of the 6,820 structureless isoforms,
  // 4,990 carry a non-standard residue, 1,168 are too short and 534 too long; only 128 are
  // genuinely eligible-but-not-yet-predicted.
  const AF3_MIN_LEN = 30;
  const AF3_MAX_LEN = 2700;
  const nonStandardAA = Array.from(new Set(seq.match(/[^ACDEFGHIKLMNPQRSTVWY]/g) ?? []));
  const structureBlockedBy: "nonstandard_aa" | "too_long" | "too_short" | null =
    nonStandardAA.length > 0 ? "nonstandard_aa"
    : seq.length > AF3_MAX_LEN ? "too_long"
    : seq.length < AF3_MIN_LEN ? "too_short"
    : null;

  const canonSeq = data.gene?.canonical_sequence;
  let dynamicJunction: [number, number] | undefined = undefined;
  
  if (seq && canonSeq && seq !== canonSeq) {
    let start = 0;
    while (start < seq.length && start < canonSeq.length && seq[start] === canonSeq[start]) {
      start++;
    }
    let endSeq = seq.length - 1;
    let endCanon = canonSeq.length - 1;
    while (endSeq >= start && endCanon >= 0 && seq[endSeq] === canonSeq[endCanon]) {
      endSeq--;
      endCanon--;
    }
    dynamicJunction = [Math.max(1, start - 1), Math.min(seq.length, endSeq + 3)];
  }

  const mid = Math.floor(seq.length / 2);
  const junctionRange: [number, number] = dynamicJunction || [Math.max(1, mid - 10), Math.min(seq.length, mid + 10)];

  const avgPlddt = calcPlddt ?? struct?.avg_plddt;
  let plddtColor = "bg-muted text-foreground/90 border-border";
  if (avgPlddt !== undefined) {
     if (avgPlddt > 90) plddtColor = "bg-[#0053D6]/20 text-[#65CBF3] border-[#0053D6]";
     else if (avgPlddt > 70) plddtColor = "bg-[#65CBF3]/20 text-[#65CBF3] border-[#65CBF3]";
     else if (avgPlddt > 50) plddtColor = "bg-[#FFDB13]/20 text-[#FFDB13] border-[#FFDB13]";
     else plddtColor = "bg-[#FF7D45]/20 text-[#FF7D45] border-[#FF7D45]";
  }

  const calculatedDelta = data.sequence && (data.metadata_record?.canonical_protein_length || data.gene?.canonical_length)
     ? data.sequence.length - (data.metadata_record?.canonical_protein_length || data.gene?.canonical_length || data.sequence.length)
     : data.delta_length;

  return (
    <div className="min-h-screen bg-background text-foreground pb-24 flex flex-col">
       
       {/* Full-width Header */}
       <header className="border-b border-border bg-card/50 backdrop-blur sticky top-0 z-50">
         <div className="w-full max-w-[1800px] mx-auto px-6 py-4 flex gap-6 items-center flex-wrap justify-between">
            <div className="flex gap-4 items-center">
               <Button variant="ghost" size="icon" onClick={() => router.back()} className="hover:bg-muted h-10 w-10 shrink-0">
                 <ArrowLeft className="w-5 h-5 text-foreground/70" />
               </Button>
                <div className="flex flex-col">
                   <div className="flex items-center gap-2">
                     {data.gene_id ? (
                       <a 
                         href={`/gene/${data.gene_id}`}
                         className="text-2xl font-bold tracking-tight bg-gradient-to-r from-[#65CBF3] to-[#0053D6] bg-clip-text text-transparent hover:opacity-80 transition-opacity"
                       >
                         {data.gene?.symbol || "Unknown Gene"}
                       </a>
                     ) : (
                       <span className="text-2xl font-bold tracking-tight bg-gradient-to-r from-[#65CBF3] to-[#0053D6] bg-clip-text text-transparent">
                         {data.gene?.symbol || "Unknown Gene"}
                       </span>
                     )}
                     <Badge variant="outline" className="text-[10px] uppercase tracking-tighter border-border text-muted-foreground py-0 h-4">
                       {data.gene_id ? "Gene" : "Gene Symbol"}
                     </Badge>
                   </div>
                  <div className="text-sm font-medium text-foreground/80 flex items-center gap-2">
                    <span className="font-mono text-xs opacity-70">{data.transcript_name || data.transcript_id}</span>
                    {data.transcript_name && data.transcript_name !== data.transcript_id && (
                      <span className="font-mono text-[10px] opacity-40">{data.transcript_id}</span>
                    )}
                    {data.metadata_record?.subcategory && (
                        <span className="capitalize opacity-60">• {data.metadata_record.subcategory.replace(/_/g, " ")}</span>
                    )}
                    {data.metadata_record?.uniprot_id && (
                        <span className="opacity-60">• UniProt: {data.metadata_record.uniprot_id}</span>
                    )}
                  </div>
                </div>
            </div>

            {/* Tissue Visualizer & Global Confidence Badge */}
            <div className="flex items-stretch gap-8 py-1">
               {/* No left border: the organ labels fan out to the left of the body figure and
                   collided with the divider. Padding is kept for spacing. */}
               <TissueVisualizer cohortTissues={cohortTissues} hpaTissues={hpaTissues} className="hidden lg:flex pl-12 pr-6 h-full items-center" />
               <div className="flex flex-col items-center border-l border-border pl-8 h-full justify-center">
                  <span className="text-xs font-semibold uppercase tracking-widest text-muted-foreground mb-1">Global pLDDT</span>
                  <div className={`w-12 h-12 flex items-center justify-center rounded-full border-2 font-bold text-lg shadow-lg ${plddtColor}`}>
                     {avgPlddt ? Math.round(avgPlddt) : 'N/A'}
                  </div>
               </div>
            </div>
         </div>
       </header>

       {/* CSS Grid Main Content */}
       <div className="flex-1 w-full p-4 lg:p-6 grid grid-cols-1 lg:grid-cols-12 gap-6 mx-auto max-w-[1800px] items-start">
          
          {/* Left Pane: 5 Columns, Sticky 3D Viewer */}
          <div className="lg:col-span-5 flex flex-col aspect-square min-h-[400px] lg:sticky lg:top-[8.5rem] z-10">
            {data.metadata_record?.folded_orf_mismatch && (
               <div className="mb-4 flex items-start gap-2 rounded-lg border border-orange-500/40 bg-orange-500/5 px-4 py-3 text-sm text-orange-400">
                 <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                 <span>
                   <span className="font-semibold">Folded ORF differs from the gene’s canonical protein.</span>{" "}
                   The sequence and AlphaFold&nbsp;3 model shown are the author-curated ORF folded for this transcript. Its protein
                   {data.metadata_record.folded_best_ref
                     ? <> best-matches UniProt <span className="font-mono text-foreground/90">{data.metadata_record.folded_best_ref.split("|").pop()}</span>
                        {data.metadata_record.folded_best_ref_identity != null && ` (${data.metadata_record.folded_best_ref_identity.toFixed(0)}%)`}</>
                     : <> does not match a reviewed protein</>}
                   , not canonical {data.gene?.symbol}. The transcript maps to the {data.gene?.symbol} locus, but the model depicts that ORF’s product — kept for faithfulness to the source data and flagged here for transparency.
                 </span>
               </div>
            )}
            {!isCanonical && data.gene?.canonical_sequence && data.sequence && sequenceSimilarityToCanonical !== null && sequenceSimilarityToCanonical < 0.05 && (
               <div className="mb-4 flex items-start gap-2 rounded-lg border border-amber-500/30 bg-amber-500/5 px-4 py-3 text-sm text-amber-400">
                 <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                 <span>
                   <span className="font-semibold">Very low sequence similarity with canonical</span>
                   {" "}({`${(sequenceSimilarityToCanonical * 100).toFixed(0)}%`} shared 5-mers).
                   This isoform likely represents a distinct protein product of {data.gene?.symbol} with a different reading frame or exon composition.
                   3D structure superposition is suppressed.
                 </span>
               </div>
            )}
            
            {cifUrl ? (
              <StructureViewer
                 cifUrl={cifUrl}
                 canonicalCifUrl={!isCanonical && canShowOverlay ? canonicalCifUrl : undefined}
                 hoveredResidues={hoveredResidues}
                 selectedResidues={selectedResidues}
                 junctionRange={junctionRange}
                 sequence={data.sequence}
                 canonicalSequence={(data.gene as any)?.canonical_sequence || undefined}
                 onPlddtCalculated={setCalcPlddt}
                 primaryLabel={isCanonical ? "Canonical Isoform" : "Novel Isoform"}
              />
            ) : structureBlockedBy === "nonstandard_aa" ? (
              <div className="flex-1 flex flex-col items-center justify-center text-center p-12 bg-background/40 backdrop-blur-md rounded-xl border border-border">
                <AlertCircle className="w-16 h-16 text-amber-500 mb-6" />
                <h4 className="text-xl font-bold text-foreground/90 mb-2">Cannot Be Predicted — Non-Standard Residues</h4>
                <p className="text-muted-foreground max-w-sm">
                  This sequence contains{" "}
                  <span className="font-mono text-foreground/80">{nonStandardAA.join(", ")}</span>
                  {nonStandardAA.includes("X") ? " (an undetermined residue)" : ""}, which AlphaFold 3
                  cannot model. It is <span className="text-foreground/80">not queued</span> for
                  prediction, and no structure will become available unless the underlying sequence
                  is resolved.
                </p>
              </div>
            ) : structureBlockedBy === "too_long" ? (
              <div className="flex-1 flex flex-col items-center justify-center text-center p-12 bg-background/40 backdrop-blur-md rounded-xl border border-border">
                <AlertCircle className="w-16 h-16 text-amber-500 mb-6" />
                <h4 className="text-xl font-bold text-foreground/90 mb-2">Sequence Too Long for AF3</h4>
                <p className="text-muted-foreground max-w-sm">
                  At {seq.length.toLocaleString()} amino acids this isoform exceeds the {AF3_MAX_LEN.toLocaleString()} AA
                  limit of our AlphaFold 3 pipeline, so it is{" "}
                  <span className="text-foreground/80">not queued</span> for prediction.
                </p>
              </div>
            ) : structureBlockedBy === "too_short" ? (
              <div className="flex-1 flex flex-col items-center justify-center text-center p-12 bg-background/40 backdrop-blur-md rounded-xl border border-border">
                <AlertCircle className="w-16 h-16 text-amber-500 mb-6" />
                <h4 className="text-xl font-bold text-foreground/90 mb-2">Sequence Too Short for AF3</h4>
                <p className="text-muted-foreground max-w-sm">
                  At {seq.length} amino acids this isoform falls below the {AF3_MIN_LEN} AA minimum of our
                  AlphaFold 3 pipeline, so it is <span className="text-foreground/80">not queued</span> for
                  prediction.
                </p>
              </div>
            ) : (
              <div className="flex-1 flex flex-col items-center justify-center text-center p-12 bg-background/40 backdrop-blur-md rounded-xl border border-border">
                <Box className="w-16 h-16 text-slate-700 mb-6 animate-pulse" />
                <h4 className="text-xl font-bold text-foreground/90 mb-2">Structural Model Pending</h4>
                <p className="text-muted-foreground max-w-sm">
                  This transcript is eligible for AlphaFold 3 prediction and is queued for an upcoming
                  batch. The 3D model will appear here once it has been computed.
                </p>
              </div>
            )}
          </div>

          {/* Right Pane: 7 Columns, Vertical Dashboard */}
          <div className="lg:col-span-7 flex flex-col gap-4">
             <Accordion 
                multiple 
                defaultValue={["details", "sequence", "domains", "dssp", "ptm", "neoantigen"]} 
                className="w-full space-y-4"
             >
                {/* 1. Details & Provenance */}
                <AccordionItem value="details" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Overview</AccordionTrigger>
                   <AccordionContent className="pb-4">
                     <MetadataDashboard metadata={data.metadata_record} fallbackId={finalUniprotId} isoform={data} />
                   </AccordionContent>
                </AccordionItem>

                {/* 2. Sequence & Topology */}
                <AccordionItem value="sequence" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">
                     Sequence & Domain Topology
                   </AccordionTrigger>
                   <AccordionContent className="pb-4">
                     <div className="flex items-center justify-start w-full mb-2">
                       <Badge className="bg-primary/20 text-primary hover:bg-primary/30 border-none">Hover/Click sequence to Focus 3D</Badge>
                     </div>
                     <SequenceBar 
                       sequence={seq} 
                       junctionEvents={data.junction_events}
                       hoveredResidues={hoveredResidues}
                       selectedResidues={selectedResidues}
                       onHover={setHoveredResidues}
                       onSelectionChange={setSelectedResidues}
                     />
                     <div className="mt-4">
                        <SequenceMiniMap sequence={seq} junctionEvents={data.junction_events} geneData={geneData} currentIsoformId={data.id} />
                     </div>
                     {!isCanonical && data.gene?.canonical_sequence && data.sequence && (
                       <>
                         <div className="mt-4">
                           <IsoformCanonicalAlignment
                             canonicalSequence={data.gene.canonical_sequence}
                             isoformSequence={data.sequence}
                             canonicalLabel={data.gene?.symbol || "Canonical"}
                             isoformLabel={data.transcript_name || data.transcript_id}
                             hoveredResidues={hoveredResidues}
                             selectedResidues={selectedResidues}
                             onHover={setHoveredResidues}
                             onSelectionChange={setSelectedResidues}
                           />
                         </div>
                       </>
                     )}

                     {/* Moved Exon Architecture to the bottom of the Sequence & Topology panel */}
                     {geneData && data.id && (
                       <div className="mt-8 border-t border-border/50 pt-4 -mb-2">
                         <GeneExonArchitecture 
                             gene={geneData} 
                             selectedIsoformIds={[geneData.isoforms?.find((i: any) => i.is_canonical)?.id, data.id].filter(Boolean)} 
                             hideHeader={true} 
                             hideCardBackground={true} 
                         />
                       </div>
                     )}
                   </AccordionContent>
                </AccordionItem>

                {/* 3. Functional Domains */}
                <AccordionItem value="domains" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Functional Domains</AccordionTrigger>
                   <AccordionContent className="pb-4">
                      <DomainsPanel 
                        pfamDomains={data.pfam_domains} 
                        tedDomains={data.ted_domains} 
                        seqLength={data.sequence?.length || 0}
                        canonicalPfamDomains={data.canonical_pfam_domains}
                        canonicalSeqLength={data.gene?.canonical_length}
                        isoformSeq={data.sequence}
                        canonicalSeq={data.gene?.canonical_sequence}
                        hoveredResidues={hoveredResidues}
                        selectedResidues={selectedResidues}
                        onHover={setHoveredResidues}
                        onSelect={setSelectedResidues}
                      />
                   </AccordionContent>
                </AccordionItem>

                {/* 4. Structural Landscape */}
                <AccordionItem value="dssp" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Structural Landscape</AccordionTrigger>
                   <AccordionContent className="pb-4">
                      <StructuralLandscapePanel 
                        isoformId={data.id} 
                        hoveredResidues={hoveredResidues}
                        selectedResidues={selectedResidues}
                        onHover={setHoveredResidues}
                        onSelect={setSelectedResidues}
                      />
                   </AccordionContent>
                </AccordionItem>

                {/* 5. PTM Accessibility */}
                <AccordionItem value="ptm" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">PTM Accessibility</AccordionTrigger>
                   <AccordionContent className="pb-4">
                     <PtmAccessibilityPanel
                       isoformId={data.id}
                       isCanonical={isCanonical}
                       canonicalSeq={!isCanonical ? (data.gene?.canonical_sequence ?? undefined) : undefined}
                       isoformSeq={!isCanonical ? (data.sequence ?? undefined) : undefined}
                       hoveredResidues={hoveredResidues}
                       selectedResidues={selectedResidues}
                       onHover={setHoveredResidues}
                       onSelect={setSelectedResidues}
                     />
                   </AccordionContent>
                </AccordionItem>

                {/* 6. Protein Interactions */}
                <AccordionItem value="interactions" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Protein Interactions</AccordionTrigger>
                   <AccordionContent className="pb-4">
                     <InteractionsPanel
                       homodimer={data.homodimer_analysis}
                       heterodimers={data.heterodimer_analyses}
                       stringInteractions={data.string_interactions}
                       isoformId={data.id}
                       geneSymbol={(data.gene as any)?.symbol ?? undefined}
                       isCanonical={isCanonical}
                       canonicalSequence={(data.gene as any)?.canonical_sequence ?? undefined}
                       isoformSequence={data.sequence ?? undefined}
                       onSelectResidues={setSelectedResidues}
                     />
                   </AccordionContent>
                </AccordionItem>

                {/* 7. Drug Pockets — only when this isoform has drug-fate data */}
                {drugFate.length > 0 && (
                  <AccordionItem value="drugs" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                     <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Drug Pockets</AccordionTrigger>
                     <AccordionContent className="pb-4">
                       <IsoformDrugPanel
                         fates={drugFate}
                         geneSymbol={geneData?.symbol ?? data.gene?.symbol}
                         onSelect={setSelectedResidues}
                       />
                     </AccordionContent>
                  </AccordionItem>
                )}

                {/* 8. Neoantigens */}
                <AccordionItem value="neoantigen" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Neoantigens</AccordionTrigger>
                   <AccordionContent className="pb-4">
                      <NeoantigenPanel
                        isoformId={data.id}
                        isCanonical={isCanonical}
                        hoveredResidues={hoveredResidues}
                        selectedResidues={selectedResidues}
                        onHover={setHoveredResidues}
                        onSelect={setSelectedResidues}
                      />
                   </AccordionContent>
                </AccordionItem>

                {/* 8. NMD Prediction */}
                <AccordionItem value="nmd" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">NMD Prediction</AccordionTrigger>
                   <AccordionContent className="pb-4">
                     <NmdPredictionPanel isoformId={data.id} />
                   </AccordionContent>
                </AccordionItem>

                {/* 9. Expression & Tissues */}
                <AccordionItem value="expression" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Expression & Tissues</AccordionTrigger>
                   <AccordionContent className="pb-4">
                     <CellTypePanel
                       geneId={data.gene_id}
                       expressionData={data.expression_data ?? []}
                       tissueFlags={{
                         in_normal: data.in_normal,
                         in_primary_tumor: data.in_primary_tumor,
                         in_metastasis: data.in_metastasis,
                       }}
                       sourceName={data.source?.name}
                       sourceMemberships={data.dataset_presence}
                     />
                   </AccordionContent>
                </AccordionItem>

                {/* 10. Predicted Aligned Error (PAE) */}
                <AccordionItem value="pae" className="bg-card/50 border border-border rounded-xl px-6 border-b-0 shadow-sm overflow-hidden">
                   <AccordionTrigger className="text-lg font-semibold hover:no-underline py-4">Predicted Aligned Error (PAE)</AccordionTrigger>
                   <AccordionContent className="pb-4">
                     {paeUrl ? (
                        <div className="w-full aspect-square max-w-[640px] mx-auto">
                          <PAEMatrix
                             paeUrl={paeUrl} 
                             sequenceLength={seq.length}
                             onHover={(res) => {
                                setHoveredResidues(res);
                                setSelectedResidues(res);
                             }} 
                          />
                        </div>
                     ) : (
                        <div className="w-full h-full flex flex-col items-center justify-center bg-card/30 rounded-lg border border-border text-muted-foreground p-8 text-center min-h-[300px]">
                          <AlertCircle className="w-8 h-8 mb-2 opacity-50" />
                          <p>PAE Matrix Data Unavailable for this Structure</p>
                        </div>
                     )}
                   </AccordionContent>
                </AccordionItem>

             </Accordion>
          </div>
       </div>
    </div>
  );
}
