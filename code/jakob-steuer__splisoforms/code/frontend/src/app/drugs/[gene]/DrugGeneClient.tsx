"use client";

import { useMemo, useState } from "react";
import Link from "next/link";
import dynamic from "next/dynamic";
import { ArrowLeft, Pill, FlaskConical, Beaker, ImageOff, Network, Info, ExternalLink, HelpCircle } from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { cn } from "@/lib/utils";

function InfoTip({ text }: { text: string }) {
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger>
          <HelpCircle className="w-3 h-3 text-muted-foreground/50 hover:text-muted-foreground transition-colors cursor-help" />
        </TooltipTrigger>
        <TooltipContent className="bg-popover border-border text-foreground text-xs max-w-xs shadow-xl font-normal normal-case tracking-normal">
          {text}
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}
import DrugSummaryPanel, { type DrugSummary } from "@/components/DrugSummaryPanel";
import DrugGeneOverviewPanel from "@/components/DrugGeneOverviewPanel";

const StructureViewer = dynamic(() => import("@/components/StructureViewer"), {
  ssr: false,
  loading: () => (
    <div className="h-full w-full min-h-[320px] bg-muted/40 animate-pulse rounded-xl border border-border/50 flex items-center justify-center text-xs text-muted-foreground">
      Loading 3D viewer…
    </div>
  ),
});

// ─── types (loose — mirror API JSON) ─────────────────────────────────────────
interface Complex {
  id: number; job_name: string; drug: DrugSummary;
  iptm?: number; holo_grade?: string; cif_url?: string;
  contact_residues_json?: string; n_contact_residues?: number;
  target_isoform_id?: number; pocket_plddt?: number;
}
interface Fate {
  drug_complex_id: number; isoform_id: number; drug_name?: string;
  n_pocket_residues?: number; n_retained?: number; frac_retained?: number;
  fate?: string; mechanism_class?: string; decidable?: boolean;
  pocket_plddt?: number; retained_residues_json?: string;
}
interface Isoform {
  id: number; transcript_id: string; transcript_name?: string | null; sequence?: string;
  is_canonical?: boolean; is_mane_select?: boolean;
  structures?: { af3_cif_url?: string }[];
}

/** Readable handle (MET-206-NOV-a6d00a), falling back to the raw source id. */
const isoLabel = (i?: Isoform | null) => i?.transcript_name || i?.transcript_id || "isoform";

const FATE_STYLE: Record<string, { label: string; cls: string; dot: string }> = {
  lost:       { label: "Lost",       cls: "bg-red-500/15 text-red-400 border-red-500/30",       dot: "bg-red-500" },
  disrupted:  { label: "Disrupted",  cls: "bg-amber-500/15 text-amber-400 border-amber-500/30",  dot: "bg-amber-500" },
  preserved:  { label: "Preserved",  cls: "bg-emerald-500/15 text-emerald-400 border-emerald-500/30", dot: "bg-emerald-500" },
  no_call:    { label: "No call",    cls: "bg-slate-500/15 text-slate-400 border-slate-500/30",  dot: "bg-slate-500" },
};
const gradeCls = (g?: string) =>
  g === "A" ? "bg-emerald-500/15 text-emerald-400 border-emerald-500/30"
  : g === "B" ? "bg-amber-500/15 text-amber-400 border-amber-500/30"
  : "bg-slate-500/15 text-slate-400 border-slate-500/30";

// pLDDT confidence circle — mirrors the isoform page's Global pLDDT badge
// (label above, coloured circle with the rounded value).
function PlddtCircle({ label, value }: { label: string; value?: number | null }) {
  const v = value == null ? null : Math.round(value);
  const color = v == null ? "bg-muted text-foreground/70 border-border"
    : v > 90 ? "bg-[#0053D6]/20 text-[#65CBF3] border-[#0053D6]"
    : v > 70 ? "bg-[#65CBF3]/20 text-[#65CBF3] border-[#65CBF3]"
    : v > 50 ? "bg-[#FFDB13]/20 text-[#FFDB13] border-[#FFDB13]"
    : "bg-[#FF7D45]/20 text-[#FF7D45] border-[#FF7D45]";
  return (
    <div className="flex flex-col items-center">
      <span className="text-[10px] font-semibold uppercase tracking-widest text-muted-foreground mb-1 whitespace-nowrap">{label}</span>
      <div className={cn("w-11 h-11 flex items-center justify-center rounded-full border-2 font-bold text-base shadow-md", color)}>
        {v ?? "N/A"}
      </div>
    </div>
  );
}

export default function DrugGeneClient({ gene, drugs }: { gene: any; drugs: any }) {
  const complexes: Complex[] = drugs.complexes ?? [];
  const fates: Fate[] = drugs.fates ?? [];

  // complexes sorted by grade (A>B>C) then drug name
  const sortedComplexes = useMemo(
    () => [...complexes].sort((a, b) =>
      (a.holo_grade ?? "Z").localeCompare(b.holo_grade ?? "Z") ||
      a.drug.name.localeCompare(b.drug.name)),
    [complexes]);

  // the drugged (canonical) transcript(s); siblings = everything else
  const targetIds = useMemo(() => new Set(complexes.map(c => c.target_isoform_id)), [complexes]);
  const siblings: Isoform[] = useMemo(
    () => (gene.isoforms ?? [])
      .filter((i: Isoform) => !targetIds.has(i.id))
      .sort((a: Isoform, b: Isoform) => (b.sequence?.length ?? 0) - (a.sequence?.length ?? 0)),
    [gene.isoforms, targetIds]);

  const [complexId, setComplexId] = useState<number>(sortedComplexes[0]?.id);
  const selectedComplex = sortedComplexes.find(c => c.id === complexId) ?? sortedComplexes[0];

  // fates for the selected drug, keyed by isoform
  const fateByIso = useMemo(() => {
    const m = new Map<number, Fate>();
    for (const f of fates) if (f.drug_complex_id === selectedComplex?.id) m.set(f.isoform_id, f);
    return m;
  }, [fates, selectedComplex]);

  // default isoform: prefer a subtractive one so the view shows an effect
  const defaultIso = useMemo(() => {
    const sub = siblings.find(i => ["lost", "disrupted"].includes(fateByIso.get(i.id)?.fate ?? ""));
    return sub ?? siblings[0];
  }, [siblings, fateByIso]);
  const [isoId, setIsoId] = useState<number | undefined>(defaultIso?.id);
  const selectedIso = siblings.find(i => i.id === isoId) ?? defaultIso;
  const fate = selectedIso ? fateByIso.get(selectedIso.id) : undefined;

  // global pLDDT (from the viewers) + the drugged/canonical sequence for overlay alignment
  const [canonPlddt, setCanonPlddt] = useState<number | null>(null);
  const [isoPlddt, setIsoPlddt] = useState<number | null>(null);
  const canonicalSeq = useMemo(
    () => (gene.isoforms ?? []).find((i: Isoform) => i.id === selectedComplex?.target_isoform_id)?.sequence,
    [gene.isoforms, selectedComplex]);

  const pocketResidues: number[] = useMemo(() => {
    try { return JSON.parse(selectedComplex?.contact_residues_json ?? "[]"); } catch { return []; }
  }, [selectedComplex]);
  const retainedResidues: number[] = useMemo(() => {
    try { return JSON.parse(fate?.retained_residues_json ?? "[]"); } catch { return []; }
  }, [fate]);

  const isoCif = selectedIso?.structures?.[0]?.af3_cif_url;
  const fs = FATE_STYLE[fate?.fate ?? "no_call"] ?? FATE_STYLE.no_call;

  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-[1700px] mx-auto p-5 md:p-8">
        {/* Header */}
        <header className="flex items-start gap-3 mb-6">
          <Link href="/explore/drugs">
            <Button variant="ghost" size="icon" className="hover:bg-emerald-500/10 shrink-0">
              <ArrowLeft className="w-5 h-5" />
            </Button>
          </Link>
          <Pill className="w-6 h-6 text-emerald-500 shrink-0 mt-1.5" />
          <div className="min-w-0">
            <h1 className="text-3xl font-bold tracking-tight">
              {gene.symbol} <span className="text-muted-foreground font-normal text-xl">· drug pockets</span>
            </h1>

            {gene.protein_name && (
              <p className="text-sm text-foreground/70 leading-snug mt-0.5">{gene.protein_name}</p>
            )}

            {/* Counts live in the Overview panel now — the header keeps only identity. */}
            <p className="text-sm text-muted-foreground mt-1">
              <Link href={`/gene/${gene.id}`} className="text-primary hover:underline">gene overview →</Link>
            </p>

            {/* External refs, mirroring the gene page header */}
            <div className="flex items-center gap-3 flex-wrap mt-2 text-xs">
              {gene.uniprot_id && (
                <a href={`https://www.uniprot.org/uniprotkb/${gene.uniprot_id}/entry`}
                   target="_blank" rel="noopener noreferrer"
                   className="inline-flex items-center gap-1 text-muted-foreground hover:text-primary transition-colors">
                  UniProt <span className="font-mono">{gene.uniprot_id}</span>
                  <ExternalLink className="w-2.5 h-2.5 opacity-60" />
                </a>
              )}
              {gene.ensembl_gene_id && (
                <a href={`https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${gene.ensembl_gene_id}`}
                   target="_blank" rel="noopener noreferrer"
                   className="inline-flex items-center gap-1 text-muted-foreground hover:text-primary transition-colors">
                  Ensembl <span className="font-mono">{gene.ensembl_gene_id}</span>
                  <ExternalLink className="w-2.5 h-2.5 opacity-60" />
                </a>
              )}
            </div>

          </div>
        </header>

        {/* Two aligned ROWS rather than two long columns.
            Row 1 pairs the drugged-complex viewer with the panels that describe the DRUG
            (Overview / selector / chemistry); row 2 pairs the isoform viewer with the panels
            that describe THAT ISOFORM's fate. Because each row is its own grid, "Isoform
            fate" is guaranteed to start level with the isoform viewer.
            The viewer cards stretch to their row's height (flex-1 content instead of a fixed
            360px): the right-hand stack is the taller side, so a fixed height left ~870px of
            dead space under the left column. Filling it makes the structures bigger, which is
            the point of the page. */}
        <div className="grid grid-cols-1 lg:grid-cols-[minmax(0,1.4fr)_minmax(0,1fr)] gap-5 items-stretch mb-5">
          {/* ── ROW 1 LEFT: the drugged complex ───────────────────────────── */}
          <Card className="overflow-hidden flex flex-col">
              <CardHeader className="py-4">
                <div className="flex items-start justify-between gap-3">
                  <CardTitle className="text-lg font-semibold flex items-center gap-2 flex-wrap">
                    <Beaker className="w-4 h-4 text-emerald-500" />
                    Canonical + {selectedComplex?.drug.name}
                    <Badge className={cn("ml-1 border", gradeCls(selectedComplex?.holo_grade))}>
                      grade {selectedComplex?.holo_grade ?? "–"} · ipTM {selectedComplex?.iptm ?? "–"}
                    </Badge>
                    <InfoTip text="Co-fold confidence from the interface ipTM: A ≥ 0.75, B 0.50–0.75, C < 0.50." />
                  </CardTitle>
                  <div className="flex items-start gap-4 shrink-0">
                    <PlddtCircle label="Global pLDDT" value={canonPlddt} />
                    <PlddtCircle label="Pocket pLDDT" value={selectedComplex?.pocket_plddt} />
                  </div>
                </div>
              </CardHeader>
              <CardContent className="p-0 flex-1 min-h-[360px]">
                {selectedComplex?.cif_url ? (
                  <StructureViewer
                    key={`c-${selectedComplex.id}`}
                    cifUrl={`${selectedComplex.cif_url}`}
                    ligandBallAndStick
                    enablePocketToggle
                    pocketResidues={pocketResidues}
                    primaryLabel="Drugged Transcript"
                    onPlddtCalculated={setCanonPlddt}
                  />
                ) : (
                  <div className="h-full flex items-center justify-center text-sm text-muted-foreground">
                    No structure for this complex
                  </div>
                )}
              </CardContent>
              <div className="px-4 py-2 text-xs text-muted-foreground border-t border-border/50 shrink-0">
                Drug shown in ball &amp; stick. “Show binding pocket” highlights {pocketResidues.length} drug-contact
                residues (≤4.5 Å). Click the ligand to show the hydrogen bonds between the drug and protein.
              </div>
          </Card>

          {/* ── ROW 1 RIGHT: everything about the DRUG ────────────────────── */}
          <div className="flex flex-col gap-5">
            {/* Pooled counts + how splicing treats this gene's druggable pockets */}
            <DrugGeneOverviewPanel
              nComplexes={complexes.length}
              nDrugs={new Set(complexes.map(c => c.drug.pubchem_cid)).size}
              nIsoforms={siblings.length}
              nPockets={complexes.filter(c => (c.n_contact_residues ?? 0) > 0).length}
              fates={fates}
            />

            {/* Drug selector */}
            <Card>
              <CardHeader className="py-4"><CardTitle className="text-lg font-semibold">Drug ({sortedComplexes.length})</CardTitle></CardHeader>
              <CardContent className="max-h-[150px] overflow-y-auto space-y-1 pt-0">
                {sortedComplexes.map(c => (
                  <button key={c.id} onClick={() => setComplexId(c.id)}
                    className={cn("w-full text-left px-3 py-1.5 rounded-md text-sm flex items-center justify-between gap-2 transition-colors",
                      c.id === complexId ? "bg-emerald-500/15 text-emerald-300" : "hover:bg-muted/60")}>
                    <span className="truncate">{c.drug.name}</span>
                    <Badge className={cn("border text-[10px] shrink-0", gradeCls(c.holo_grade))}>{c.holo_grade ?? "–"}</Badge>
                  </button>
                ))}
              </CardContent>
            </Card>

            {/* Summary of the selected drug: skeleton + bound conformer + chemistry */}
            <DrugSummaryPanel
              drug={selectedComplex?.drug}
              complexCifUrl={selectedComplex?.cif_url}
            />
          </div>
        </div>

        {/* ── ROW 2: the selected isoform, and what happens to its pocket ── */}
        <div className="grid grid-cols-1 lg:grid-cols-[minmax(0,1.4fr)_minmax(0,1fr)] gap-5 items-stretch">
          <Card className="overflow-hidden flex flex-col">
              <CardHeader className="py-4">
                <div className="flex items-start justify-between gap-3">
                  <CardTitle className="text-lg font-semibold flex items-center gap-2 flex-wrap">
                    <span className={cn("w-2.5 h-2.5 rounded-full", fs.dot)} />
                    <span className="font-mono">{isoLabel(selectedIso)}</span>
                    {selectedIso?.transcript_name && (
                      <span className="text-muted-foreground/70 text-[11px] font-mono font-normal">
                        {selectedIso.transcript_id}
                      </span>
                    )}
                    {selectedIso?.is_mane_select && <Badge variant="secondary" className="text-[10px]">MANE</Badge>}
                    {fate && <Badge className={cn("border", fs.cls)}>{fs.label}</Badge>}
                  </CardTitle>
                  <div className="flex items-start gap-4 shrink-0">
                    <PlddtCircle label="Global pLDDT" value={isoPlddt} />
                    <PlddtCircle label="Pocket pLDDT" value={fate?.pocket_plddt} />
                  </div>
                </div>
              </CardHeader>
              <CardContent className="p-0 flex-1 min-h-[360px]">
                {isoCif ? (
                  <StructureViewer
                    key={`i-${selectedIso?.id}`}
                    cifUrl={`${isoCif}`}
                    canonicalCifUrl={selectedComplex?.cif_url}
                    sequence={selectedIso?.sequence}
                    canonicalSequence={canonicalSeq}
                    pocketResidues={pocketResidues}
                    pocketResiduesPrimary={retainedResidues}
                    enablePocketToggle={retainedResidues.length > 0}
                    ligandBallAndStick
                    focusPocketOnCanonical
                    canonicalLabel="canonical + drug"
                    primaryLabel={isoLabel(selectedIso)}
                    onPlddtCalculated={setIsoPlddt}
                  />
                ) : (
                  <div className="h-full flex items-center justify-center text-sm text-muted-foreground">
                    No structure for this isoform
                  </div>
                )}
              </CardContent>
              <div className="px-4 py-2 text-xs text-muted-foreground border-t border-border/50 space-y-1 shrink-0">
                <p>
                  {retainedResidues.length} pocket residues retained (aligned register). “Show canonical + drug”
                  ghosts the canonical, keeps the drug solid, and shows the pocket as ball-and-stick with H-bonds.
                </p>
                <p>
                  The two structures are superposed on the retained pocket Cα (exact-match anchors from a
                  splice-aware sequence alignment) so the binding site overlays cleanly; when ≤25% of the pocket
                  survives it falls back to a whole-chain fit.
                </p>
              </div>
          </Card>

          {/* ── ROW 2 RIGHT: fate of the selected isoform ─────────────────── */}
          <div className="flex flex-col gap-5">
            {/* Isoform selector with fate for the chosen drug */}
            <Card>
              <CardHeader className="py-4"><CardTitle className="text-lg font-semibold">Isoform fate · {selectedComplex?.drug.name}</CardTitle></CardHeader>
              <CardContent className="max-h-[190px] overflow-y-auto space-y-1 pt-0">
                {siblings.map(i => {
                  const f = fateByIso.get(i.id);
                  const st = FATE_STYLE[f?.fate ?? "no_call"] ?? FATE_STYLE.no_call;
                  return (
                    <button key={i.id} onClick={() => setIsoId(i.id)}
                      className={cn("w-full text-left px-3 py-1.5 rounded-md text-sm flex items-center justify-between gap-2 transition-colors",
                        i.id === isoId ? "bg-muted" : "hover:bg-muted/60")}>
                      <span className="truncate flex items-center gap-2">
                        <span className={cn("w-2 h-2 rounded-full shrink-0", st.dot)} />
                        <span className="font-mono">{isoLabel(i)}</span>
                        {i.transcript_name && (
                          <span className="text-muted-foreground/60 text-[10px] font-mono shrink-0">{i.transcript_id}</span>
                        )}
                        <span className="text-muted-foreground text-xs">{i.sequence?.length ?? "?"} aa</span>
                      </span>
                      {f && <span className="text-xs text-muted-foreground shrink-0">{Math.round((f.frac_retained ?? 0) * 100)}%</span>}
                    </button>
                  );
                })}
                {siblings.length === 0 && <p className="text-xs text-muted-foreground px-2">No sibling isoforms.</p>}
              </CardContent>
            </Card>

            {/* Pocket-fate analysis */}
            <Card>
              <CardHeader className="py-4"><CardTitle className="text-lg font-semibold flex items-center gap-2"><FlaskConical className="w-4 h-4 text-emerald-500" />Pocket fate</CardTitle></CardHeader>
              <CardContent className="pt-0 space-y-3 text-sm">
                {fate ? (
                  <>
                    <div className="flex items-center justify-between">
                      <span className="text-muted-foreground">Retention</span>
                      <span className="font-semibold">{fate.n_retained}/{fate.n_pocket_residues} residues · {Math.round((fate.frac_retained ?? 0) * 100)}%</span>
                    </div>
                    <div className="h-2 w-full rounded-full bg-muted overflow-hidden">
                      <div className="h-full bg-gradient-to-r from-amber-500 to-emerald-500"
                           style={{ width: `${(fate.frac_retained ?? 0) * 100}%` }} />
                    </div>
                    <div className="flex flex-wrap gap-2 pt-1">
                      <Badge className={cn("border", fs.cls)}>{fs.label}</Badge>
                      <Badge variant="secondary">{fate.mechanism_class ?? "—"}</Badge>
                      {!fate.decidable && (
                        <Badge className="bg-slate-500/15 text-slate-400 border-slate-500/30 border" title="isoform pocket-region pLDDT < 70">
                          not decidable (low pLDDT)
                        </Badge>
                      )}
                      {fate.pocket_plddt != null && <Badge variant="outline">pocket pLDDT {fate.pocket_plddt}</Badge>}
                    </div>
                    <p className="text-xs text-muted-foreground leading-relaxed pt-1">
                      Pocket defined from real drug contacts on the co-folded canonical structure, carried onto this
                      isoform by sequence alignment. A per-pocket pLDDT ≥ 70 gate marks low-confidence calls as
                      <span className="italic"> no&nbsp;call</span> rather than guessing.
                    </p>
                  </>
                ) : (
                  <p className="text-xs text-muted-foreground">Select an isoform to see its pocket fate.</p>
                )}
              </CardContent>
            </Card>

            {/* Placeholder: LigPlot+ (bachelor project — data pending) */}
            <Card className="border-dashed opacity-80">
              <CardHeader className="py-4"><CardTitle className="text-lg font-semibold flex items-center gap-2"><ImageOff className="w-4 h-4 text-muted-foreground" />2D interaction diagram (LigPlot+)</CardTitle></CardHeader>
              <CardContent className="pt-0">
                <div className="h-24 rounded-md border border-dashed border-border/60 flex items-center justify-center text-xs text-muted-foreground text-center px-4">
                  LigPlot+ binding-pocket diagram — pending from the interaction-analysis pipeline.
                </div>
              </CardContent>
            </Card>

            {/* Placeholder: ChimeraX contacts (data pending) */}
            <Card className="border-dashed opacity-80">
              <CardHeader className="py-4"><CardTitle className="text-lg font-semibold flex items-center gap-2"><Network className="w-4 h-4 text-muted-foreground" />H-bond / non-bonded contacts (ChimeraX)</CardTitle></CardHeader>
              <CardContent className="pt-0">
                <div className="h-20 rounded-md border border-dashed border-border/60 flex items-center justify-center text-xs text-muted-foreground text-center px-4">
                  Per-residue H-bond and non-bonded contact table — pending; residues will be cross-checked against each isoform.
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </div>
  );
}
