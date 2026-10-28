"use client";

import { ExternalLink, Box, Columns2, Dna, ArrowLeft, Scissors, Activity, Zap } from "lucide-react";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Badge } from "@/components/ui/badge";
import Link from "next/link";

export default function ProteinStructureDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4">
        <div className="flex items-center gap-3 text-primary">
          <Dna className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">Protein Structure</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          The 3D models, confidence signals, structural impacts, and alignment logic that project splicing events into physical space.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-8 text-foreground">
        <section className="space-y-4">
          <p className="leading-relaxed">
            Every isoform in SPLISOFORMS is anchored to a three-dimensional protein model predicted with{" "}
            <a href="https://github.com/google-deepmind/alphafold3" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">AlphaFold 3 (AF3)</a>.
            These models, together with their per-residue confidence signals, are the substrate on which every
            downstream panel — domains, structural impact, disorder, PTM accessibility, interfaces — is computed.
            Canonical and novel isoforms are all modelled under the same AF3 workflow, so their structures are
            directly comparable.
          </p>
        </section>

        {/* ── Confidence Metrics ─────────────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Model Confidence — pLDDT &amp; PAE</h2>
          <p className="leading-relaxed">
            Every predicted structure carries two AlphaFold confidence signals that gate the rest of the
            analysis stack — low-confidence regions are explicitly flagged rather than analysed as if reliable:
          </p>

          <div className="grid gap-4 md:grid-cols-2">
            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 space-y-2">
              <h3 className="font-bold text-foreground">pLDDT — per-residue confidence</h3>
              <p className="text-xs text-muted-foreground leading-relaxed">
                A score in <strong>[0, 100]</strong> stored in the B-factor column of the mmCIF; the
                per-residue value and the whole-chain mean are read from the Cα atom. Higher = better.
                We follow standard bands:
              </p>
              <ul className="text-xs text-muted-foreground list-disc list-inside space-y-1">
                <li><span className="font-mono text-[#0053D6]">≥ 90</span> — very high confidence (atomic detail trustworthy)</li>
                <li><span className="font-mono text-[#63C7F4]">70–90</span> — confident backbone, plausible side chains</li>
                <li><span className="font-mono text-[#FFD600]">50–70</span> — low confidence; backbone may be roughly right but side chains drift</li>
                <li><span className="font-mono text-[#FF8A00]">&lt; 50</span> — very low confidence; often disordered or modelled as a guess</li>
              </ul>
              <p className="text-xs text-muted-foreground leading-relaxed">
                The PTM, DSSP, and structural-landscape panels apply pLDDT-confidence
                gates so they mark low-confidence regions as
                <code className="font-mono mx-1">undetermined</code> /
                <code className="font-mono mx-1">low_confidence</code> rather than emitting structural state
                changes that ride on noise.
              </p>
            </div>

            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 space-y-2">
              <h3 className="font-bold text-foreground">PAE — predicted aligned error</h3>
              <p className="text-xs text-muted-foreground leading-relaxed">
                A residue × residue matrix where entry (i, j) is the model's expected position error of
                residue j when the structure is aligned on residue i. Low PAE between two residues =
                "their relative position is confident" = <strong>rigid pair</strong>.
              </p>
              <p className="text-xs text-muted-foreground leading-relaxed">
                We summarise the matrix as <strong>rigid blocks</strong> — contiguous regions whose internal
                mean PAE stays below 6 Å, meaning the model expects them to move as a single rigid unit. These
                blocks drive the per-isoform structural landscape view and flag which regions are confident
                enough for downstream impact assessment. Mean PAE and rigid blocks are computed for every
                one of the 124,734 AlphaFold 3 structures (131,554 isoforms in total; the remainder are
                reference proteins without a model).
              </p>
            </div>
          </div>
        </section>

        {/* ── Per-residue Structural Impact (DSSP) ───────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Scissors className="w-5 h-5 text-rose-500" />
            Structural Impact (DSSP)
          </h2>
          <p className="leading-relaxed">
            For every isoform structure we extract the N/CA/C/O backbone coordinates and run{" "}
            <a href="https://github.com/ShintaroMinami/PyDSSP" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">PyDSSP</a>
            {" "}(a vectorized, simplified reimplementation of DSSP) to assign per-residue secondary structure in
            three states — helix, strand, and coil. Pairing canonical and novel residues through the
            splice-aware aligner then lets us track the fold consequences of each splice event as a change of
            secondary-structure state. The resulting per-residue <code className="font-mono">structural_change</code> is one of:
          </p>

          <div className="grid gap-3">
            <div className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
              <Badge className="bg-emerald-500/15 text-emerald-400 hover:bg-emerald-500/15 text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap">same</Badge>
              <p className="text-sm leading-relaxed text-muted-foreground">DSSP assignment (Helix/Sheet/Coil) is unchanged between novel and canonical at this aligned position.</p>
            </div>
            <div className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
              <Badge className="bg-rose-500/15 text-rose-400 hover:bg-rose-500/15 text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap">coil_to_helix · helix_to_coil · …</Badge>
              <p className="text-sm leading-relaxed text-muted-foreground">Confident structural transitions: both canonical and novel pLDDT ≥ 70. Six pairwise transitions between Helix, Sheet, and Coil are tracked.</p>
            </div>
            <div className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
              <Badge className="bg-amber-500/15 text-amber-400 hover:bg-amber-500/15 text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap">low_confidence</Badge>
              <p className="text-sm leading-relaxed text-muted-foreground">
                At least one of canonical / novel pLDDT &lt; 70. Side-chain placement is unreliable there,
                so the structural transition direction cannot be trusted. 
              </p>
            </div>
            <div className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
              <Badge className="bg-blue-500/15 text-blue-400 hover:bg-blue-500/15 text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap">novel_residue</Badge>
              <p className="text-sm leading-relaxed text-muted-foreground">
                The novel residue has no canonical counterpart in the splice-aware alignment — splicing
                introduced new sequence at this position. 
              </p>
            </div>
          </div>
        </section>

        {/* ── Predicted disorder & Binding ───────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Activity className="w-5 h-5 text-amber-500" />
            Disorder &amp; Functional Binding (AIUPred)
          </h2>
          <p className="leading-relaxed">
            Splicing frequently rewires intrinsically disordered regions (IDRs), which act as flexible hubs for
            protein–protein interactions. For each isoform we predict per-residue intrinsic disorder and
            binding propensity with{" "}
            <a href="https://aiupred.elte.hu/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">AIUPred and AIUPred-binding</a>{" "}
            (<a href="https://doi.org/10.1093/nar/gkae385" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">Erdős &amp; Dosztányi, <em>Nucleic Acids Res.</em> 2024</a>),
            a deep-learning predictor that pairs an energy-estimation model with a transformer sequence
            embedding. Profiles exist for 130,501 of the 131,554 isoforms; the remainder are sequences
            shorter than 25 or longer than 4,000 residues, which the predictor does not accept.
          </p>
          <div className="grid gap-4 md:grid-cols-2">
            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 space-y-2">
              <h3 className="font-bold text-foreground flex items-center gap-2">
                 <Activity className="w-4 h-4 text-amber-500" />
                 AIUPred Disorder Propensity
              </h3>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Estimates the probability of a residue lacking a well-defined 3D structure, from the
                sequence alone.
                <br /><br />
                <strong>Scores ≥ 0.5</strong> mark residues that are highly likely to be flexible/disordered in vivo. Splicing within these regions often alters linear motifs without destabilizing a structured domain.
              </p>
            </div>

            <div className="rounded-xl border border-border/40 bg-muted/10 p-4 space-y-2">
              <h3 className="font-bold text-foreground flex items-center gap-2">
                 <Zap className="w-4 h-4 text-emerald-500" />
                 AIUPred-binding Functional Binding
              </h3>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Predicts protein binding regions located within intrinsically disordered sequences. AIUPred-binding scores
                the probability of an unstructured region undergoing a disorder-to-order transition upon binding a partner.
                <br /><br />
                <strong>Scores ≥ 0.5</strong> indicate a likely binding-coupled IDR. Splicing out these residues can cleanly ablate a specific protein interaction.
              </p>
            </div>
          </div>
        </section>

        {/* Structure Viewer */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Box className="w-5 h-5 text-indigo-400" />
            Structure Viewer
          </h2>
          <p className="text-muted-foreground">
            Every gene page includes an interactive <strong>Structure Viewer</strong>, powered by{" "}
            <a href="https://molstar.org/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">Mol*</a>,
            that overlays selected isoform structures onto the canonical reference.
          </p>
          <div className="grid gap-4 sm:grid-cols-2">
            <div className="p-4 rounded-lg bg-muted/20 border border-border/40 space-y-2">
              <h4 className="font-bold text-foreground text-sm">Superposition</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Isoforms are superposed onto the canonical using only the shared-exon residues from the same
                splice-aware alignment shown in the Sequence Alignment panel below — so novel exon sequence
                never skews the fit.
              </p>
            </div>
            <div className="p-4 rounded-lg bg-muted/20 border border-border/40 space-y-2">
              <h4 className="font-bold text-foreground text-sm">RMSD</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Root-Mean-Square Deviation (Å) is computed for each overlaid isoform and shown in the
                legend panel. Lower RMSD indicates a more structurally similar isoform.
              </p>
            </div>
            <div className="p-4 rounded-lg bg-muted/20 border border-border/40 space-y-2">
              <h4 className="font-bold text-foreground text-sm">Color Coding</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                The canonical reference is rendered in neutral slate-silver, and each overlaid isoform gets
                its own distinct colour, matched to its entry in the legend.
              </p>
            </div>
            <div className="p-4 rounded-lg bg-muted/20 border border-border/40 space-y-2">
              <h4 className="font-bold text-foreground text-sm">Multi-Select</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Several isoforms can be overlaid at once. Use the Mapped Transcripts table checkboxes
                (or Select All) to add or remove isoforms from the 3D view.
              </p>
            </div>
          </div>
        </section>

        {/* MSA */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Columns2 className="w-5 h-5 text-emerald-500" />
            Sequence Alignment Panel (MSA)
          </h2>
          <p className="text-muted-foreground">
            The gene page also displays a multi-sequence alignment (MSA) panel comparing all selected isoforms
            to the canonical. Alignment is performed by a custom <strong>splice-aware algorithm</strong> rather
            than standard Needleman-Wunsch, ensuring biologically correct coloring:
          </p>
          <div className="grid gap-4 sm:grid-cols-3">
            <div className="p-4 rounded-lg bg-emerald-500/5 border border-emerald-500/10 space-y-1">
              <h4 className="font-bold text-emerald-400 text-sm">Green columns</h4>
              <p className="text-xs text-muted-foreground">Novel exon insertions — sequence present only in the isoform, not in canonical.</p>
            </div>
            <div className="p-4 rounded-lg bg-rose-500/5 border border-rose-500/10 space-y-1">
              <h4 className="font-bold text-rose-400 text-sm">Pink dashes</h4>
              <p className="text-xs text-muted-foreground">Skipped canonical exons — regions present in canonical but absent in the isoform.</p>
            </div>
            <div className="p-4 rounded-lg bg-amber-500/5 border border-amber-500/10 space-y-1">
              <h4 className="font-bold text-amber-400 text-sm">Amber columns</h4>
              <p className="text-xs text-muted-foreground">Short substitutions (≤ 3 aa) — likely SNPs or minor junction differences, not splice events.</p>
            </div>
          </div>

          <Alert className="bg-indigo-500/5 border-indigo-500/20">
            <Dna className="h-4 w-4 text-indigo-400" />
            <AlertTitle className="text-indigo-400 font-bold">Algorithm: Splice-Aware k-mer Anchoring</AlertTitle>
            <AlertDescription className="text-muted-foreground text-xs leading-relaxed mt-1">
              The aligner builds a k-mer index on the canonical sequence (k=5), finds exact-match seeds in the
              isoform, extends them to maximal shared blocks, and chains them using a longest-weight increasing
              subsequence. Only genuinely shared exon sequence is aligned as identity; novel/alternative exon
              sequence appears as green insertions rather than misleading amber substitutions. The same aligner
              underpins the Structural Impact, PTM, and neoantigen panels.
            </AlertDescription>
          </Alert>
        </section>

        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "AlphaFold 3", note: "structure prediction", href: "https://github.com/google-deepmind/alphafold3" },
              { label: "Mol*", note: "interactive 3D viewer", href: "https://molstar.org/" },
              { label: "PyDSSP", note: "secondary-structure assignment", href: "https://github.com/ShintaroMinami/PyDSSP" },
              { label: "AIUPred / AIUPred-binding", note: "disorder & binding prediction", href: "https://aiupred.elte.hu/" },
            ].map((s) => (
              <li key={s.label}>
                <a href={s.href} target="_blank" rel="noopener noreferrer"
                   className="flex items-center gap-1.5 text-sm text-muted-foreground hover:text-primary transition-colors">
                  <ExternalLink className="w-3.5 h-3.5 shrink-0" />
                  <span className="font-semibold text-foreground">{s.label}</span>
                  <span className="text-xs">— {s.note}</span>
                </a>
              </li>
            ))}
          </ul>
        </section>

        {/* Footer nav */}
        <div className="pt-10 border-t border-border/20 flex items-center justify-between gap-4">
          <Link href="/documentation" className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2">
            <ArrowLeft className="w-4 h-4" /> Overview
          </Link>
          <div className="flex flex-col gap-1 items-end">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest">Next Section</span>
            <Link href="/documentation/functional-domains" className="text-lg font-bold text-primary hover:underline flex items-center gap-2">
              Functional Domains <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
