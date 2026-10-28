"use client";

import { Book, Dna, Puzzle, OctagonAlert, Shield, Target, Compass, Scissors, Network, ArrowRight, Database, FlaskConical, Microscope, Pill } from "lucide-react";
import { Button } from "@/components/ui/button";
import Link from "next/link";
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from "@/components/ui/card";

const modules = [
  {
    name: "Protein Structure & Impact",
    description: "How splicing reshapes the 3D fold, compared residue by residue against the canonical protein.",
    href: "/documentation/protein-structure",
    icon: Dna,
    color: "text-blue-500",
    bg: "bg-blue-500/10"
  },
  {
    name: "Functional Domains",
    description: "Which protein domains a splice event keeps, truncates, or removes.",
    href: "/documentation/functional-domains",
    icon: Puzzle,
    color: "text-emerald-500",
    bg: "bg-emerald-500/10"
  },

  {
    name: "NMD & Transcript Fate",
    description: "Whether a transcript is likely degraded by nonsense-mediated decay or translated into protein.",
    href: "/documentation/nmd-transcript-fate",
    icon: OctagonAlert,
    color: "text-amber-500",
    bg: "bg-amber-500/10"
  },
  {
    name: "Neoantigens",
    description: "Splicing-derived peptides scored for MHC-I presentation and immunogenicity, ranked into P1–P4 priority tiers.",
    href: "/documentation/neoantigens",
    icon: Shield,
    color: "text-purple-500",
    bg: "bg-purple-500/10"
  },
  {
    name: "PTM Accessibility",
    description: "How splicing exposes, buries, or removes post-translational modification sites.",
    href: "/documentation/ptm-accessibility",
    icon: Target,
    color: "text-rose-500",
    bg: "bg-rose-500/10"
  },
  {
    name: "Interactions & Expression",
    description: "Protein interactions and dimer interfaces disrupted by splicing, plus tissue and cell-type expression.",
    href: "/documentation/interactions-expression",
    icon: Network,
    color: "text-sky-500",
    bg: "bg-sky-500/10"
  },
  {
    name: "Discovery & Browsing",
    description: "Rank and explore the most consequential isoforms through the impact score and interactive browsers.",
    href: "/documentation/discovery",
    icon: Compass,
    color: "text-indigo-500",
    bg: "bg-indigo-500/10"
  },
  {
    name: "Drug Interactions",
    description: "Whether each drug's binding pocket survives splicing — preserved, disrupted, or lost.",
    href: "/documentation/drug-interactions",
    icon: Pill,
    color: "text-teal-500",
    bg: "bg-teal-500/10"
  }
];

export default function DocumentationPage() {
  return (
    <div className="space-y-12">
      {/* Hero Section */}
      <div className="space-y-4 max-w-2xl px-3">
        <h1 className="text-4xl font-extrabold tracking-tight text-foreground sm:text-5xl">
          SPLISOFORMS Documentation
        </h1>
        <p className="text-lg text-muted-foreground leading-relaxed">
          The SPLISOFORMS splicing proteomics knowledgebase provides a unified view of how alternative splicing events reshape the protein landscape.
          Use these guides to understand the biological logic and technical algorithms powering our analytical discovery engine.
        </p>
      </div>

      <div className="grid gap-6 sm:grid-cols-2">
        {modules.map((module) => (
          <Link key={module.name} href={module.href} className="group transition-all duration-200 hover:-translate-y-1">
            <Card className="h-full border-border/40 bg-muted/20 hover:bg-muted/40 transition-colors">
              <CardHeader>
                <div className={`w-10 h-10 rounded-lg ${module.bg} flex items-center justify-center mb-2`}>
                  <module.icon className={`w-5 h-5 ${module.color}`} />
                </div>
                <CardTitle className="group-hover:text-primary transition-colors">{module.name}</CardTitle>
                <CardDescription className="text-sm leading-relaxed pt-1">
                  {module.description}
                </CardDescription>
              </CardHeader>
              <CardContent>
                <div className="flex items-center text-xs font-bold text-primary opacity-0 group-hover:opacity-100 transition-opacity">
                  LEARN MORE <ArrowRight className="ml-1 w-3 h-3" />
                </div>
              </CardContent>
            </Card>
          </Link>
        ))}
      </div>

      {/* Datasets Section */}
      <div className="space-y-6">
        <div className="space-y-1">
          <h2 className="text-2xl font-bold tracking-tight">Datasets</h2>
          <p className="text-muted-foreground text-sm leading-relaxed">
            The foundation of SPLISOFORMS is the GRCh38 reference proteome, supplemented by additional 
            long-read transcriptomes from multiple experimental contexts. All models are processed through
            the same structural prediction and analysis pipeline.
          </p>
        </div>

        <div className="grid gap-4 sm:grid-cols-1">
          {/* Reference */}
          <div className="border border-border/40 bg-muted/20 rounded-xl p-5 space-y-3">
            <div className="flex items-center gap-3">
              <div className="w-9 h-9 rounded-lg bg-emerald-500/10 flex items-center justify-center flex-shrink-0">
                <Database className="w-4 h-4 text-emerald-500" />
              </div>
              <div>
                <p className="font-semibold text-sm text-foreground">Reference Proteome</p>
                <p className="text-[10px] uppercase tracking-widest text-muted-foreground font-bold">GRCh38 · GENCODE v47</p>
              </div>
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed">
              All protein-coding transcripts of GENCODE v47 on GRCh38, with structures from the
              AlphaFold Database and the MANE Select transcript designated as each gene's canonical
              isoform. This is the baseline for every comparative calculation — structural impact,
              domain loss, PTM accessibility, and neoepitope novelty.
            </p>
            <div className="flex flex-wrap gap-1.5 pt-1">
              {["GENCODE v47", "MANE Select", "GRCh38"].map(t => (
                <span key={t} className="px-2 py-0.5 bg-emerald-500/10 text-emerald-400 rounded text-[10px] font-medium">{t}</span>
              ))}
            </div>
          </div>
        </div>

        <div className="space-y-4 pt-4 border-t border-border/40">
          <div className="space-y-1">
            <h3 className="text-xl font-bold tracking-tight">Long-Read Datasets</h3>
            <p className="text-muted-foreground text-sm leading-relaxed">
              Experimental datasets containing both annotated transcripts and novel isoforms discovered through long-read sequencing.
            </p>
          </div>
          <div className="grid gap-4 sm:grid-cols-1 md:grid-cols-2">
            {/* ccRCC */}
            <div className="border border-border/40 bg-muted/20 rounded-xl p-5 space-y-3">
              <div className="flex items-center gap-3">
                <div className="w-9 h-9 rounded-lg bg-blue-500/10 flex items-center justify-center flex-shrink-0">
                  <FlaskConical className="w-4 h-4 text-blue-500" />
                </div>
                <div>
                  <p className="font-semibold text-sm text-foreground">ccRCC Organoids</p>
                  <p className="text-[10px] uppercase tracking-widest text-muted-foreground font-bold">Karakulak et al. 2025</p>
                </div>
              </div>
              <p className="text-xs text-muted-foreground leading-relaxed">
                PacBio MAS-seq long-read sequencing of patient-derived clear cell renal cell
                carcinoma organoids (tumour organoid lines with a matched normal control).
                Isoform models derived from SQANTI3-classified assemblies. <a href="https://doi.org/10.1101/gr.279345.124" target="_blank" rel="noopener noreferrer" className="text-blue-500 hover:underline">View Publication</a>
              </p>
              <div className="flex flex-wrap gap-1.5 pt-1">
                {["Organoids", "Tumour + Normal", "SQANTI3", "PacBio MAS-seq"].map(t => (
                  <span key={t} className="px-2 py-0.5 bg-blue-500/10 text-blue-400 rounded text-[10px] font-medium">{t}</span>
                ))}
              </div>
            </div>

            {/* Breast Cancer */}
            <div className="border border-border/40 bg-muted/20 rounded-xl p-5 space-y-3">
              <div className="flex items-center gap-3">
                <div className="w-9 h-9 rounded-lg bg-rose-500/10 flex items-center justify-center flex-shrink-0">
                  <Microscope className="w-4 h-4 text-rose-500" />
                </div>
                <div>
                  <p className="font-semibold text-sm text-foreground">Breast Cancer</p>
                  <p className="text-[10px] uppercase tracking-widest text-muted-foreground font-bold">Head et al. 2026</p>
                </div>
              </div>
              <p className="text-xs text-muted-foreground leading-relaxed">
                Nanopore long-read ESPRESSO assemblies from healthy breast tissue, primary breast tumour,
                and cultured fibroblasts. Includes both annotated GENCODE transcripts and novel isoforms
                not present in reference databases. <a href="https://www.biorxiv.org/content/10.64898/2026.03.22.713514v2" target="_blank" rel="noopener noreferrer" className="text-rose-500 hover:underline">View Publication</a>
              </p>
              <div className="flex flex-wrap gap-1.5 pt-1">
                {["Healthy Breast", "Tumour", "Fibroblasts", "Nanopore · ESPRESSO"].map(t => (
                  <span key={t} className="px-2 py-0.5 bg-rose-500/10 text-rose-400 rounded text-[10px] font-medium">{t}</span>
                ))}
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* System Overview */}
      <div className="bg-muted/30 border border-border/40 rounded-xl p-8 space-y-4">
        <div className="flex items-center gap-3 mb-2">
           <Book className="w-5 h-5 text-primary" />
           <h3 className="font-bold text-lg">System Overview</h3>
        </div>
        <div className="prose prose-sm prose-invert max-w-none text-muted-foreground leading-relaxed space-y-4">
          <p>
            SPLISOFORMS is a splicing proteomics knowledgebase that transforms genomic transcriptomic data
            into three-dimensional conformational, immunological, and regulatory insights.
            By integrating structural predictions from custom <strong>AlphaFold 3</strong> models
            across multiple tissue contexts, the platform identifies high-impact perturbations that affect
            protein stability, functional domains, neoantigen potential, and post-translational regulation.
          </p>
          <p>
            Every protein-coding transcript and novel isoform is modelled with the same
            <strong> AlphaFold 3</strong> workflow for cross-comparable structures, on the{" "}
            <strong>GENCODE v47</strong> reference. Identical proteins from different transcripts or
            datasets are collapsed to one isoform by protein sequence hash; when a MANE Select
            transcript is among them, it is the one kept.
          </p>
          <p>
            <strong>One canonical per gene.</strong> Each gene has exactly one designated canonical
            isoform: its <strong>MANE Select</strong> transcript where one exists, otherwise the
            GENCODE/Ensembl canonical reference transcript, otherwise the longest reference protein.
            Genes with no GENCODE protein at all (lncRNA hosts, novel loci — 2,523 genes) take their
            longest novel ORF. Every canonical-relative analysis — domain loss, PTM accessibility,
            NMD truncation, neoantigen novelty, dimer interfaces — is measured against that
            isoform's sequence, and against nothing else.
          </p>
          <div className="not-prose rounded-xl border border-border/40 bg-background/40 p-4">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              Release snapshot · October 2026 (v2)
            </p>
            <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
              {[
                ["131,554", "isoforms"],
                ["22,075", "genes"],
                ["124,734", "AlphaFold 3 structures"],
                ["19,068", "MANE Select isoforms"],
              ].map(([v, l]) => (
                <div key={l}>
                  <div className="text-xl font-bold font-mono text-foreground">{v}</div>
                  <div className="text-xs text-muted-foreground">{l}</div>
                </div>
              ))}
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              Three datasets: GENCODE v47 reference (AlphaFold Database models), ccRCC organoids
              (PacBio), breast cancer (Nanopore). Per-layer counts are given on each module page.
            </p>
          </div>
          <p>
            <strong>October 2026 methods revision.</strong> Following an independent review of the
            pipeline, the database was rebuilt with the fixes described on the module pages: MANE-first
            resolution of identical proteins, clan-de-duplicated Pfam hits with 1-based coordinates, the
            full human TED set, the NMD rule applied only where a stop codon is annotated (GENCODE
            tags honoured), PTM sites checked against the canonical sequence, a complete three-tool
            neoantigen table with a proteome-uniqueness flag, per-sample expression joined within
            the right sample, and the repaired homodimer layer. Each module page states what changed.
          </p>
          <p>
            Live database counts are shown on the <Link href="/explore/impact" className="text-primary hover:underline">Isoform Explorer</Link>,
            {" "}<Link href="/explore/neoepitopes" className="text-primary hover:underline">Neoepitope Atlas</Link>, and
            {" "}<Link href="/explore/drugs" className="text-primary hover:underline">Drug Explorer</Link>.
          </p>
          <p>
            SPLISOFORMS is developed in the open. Questions, corrections, and feature requests are
            handled through GitHub — please{" "}
            <a href="https://github.com/jakob-steuer/splisoforms/issues" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">open an issue</a>{" "}
            on the <a href="https://github.com/jakob-steuer/splisoforms" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">source repository</a>.
          </p>
          <p>
            The resource and its methods are described in the preprint: Steuer J &amp; Kahraman A,{" "}
            <em>SPLISOFORMS: a Structure-Resolved Knowledge Base of Alternative Splicing Isoforms</em>,
            bioRxiv 2026,{" "}
            <a href="https://doi.org/10.64898/2026.08.06.743218" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">doi:10.64898/2026.08.06.743218</a>.
            Please cite it when you use SPLISOFORMS data.
          </p>
        </div>
      </div>
    </div>
  );
}
