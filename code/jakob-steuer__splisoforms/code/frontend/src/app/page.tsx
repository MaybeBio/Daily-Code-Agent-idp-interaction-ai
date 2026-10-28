"use client";

import { Dna, Scissors, Shield, OctagonAlert, Target, Layers, Puzzle, Pill, Network } from "lucide-react";
import { Card, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { SearchCatalog } from "@/components/SearchCatalog";
import Link from "next/link";

import { StatsBanner } from "@/components/StatsBanner";

// Landing feature cards — a 3×3 tour of the knowledgebase. Explore lenses (the interactive
// tools) come first, deep-dives link into the docs. `hover`/`icon` are full literal Tailwind
// class strings, NOT interpolated: the JIT compiler only keeps classes it sees verbatim in
// source, so a template like `hover:border-${c}` would be purged and never render.
const FEATURES = [
  {
    href: "/structure", icon: Dna, iconCls: "text-primary", hover: "hover:border-primary/50",
    title: "Structure Explorer", titleHover: "group-hover:text-primary",
    desc: "AlphaFold 3 models for every long-read isoform — per-residue pLDDT confidence and interactive Mol* visualization.",
  },
  {
    href: "/explore/impact", icon: Scissors, iconCls: "text-secondary", hover: "hover:border-secondary/50",
    title: "Splicing Impact", titleHover: "group-hover:text-secondary",
    desc: "How splice variants reshape the fold: secondary-structure shifts (PyDSSP), solvent-accessibility changes, and AIUPred disorder — scored per residue against the canonical.",
  },
  {
    href: "/explore/drugs", icon: Pill, iconCls: "text-teal-500", hover: "hover:border-teal-500/50",
    title: "Drug Pockets", titleHover: "group-hover:text-teal-400",
    desc: "AlphaFold 3 co-folded drug–protein complexes. Track whether each binding pocket is preserved, disrupted, or lost in every splice isoform.",
  },
  {
    href: "/explore/neoepitopes", icon: Shield, iconCls: "text-purple-500", hover: "hover:border-purple-500/50",
    title: "Neoepitope Atlas", titleHover: "group-hover:text-purple-400",
    desc: "Splice-junction neoepitopes scored by a three-tool MHC ensemble (MHCflurry · NetMHCpan-4.2 · BigMHC), ranked into P1–P4 priority tiers, and flagged when the peptide is absent from the whole reference proteome.",
  },
  {
    href: "/documentation/functional-domains", icon: Puzzle, iconCls: "text-emerald-500", hover: "hover:border-emerald-500/50",
    title: "Functional Domains", titleHover: "group-hover:text-emerald-400",
    desc: "Pfam HMM mapping and TED / Chainsaw structural segmentation, classifying every domain as gained, lost, or truncated relative to the gene's canonical isoform.",
  },
  {
    href: "/documentation/interactions-expression", icon: Network, iconCls: "text-sky-500", hover: "hover:border-sky-500/50",
    title: "Protein Interactions", titleHover: "group-hover:text-sky-400",
    desc: "STRING seven-channel interaction evidence, plus AlphaFold-derived homodimer and heterodimer interface disruption — quantifying which contacts a splice event breaks.",
  },
  {
    href: "/documentation/nmd-transcript-fate", icon: OctagonAlert, iconCls: "text-amber-500", hover: "hover:border-amber-500/50",
    title: "NMD & Transcript Fate", titleHover: "group-hover:text-amber-400",
    desc: "Nonsense-mediated decay prediction via the 50-nucleotide rule on annotated CDS coordinates, honouring GENCODE's NMD tags, with truncation severity for every isoform.",
  },
  {
    href: "/documentation/ptm-accessibility", icon: Target, iconCls: "text-rose-500", hover: "hover:border-rose-500/50",
    title: "PTM Accessibility", titleHover: "group-hover:text-rose-400",
    desc: "How splicing reshapes post-translational modification sites — ΔSASA exposure shifts for phosphorylation, ubiquitination, and acetylation, gated by structural confidence.",
  },
  {
    href: "/documentation", icon: Layers, iconCls: "text-indigo-500", hover: "hover:border-indigo-500/50",
    title: "Multi-Tissue Coverage", titleHover: "group-hover:text-indigo-400",
    desc: "Long-read transcriptomes from ccRCC organoids and breast-cancer tissue, anchored to the GENCODE v47 / GRCh38 reference proteome.",
  },
];

export default function Home() {

  return (
    <div className="min-h-screen bg-background text-foreground flex flex-col items-center selection:bg-primary/30 relative overflow-hidden">

      {/* Dark-mode: protein hero background.
          opacity dims the saturated structure without greying-out the solid
          black bg (black pixel * 0.5 over black page bg = still black).
          saturate(0.7) takes the colour intensity down a notch. */}
      <div
        className="absolute inset-0 z-0 pointer-events-none hidden dark:block"
        style={{
          backgroundImage: 'url("/background-dark.webp")',
          backgroundSize: 'cover',
          backgroundPosition: 'center',
          backgroundRepeat: 'no-repeat',
          opacity: 0.55,
          filter: 'saturate(0.7)',
        }}
      />
      <div
        className="absolute inset-0 z-0 pointer-events-none hidden dark:block"
        style={{
          background: 'radial-gradient(circle at center, transparent 0%, rgba(var(--background-rgb), 0.4) 40%, var(--background) 85%)'
        }}
      />

      {/* Light-mode: protein hero background.
          More aggressive dim/desaturate — saturated structure pops harder
          against a light page bg, so we push it further back. */}
      <div
        className="absolute inset-0 z-0 pointer-events-none dark:hidden"
        style={{
          backgroundImage: 'url("/background-light.webp")',
          backgroundSize: 'cover',
          backgroundPosition: 'center',
          backgroundRepeat: 'no-repeat',
          opacity: 0.4,
          filter: 'saturate(0.6)',
        }}
      />
      <div
        className="absolute inset-0 z-0 pointer-events-none dark:hidden"
        style={{
          background: 'radial-gradient(circle at center, transparent 0%, rgba(var(--background-rgb), 0.4) 40%, var(--background) 85%)'
        }}
      />
      
      {/* Navigation moved to layout.tsx */}

      {/* Hero Section */}
      <main className="w-full max-w-7xl mx-auto px-6 py-24 flex flex-col items-center text-center flex-1 relative z-10">
        
        {/* Release pill — flat, no blur-glow halo (the glow was the worst
            offender for the production paint-thrash that caused the
            top-of-page flicker reported on the FHNW server). */}
        <div className="relative border border-primary/30 bg-background/50 backdrop-blur-sm px-4 py-1.5 rounded-full text-sm font-medium text-primary mb-8 inline-flex items-center gap-2">
          <span className="relative flex items-center justify-center w-2 h-2">
            <span className="absolute inset-0 rounded-full bg-primary animate-ping opacity-60" />
            <span className="relative w-1.5 h-1.5 rounded-full bg-primary" />
          </span>
          v1.3.0 · Live
        </div>

        <h1 className="text-6xl md:text-8xl font-extrabold tracking-tight mb-8 max-w-4xl leading-tight">
          <span className="block mb-4 text-transparent bg-clip-text bg-gradient-to-r from-primary via-primary/80 to-secondary drop-shadow-sm">
            SPLISOFORMS
          </span>
          <span className="block text-3xl md:text-5xl text-foreground font-medium tracking-wide">
            Splicing Proteomics Knowledgebase
          </span>
        </h1>

        <p className="text-lg md:text-xl text-muted-foreground max-w-2xl mb-12 leading-relaxed">
          Translating alternative splicing into structural, immunological, and regulatory insights — across cancer types and tissues.
        </p>

        <StatsBanner />

        <div className="w-full max-w-2xl flex flex-col sm:flex-row gap-4 mb-24 relative z-10">
          <SearchCatalog />
        </div>

        {/* Feature Cards — 3×3 tour of the knowledgebase */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 w-full text-left relative z-0">
          <div className="absolute -inset-x-20 top-1/2 -translate-y-1/2 h-64 bg-secondary/10 blur-[100px] -z-10 rounded-[100%]" />

          {FEATURES.map((f) => {
            const Icon = f.icon;
            return (
              <Link key={f.href} href={f.href} className="block group">
                <Card className={`h-full bg-card/40 backdrop-blur-sm border-border/50 ${f.hover} hover:-translate-y-1 transition-all`}>
                  <CardHeader>
                    <Icon className={`w-10 h-10 ${f.iconCls} mb-4`} />
                    <CardTitle className={`${f.titleHover} transition-colors`}>{f.title}</CardTitle>
                    <CardDescription>{f.desc}</CardDescription>
                  </CardHeader>
                </Card>
              </Link>
            );
          })}
        </div>

      </main>
      
    </div>
  );
}
