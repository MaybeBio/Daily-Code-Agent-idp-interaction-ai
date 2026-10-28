"use client";

import { useRouter } from "next/navigation";
import Link from "next/link";
import { ArrowLeft, Layers, Box, Columns2, Activity, Flag, GitCompare, BookOpen } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Card, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { SearchCatalog } from "@/components/SearchCatalog";

const features = [
  {
    icon: Box,
    color: "text-blue-500",
    title: "Multi-isoform overlay",
    desc: "Superpose every isoform onto the canonical reference in an interactive Mol* viewer, aligned on shared-exon residues only.",
  },
  {
    icon: Activity,
    color: "text-indigo-400",
    title: "Confidence-aware",
    desc: "Per-residue pLDDT colouring and PAE rigid-block decomposition gate what is and isn't trustworthy in each model.",
  },
  {
    icon: Columns2,
    color: "text-emerald-500",
    title: "Splice-aware alignment",
    desc: "A k-mer-anchored aligner colours novel insertions, skipped exons, and substitutions correctly — not as misleading mismatches.",
  },
];

const shortcuts = [
  { icon: Flag, href: "/flagged", label: "Flagged isoforms", desc: "Curated isoforms of interest" },
  { icon: GitCompare, href: "/compare", label: "Compare isoforms", desc: "Side-by-side structural comparison" },
  { icon: BookOpen, href: "/documentation/protein-structure", label: "How it works", desc: "Protein Structure documentation" },
];

export default function StructurePage() {
  const router = useRouter();

  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-5xl mx-auto p-6 md:p-10">
        <header className="flex items-center gap-4 mb-10">
          <Button variant="ghost" size="icon" onClick={() => router.push("/")} className="hover:bg-primary/10">
            <ArrowLeft className="w-5 h-5" />
          </Button>
          <div className="flex items-center gap-3">
            <Layers className="w-7 h-7 text-primary" />
            <div>
              <h1 className="text-4xl font-bold tracking-tight bg-clip-text text-transparent bg-gradient-to-r from-primary via-foreground/80 to-secondary">
                Structure Explorer
              </h1>
              <p className="text-sm text-muted-foreground mt-1">
                Search a gene to open its isoform structures in the 3D viewer.
              </p>
            </div>
          </div>
        </header>

        {/* Search */}
        <div className="mb-12">
          <SearchCatalog />
        </div>

        {/* What the viewer offers */}
        <div className="grid gap-6 sm:grid-cols-3 mb-12">
          {features.map((f) => (
            <Card key={f.title} className="bg-card/40 border-border/40">
              <CardHeader>
                <f.icon className={`w-8 h-8 ${f.color} mb-2`} />
                <CardTitle className="text-base">{f.title}</CardTitle>
                <CardDescription className="text-sm leading-relaxed pt-1">{f.desc}</CardDescription>
              </CardHeader>
            </Card>
          ))}
        </div>

        {/* Shortcuts */}
        <div className="grid gap-4 sm:grid-cols-3">
          {shortcuts.map((s) => (
            <Link
              key={s.href}
              href={s.href}
              className="group flex items-center gap-3 p-4 rounded-xl border border-border/40 bg-muted/20 hover:bg-muted/40 hover:border-primary/30 transition-colors"
            >
              <div className="w-9 h-9 rounded-lg bg-primary/10 flex items-center justify-center shrink-0">
                <s.icon className="w-4 h-4 text-primary" />
              </div>
              <div className="min-w-0">
                <p className="text-sm font-semibold text-foreground group-hover:text-primary transition-colors">{s.label}</p>
                <p className="text-xs text-muted-foreground truncate">{s.desc}</p>
              </div>
            </Link>
          ))}
        </div>
      </div>
    </div>
  );
}
