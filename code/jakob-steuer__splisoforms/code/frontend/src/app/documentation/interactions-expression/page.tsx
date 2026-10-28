"use client";

import { Network, ExternalLink, ArrowLeft } from "lucide-react";
import Link from "next/link";

export default function InteractionsExpressionDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4">
        <div className="flex items-center gap-3 text-primary">
          <Network className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">Interactions &amp; Expression</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          Protein interactions, dimer interfaces, cell-type expression, and per-sample isoform abundance.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-10 text-foreground">
        <section className="space-y-4">
          <p className="leading-relaxed">
            Few proteins act alone. These panels place each isoform in its interaction and
            expression context — the dimer interfaces splicing can break, the wider network of
            partners a gene participates in, and where the gene is expressed across cell types and
            tissues.
          </p>
        </section>

        {/* ── Dimer interface disruption (homo + hetero) ─────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Dimer Interface Disruption</h2>
          <p className="leading-relaxed">
            Many proteins function as dimers — bound either to a second copy of themselves
            (<strong>homodimer</strong>) or to a distinct partner protein (<strong>heterodimer</strong>).
            A splice event that removes or reshapes the residues forming that binding interface can
            weaken or abolish the interaction, even when the rest of the protein looks intact.
          </p>

          <h3 className="text-lg font-bold text-foreground">Where the interface comes from</h3>
          <p className="leading-relaxed">
            Interface residues are read from the <strong>AlphaFold Database</strong> multimer models
            of each complex. AFDB's{" "}
            <a href="https://alphafold.ebi.ac.uk/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">summary endpoint</a>{" "}
            lists every structure for a protein tagged with its oligomeric state (monomer, homodimer,
            or heterodimer); for each dimer, the interface is the set of canonical residues with any
            heavy atom within <strong>5 Å</strong> of the partner chain in that model.
          </p>

          <h3 className="text-lg font-bold text-foreground">How residues are scored</h3>
          <p className="leading-relaxed">
            Each canonical interface residue is projected onto the isoform through a sequence
            alignment of the two proteins (Python's <code className="font-mono text-xs">SequenceMatcher</code>,
            which handles exon-skip deletions well), then sorted into one of three states:
          </p>
          <div className="grid gap-3 sm:grid-cols-3">
            <div className="p-4 rounded-lg border border-emerald-500/20 bg-emerald-500/5">
              <h4 className="font-bold text-emerald-400 text-sm mb-1">Preserved</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">Present in the isoform with essentially unchanged solvent accessibility — the contact is retained.</p>
            </div>
            <div className="p-4 rounded-lg border border-amber-500/20 bg-amber-500/5">
              <h4 className="font-bold text-amber-400 text-sm mb-1">Perturbed</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">A different amino acid sits at the aligned position, or the residue is kept but its solvent accessibility (SASA) shifted by more than 20 % — the contact may be compromised.</p>
            </div>
            <div className="p-4 rounded-lg border border-rose-500/20 bg-rose-500/5">
              <h4 className="font-bold text-rose-400 text-sm mb-1">Absent</h4>
              <p className="text-xs text-muted-foreground leading-relaxed">Mapped to a gap in the isoform — the residue was spliced out and the contact is physically gone.</p>
            </div>
          </div>
          <p className="leading-relaxed">
            The <strong>disruption score</strong> summarises this as{" "}
            <span className="font-mono text-sm">(perturbed + absent) ÷ total interface residues</span>,
            running from <strong>0</strong> (interface fully preserved) to <strong>1</strong> (fully
            disrupted). The absent and perturbed counts are stored separately, so wholesale deletion
            of an interface can be told apart from subtler conformational disruption.
          </p>
          <p className="leading-relaxed">
            Homodimer and heterodimer results are reported separately, and each heterodimer names its
            partner protein — so an isoform can be read as keeping some partnerships while losing
            others. On the isoform page, the interface residues can be projected straight into the{" "}
            <strong>Mol* structure viewer</strong> (&ldquo;Visualise interface&rdquo;) to see exactly
            which contacts survive the splice event.
          </p>
          <div className="rounded-xl border border-border/40 bg-muted/10 p-4 not-prose">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              Dimer interfaces · October 2026 release
            </p>
            <div className="grid grid-cols-2 gap-3">
              <div>
                <div className="text-xl font-bold font-mono text-foreground">12,006</div>
                <div className="text-xs text-muted-foreground">homodimer isoform rows</div>
              </div>
              <div>
                <div className="text-xl font-bold font-mono text-foreground">79,628</div>
                <div className="text-xs text-muted-foreground">heterodimer isoform rows</div>
              </div>
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              Rows include each gene's canonical compared against itself (score ≈ 0). The homodimer
              layer quadrupled in October 2026: a bug had skipped every gene whose first AlphaFold
              Database model was a homodimer.
            </p>
          </div>
        </section>

        {/* ── STRING interactions ────────────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">STRING Interaction Network</h2>
          <p className="leading-relaxed">
            Beyond direct dimers, each gene's highest-confidence partners are drawn from{" "}
            <a href="https://string-db.org/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">STRING</a>{" "}
            (queried by the gene's UniProt accession, direct partners only). Partners are kept at
            STRING's <strong>medium-confidence</strong> cut-off, combined score ≥ 0.4 (shown as
            ≥ 400 on STRING's 0–1000 scale), and the score is stored on each row so the panel can be
            narrowed to high confidence (≥ 0.7). STRING's combined score is built from seven
            independent evidence channels:
          </p>
          <div className="grid gap-3 sm:grid-cols-2">
            {[
              ["Experiments", "Direct biochemical / biophysical assays."],
              ["Database", "Curated pathway & complex memberships."],
              ["Co-expression", "Correlated expression across conditions."],
              ["Text-mining", "Co-mention in the literature."],
              ["Phylogenetic profile", "Co-occurrence across genomes."],
              ["Gene fusion", "Fused orthologs in other species."],
              ["Neighborhood", "Conserved genomic proximity (prokaryotes)."],
            ].map(([k, v]) => (
              <div key={k} className="p-3 rounded-lg bg-muted/20 border border-border/40">
                <span className="font-bold text-sky-400 text-sm">{k}</span>
                <p className="text-xs text-muted-foreground leading-relaxed">{v}</p>
              </div>
            ))}
          </div>
          <div className="rounded-xl bg-amber-500/5 border border-amber-500/20 p-4 text-sm">
            <p className="text-foreground/85 leading-relaxed">
              <strong>Note:</strong> STRING evidence is <strong>gene-level</strong> and does not
              distinguish isoform-specific interactions. A high combined score may reflect
              co-expression or literature co-mention rather than direct physical binding. The table
              stores one row per isoform of the gene for fast panel lookups, so its 1,177,680 rows
              (October 2026) correspond to <strong>176,930 gene–partner pairs</strong>, 126,080 of
              them at high confidence (≥ 0.7). Count pairs, not rows. 875 genes whose UniProt
              accession is an unreviewed TrEMBL entry have no STRING match.
            </p>
          </div>
        </section>

        {/* ── Per-sample isoform expression (ccRCC long-read) ─────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Per-Sample Isoform Expression (ccRCC)</h2>
          <p className="leading-relaxed">
            The ccRCC organoid dataset (Karakulak et al. 2025) is the one source with quantitative,
            isoform-level counts. For each of its five samples (one normal, four tumour organoid
            lines) the SQANTI3 output gives every transcript's number of full-length (FL) reads;
            we store the FL count, the counts-per-million (CPM) within the sample, and the isoform's
            share of its gene's reads.
          </p>
          <div className="rounded-xl bg-amber-500/5 border border-amber-500/20 p-4 text-sm">
            <p className="text-foreground/85 leading-relaxed">
              <strong>Why the join is by protein, not by transcript id.</strong> PacBio transcript ids
              (<code className="font-mono text-xs">PB.*</code>) are numbered independently in each
              sample's SQANTI run, so the same id can name unrelated transcripts in two samples. Each
              isoform is therefore matched by its protein sequence hash to the (id, sample) pairs it
              was assembled from, and looked up only in that sample's file. GENCODE reference proteins
              that were also observed in the organoids are covered the same way. Breast-cancer isoforms
              carry presence flags per tissue but no per-sample counts.
            </p>
          </div>
          <div className="rounded-xl border border-border/40 bg-muted/10 p-4 not-prose">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              Isoform expression · October 2026 release
            </p>
            <div className="grid grid-cols-3 gap-3">
              <div>
                <div className="text-xl font-bold font-mono text-foreground">39,768</div>
                <div className="text-xs text-muted-foreground">isoform × sample rows</div>
              </div>
              <div>
                <div className="text-xl font-bold font-mono text-foreground">34,846</div>
                <div className="text-xs text-muted-foreground">isoforms with counts</div>
              </div>
              <div>
                <div className="text-xl font-bold font-mono text-foreground">3,126</div>
                <div className="text-xs text-muted-foreground">seen in more than one sample</div>
              </div>
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              Before October 2026 ids were joined across samples; about 4 % of rows then carried
              another transcript's counts and almost every isoform appeared in exactly one sample.
            </p>
          </div>
        </section>

        {/* ── Cell-type expression (HPA) + tissue map ────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Expression — Cell Types &amp; Tissues</h2>
          <p className="leading-relaxed">
            Each gene carries a single-cell expression profile drawn from the{" "}
            <a href="https://www.proteinatlas.org/" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">Human Protein Atlas</a>{" "}
            (HPA) single-cell atlas of <strong>154 cell types</strong>. We store HPA's
            "cell-type-specific nCPM" field: the expression values for the cell types in which HPA
            classes the gene as <em>enriched</em> or <em>enhanced</em>, together with the gene's
            specificity category. A gene with <em>low specificity</em> — expressed everywhere, like
            TP53 — therefore has no cell-type rows, by design. October 2026: 86,331 rows for 17,031 genes.
          </p>
          <div className="rounded-xl bg-amber-500/5 border border-amber-500/20 p-4 text-sm">
            <p className="text-foreground/85 leading-relaxed">
              <strong>Note:</strong> HPA expression is measured at the <strong>gene level</strong>;
              isoform-specific abundance cannot be inferred directly from this profile. For
              isoform-level counts see the ccRCC per-sample table above.
            </p>
          </div>

          <h3 className="text-lg font-bold text-foreground pt-2">Interactive tissue map</h3>
          <p className="leading-relaxed">
            The top bar of every isoform page carries a human-body map that combines two layers:
          </p>
          <ul className="list-disc pl-6 space-y-2 text-foreground/90">
            <li>
              <strong>Long-read datasets (coloured):</strong> the tissues the long-read datasets
              derive from (breast tissue; kidney, for the ccRCC organoids) are highlighted in bold
              colour with animated callouts.
            </li>
            <li>
              <strong>HPA baseline (grey):</strong> the HPA baseline organs are drawn statically
              underneath as spatial context for the expression profile.
            </li>
          </ul>
        </section>

        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "AlphaFold Database", note: "dimer / multimer interface models", href: "https://alphafold.ebi.ac.uk/" },
              { label: "STRING", note: "protein–protein interaction evidence", href: "https://string-db.org/" },
              { label: "Human Protein Atlas", note: "single-cell type expression", href: "https://www.proteinatlas.org/" },
              { label: "SQANTI3", note: "per-sample full-length read counts (ccRCC)", href: "https://github.com/ConesaLab/SQANTI3" },
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
          <Link href="/documentation/ptm-accessibility" className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2">
            <ArrowLeft className="w-4 h-4" /> PTM Accessibility
          </Link>
          <div className="flex flex-col gap-1 items-end">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest">Next Section</span>
            <Link href="/documentation/drug-interactions" className="text-lg font-bold text-primary hover:underline flex items-center gap-2">
              Drug Interactions <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
