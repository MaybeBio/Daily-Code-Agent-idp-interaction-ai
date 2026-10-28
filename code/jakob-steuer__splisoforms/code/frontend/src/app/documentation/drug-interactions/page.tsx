"use client";

import { Pill, ExternalLink, ArrowLeft, Beaker } from "lucide-react";
import Link from "next/link";

export default function DrugInteractionsDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4">
        <div className="flex items-center gap-3 text-teal-500">
          <Pill className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">Drug Interactions</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          Whether a drug&apos;s binding pocket survives alternative splicing.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-8 text-foreground">
        <section className="space-y-4">
          <p className="leading-relaxed">
            This layer asks a therapeutic question directly: when a gene is alternatively spliced,
            does the pocket a drug binds still exist in the resulting isoform? The compounds are
            drawn from the <strong>PRISM Repurposing</strong> library and <strong>PDBsum</strong>{" "}
            ligands; for each drug, the target proteins to model against are taken from its known
            interactions in <strong>PubChem</strong>. Every drug–target pair is then co-folded with
            AlphaFold 3, and the drug&apos;s pocket — defined from its actual bound pose — is carried
            onto every splice isoform of the target.
          </p>
        </section>

        {/* ── Co-folded complexes ────────────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Co-folded Complexes</h2>
          <p className="leading-relaxed">
            Each complex is a single AlphaFold 3 prediction of the target protein together with
            the ligand, so the binding pose is modelled in structural context rather than docked
            after the fact. Confidence is summarised by a <strong>holo grade</strong> derived
            from the interface ipTM:
          </p>
          <div className="grid gap-3 sm:grid-cols-3">
            {[
              ["Grade A", "ipTM ≥ 0.75 — high-confidence pose.", "text-emerald-400"],
              ["Grade B", "0.50 ≤ ipTM < 0.75 — moderate.", "text-amber-400"],
              ["Grade C", "ipTM < 0.50 — low; interpret with care.", "text-rose-400"],
            ].map(([k, v, c]) => (
              <div key={k} className="p-3 rounded-lg bg-muted/20 border border-border/40">
                <span className={`font-bold text-sm ${c}`}>{k}</span>
                <p className="text-xs text-muted-foreground leading-relaxed">{v}</p>
              </div>
            ))}
          </div>
          <p className="leading-relaxed">
            Each drug page shows the compound&apos;s 2D skeleton (from its SMILES) and its real
            bound conformer lifted from the co-fold, alongside chemistry pulled from PubChem —
            molecular formula, molar mass, InChIKey, and cross-references. Two per-pocket pLDDT
            values are reported: the <strong>ligand pLDDT</strong> (confidence in the drug&apos;s
            placement) and the <strong>pocket pLDDT</strong> (confidence in the contacting
            residues).
          </p>
        </section>

        {/* ── Contact residues ───────────────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Contact Residues</h2>
          <p className="leading-relaxed">
            A pocket is the set of protein residues within <strong>4.5 Å</strong> of any ligand
            atom in the co-folded structure, stored as 1-based positions on the drugged
            transcript&apos;s sequence. These residues define what splicing must preserve for the
            drug to still bind.
          </p>
        </section>

        {/* ── Pocket fate ────────────────────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Pocket Fate Across Isoforms</h2>
          <p className="leading-relaxed">
            The pocket is anchored on the drugged transcript and aligned onto each sibling
            isoform with a splice-aware (BLOSUM62) alignment — never onto an arbitrary canonical,
            so the comparison is always drug-relative. The fraction of contact residues retained
            in aligned register gives the <strong>fate</strong>:
          </p>
          <div className="grid gap-3 sm:grid-cols-3">
            {[
              ["Preserved", "≥ 75% of contact residues retained.", "text-emerald-400"],
              ["Disrupted", "25–75% retained — partial pocket.", "text-amber-400"],
              ["Lost", "< 25% retained — pocket effectively gone.", "text-rose-400"],
            ].map(([k, v, c]) => (
              <div key={k} className="p-3 rounded-lg bg-muted/20 border border-border/40">
                <span className={`font-bold text-sm ${c}`}>{k}</span>
                <p className="text-xs text-muted-foreground leading-relaxed">{v}</p>
              </div>
            ))}
          </div>
          <div className="rounded-xl bg-amber-500/5 border border-amber-500/20 p-4 text-sm">
            <p className="text-foreground/85 leading-relaxed">
              <strong>These buckets are deliberately coarse.</strong> The 25% / 75% cut-offs are
              arbitrary triage thresholds, not a biophysical model of binding — they carry no direct
              biological guarantee. A pocket scored <em>Preserved</em> can still have lost a single
              decisive contact (one broken hydrogen bond can be enough to abolish binding), while a{" "}
              <em>Disrupted</em> pocket may retain what it needs. Read the fate as a first-pass filter,
              not a verdict. A more rigorous analysis of binding-pocket preservation is planned for a
              future release.
            </p>
          </div>
          <div className="rounded-xl bg-muted/20 border border-border/40 p-4 text-sm">
            <p className="text-foreground/85 leading-relaxed">
              <strong>No-call gate:</strong> a pocket is only scored when its residues sit in a
              confident region of both structures (per-pocket pLDDT ≥ 70). Below that, the fate is
              reported as <code className="font-mono mx-1">no_call</code> rather than guessed — a
              low-confidence region cannot be said to preserve or lose a pocket. Gene-level summaries
              are computed over decidable calls only.
            </p>
          </div>
          <p className="leading-relaxed">
            The result reproduces known biology independently — for example AR-V7, the
            treatment-resistant androgen-receptor variant, is scored <strong>Lost</strong> for its
            ligand-binding pocket because the C-terminal ligand-binding domain is spliced out.
          </p>
        </section>

        {/* ── Where it appears ───────────────────────────────────────────── */}
        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Where It Appears</h2>
          <ul className="list-disc pl-6 space-y-2 text-foreground/90">
            <li>
              The <Link href="/explore/drugs" className="text-teal-500 hover:underline">Drug Pockets lens</Link> ranks
              drugged genes by how subtractive their splicing is for drug binding.
            </li>
            <li>
              Each <strong>/drugs/[gene]</strong> page pairs the drugged-complex structure with
              the selected isoform, highlighting which contact residues survive.
            </li>
            <li>
              Gene and isoform pages carry a drug panel whenever co-fold data exists for that
              target.
            </li>
          </ul>
          <div className="rounded-xl bg-sky-500/5 border border-sky-500/20 p-4 text-sm flex gap-3 items-start">
            <Beaker className="w-4 h-4 text-sky-400 shrink-0 mt-0.5" />
            <p className="text-foreground/85 leading-relaxed">
              <strong>Scope:</strong> a pocket call speaks to structural presence, not binding
              affinity — a preserved pocket means the contact residues are there, not that the
              drug binds with unchanged potency. Co-fold selections are keyed on{" "}
              <code className="font-mono mx-1">protein_md5</code>, so predictions attach to the
              exact sequence they were computed for.
            </p>
          </div>
        </section>

        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "PRISM Repurposing", note: "drug library (Broad)", href: "https://depmap.org/repurposing/" },
              { label: "PDBsum", note: "structural ligand source", href: "https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/" },
              { label: "PubChem", note: "drug–target interactions & chemistry", href: "https://pubchem.ncbi.nlm.nih.gov/" },
              { label: "AlphaFold 3", note: "drug–protein co-folding", href: "https://github.com/google-deepmind/alphafold3" },
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
          <Link href="/documentation/interactions-expression" className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2">
            <ArrowLeft className="w-4 h-4" /> Interactions &amp; Expression
          </Link>
          <div className="flex flex-col items-end gap-1">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest text-right">Next Section</span>
            <Link href="/documentation/discovery" className="text-lg font-bold text-primary hover:underline flex items-center gap-2">
              Discovery &amp; Browsing <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
