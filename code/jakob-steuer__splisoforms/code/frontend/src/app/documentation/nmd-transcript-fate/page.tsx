"use client";

import {
  ExternalLink,
  HelpCircle,
  Microscope,
  OctagonAlert,
  ShieldCheck,
  Sparkles,
  Zap,
} from "lucide-react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";

// The three NMD calls that actually occur in the database. Two come straight from
// the 50-nt rule; "Unknown" covers every case where the rule can't be applied.
// (Kept in sync with backend/app/analysis/nmd_prediction.py.)
const statuses = [
  {
    name: "NMD Target",
    icon: OctagonAlert,
    box: "border-red-500/30 bg-red-500/5",
    accent: "text-red-400",
    def: "A premature stop codon lands more than 50 nt upstream of the last exon–exon junction, outside any escape zone.",
    outcome: "Transcript likely degraded — protein not produced.",
  },
  {
    name: "NMD Escape",
    icon: ShieldCheck,
    box: "border-emerald-500/30 bg-emerald-500/5",
    accent: "text-emerald-400",
    def: "The stop sits in the last exon, within 50 nt of the last junction, or in the 5′-proximal zone — or the transcript is single-exon, with no junction to trigger decay. GENCODE's curated NMD_exception tag also forces this call.",
    outcome: "Transcript likely translated.",
  },
  {
    name: "Unknown",
    icon: HelpCircle,
    box: "border-amber-500/30 bg-amber-500/5",
    accent: "text-amber-400",
    def: "The 50-nt rule couldn't be applied — the exon or CDS coordinates are missing, the annotated CDS has no stop codon (GENCODE cds_end_NF / mRNA_end_NF, or a CDS running past the transcript end), or the ORF is a short-ORF candidate.",
    outcome: "No call is made, rather than a guess.",
  },
];

// Verdict counts in the October 2026 release (splisoforms_v2), out of 131,554 isoforms.
const verdictCounts = [
  { label: "NMD target", value: "9,262", cls: "text-red-400" },
  { label: "NMD escape", value: "78,146", cls: "text-emerald-400" },
  { label: "Unknown", value: "44,146", cls: "text-amber-400" },
];

const truncations = [
  { name: "Intact", badge: "bg-emerald-500/10 text-emerald-400 border-emerald-500/30", desc: "≥ 90% of canonical length." },
  { name: "Minor",  badge: "bg-amber-500/10 text-amber-400 border-amber-500/30",       desc: "70–90% of canonical length (e.g. an alternative transcription start)." },
  { name: "Severe", badge: "bg-rose-500/10 text-rose-400 border-rose-500/30",          desc: "< 70% of canonical length — suggests a premature stop or a major exon skip." },
];

export default function TranslationFidelityDoc() {
  return (
    <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-700">
      {/* Header */}
      <div className="space-y-4 px-3">
        <div className="flex items-center gap-3 text-amber-500">
          <OctagonAlert className="w-8 h-8" />
          <h1 className="text-4xl font-extrabold tracking-tight">NMD &amp; Transcript Fate</h1>
        </div>
        <p className="text-xl text-foreground/80 leading-relaxed italic">
          Predicting Nonsense-Mediated Decay susceptibility and truncation severity.
        </p>
      </div>

      <div className="prose prose-invert max-w-none space-y-12 text-foreground">

        <section className="space-y-4">
          <p className="leading-relaxed">
            Nonsense-mediated decay (NMD) is a surveillance pathway that degrades mRNAs whose stop
            codon is recognised as premature — typically one sitting well upstream of the last
            exon–exon junction. SPLISOFORMS predicts each isoform's NMD susceptibility from the
            established 50-nucleotide rule, and separately tracks how truncated a surviving isoform
            is relative to its canonical counterpart.
          </p>
        </section>

        <section className="space-y-8">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Zap className="w-5 h-5 text-amber-500" />
            The 50-Nucleotide Rule
          </h2>
          <p>
            An isoform is flagged as an NMD target when it carries a <strong>premature termination
            codon (PTC)</strong> — a stop lying <strong>more than 50 nt upstream</strong> of the
            final exon–exon junction (<a href="https://doi.org/10.1016/S0968-0004(98)01208-0" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline">Nagy &amp; Maquat, 1998</a>).
            The distance is measured from the first base of the stop codon to the last junction.
            Three exceptions let a transcript escape decay even with such an upstream stop:
          </p>
          <ul className="list-disc list-inside space-y-2 text-foreground/85">
            <li><strong>5′-proximal escape</strong>: the stop codon lies within 150 nt of the start codon, close enough for the ribosome to reinitiate.</li>
            <li><strong>Single-exon transcripts</strong>: with no exon–exon junctions, there is nothing to mark the stop as premature.</li>
            <li><strong>GENCODE NMD_exception</strong>: reference transcripts that GENCODE has curated as escaping NMD despite a PTC-like stop (typically a 3′-UTR intron on a MANE transcript) keep that curated call over the arithmetic.</li>
          </ul>

          <h3 className="text-lg font-bold text-foreground pt-2">Where the coordinates come from</h3>
          <p className="leading-relaxed">
            The rule needs real coding-sequence (CDS) coordinates on the transcript's exon
            model — it is never inferred from protein length alone:
          </p>
          <ul className="list-disc list-inside space-y-2 text-foreground/85">
            <li><strong>GENCODE v47</strong>: the annotated CDS, projected from the GTF onto the transcript. Transcripts tagged <code className="font-mono text-xs">cds_end_NF</code> or <code className="font-mono text-xs">mRNA_end_NF</code> have no annotated stop codon and are <em>Unknown</em>.</li>
            <li><strong>ccRCC organoids</strong>: the author-curated ORFs that were actually folded, looked up by protein sequence hash and PacBio id in the sample they came from.</li>
            <li><strong>Breast cancer</strong>: no TransDecoder ORF output is loaded yet, so every breast-cancer-only isoform is <em>Unknown</em>.</li>
          </ul>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Status Values</h2>
          <p className="leading-relaxed">
            Every isoform receives one of three calls. Two come directly from the 50-nt rule; the
            third, <strong>Unknown</strong>, is used whenever the rule cannot be applied.
          </p>
          <div className="grid gap-4 md:grid-cols-3">
            {statuses.map((s) => (
              <div key={s.name} className={`rounded-xl border ${s.box} p-4 space-y-2`}>
                <div className="flex items-center gap-2">
                  <s.icon className={`w-5 h-5 ${s.accent}`} />
                  <h3 className={`text-sm font-bold uppercase tracking-widest ${s.accent}`}>{s.name}</h3>
                </div>
                <p className="text-xs text-muted-foreground leading-relaxed">{s.def}</p>
                <p className="text-xs italic text-foreground/70 leading-relaxed">{s.outcome}</p>
              </div>
            ))}
          </div>

          <div className="rounded-xl border border-border/40 bg-muted/10 p-4 not-prose">
            <p className="text-[10px] font-bold uppercase tracking-widest text-muted-foreground mb-2">
              Verdicts · October 2026 release · 131,554 isoforms
            </p>
            <div className="grid grid-cols-3 gap-3">
              {verdictCounts.map((v) => (
                <div key={v.label}>
                  <div className={`text-xl font-bold font-mono ${v.cls}`}>{v.value}</div>
                  <div className="text-xs text-muted-foreground">{v.label}</div>
                </div>
              ))}
            </div>
            <p className="text-xs text-muted-foreground leading-relaxed mt-3">
              Verdicts exist for GENCODE and ccRCC isoforms. The <em>Unknown</em> count includes all
              23,148 breast-cancer-only isoforms and about 18,500 GENCODE transcripts whose annotated
              CDS has no stop codon; before October 2026 the latter were reported as <em>Escape</em>.
            </p>
          </div>

          <div className="rounded-xl border border-amber-500/20 bg-amber-500/5 p-4 text-sm text-muted-foreground leading-relaxed space-y-2">
            <p className="font-bold text-foreground">When information is missing, the call is Unknown</p>
            <p>
              The 50-nt rule needs both the transcript's exon structure and its CDS coordinates,
              including an annotated stop codon. When any of those is missing, the isoform is left{" "}
              <strong>Unknown</strong> rather than guessed — a CDS that runs past the end of the
              transcript is never reported as <em>Escape</em>.
            </p>
            <p>
              This currently affects the whole <strong>breast-cancer (Head et al. 2026)</strong>{" "}
              dataset: its TransDecoder ORF/CDS predictions are not yet loaded, so those isoforms
              remain Unknown until that data lands. Short translated ORFs (smORF / uORF candidates)
              are also routed to Unknown, because the 50-nt rule was never validated for them — see below.
            </p>
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground">Truncation Severity</h2>
          <p className="leading-relaxed">
            Protein-length divergence vs. the canonical sequence is reported as 
            <code className="font-mono ml-1">truncation_severity</code>. A short isoform may have a legitimate 
            alternative TSS without triggering NMD, whereas a full-length isoform could still harbor a PTC.
          </p>
          <div className="grid gap-3">
            {truncations.map((t) => (
              <div key={t.name} className="flex gap-4 p-4 rounded-xl border border-border/40 bg-muted/10 items-start">
                <Badge className={`${t.badge} text-[10px] uppercase font-bold tracking-widest px-3 py-1 ring-1 ring-border/50 whitespace-nowrap`}>
                  {t.name}
                </Badge>
                <p className="text-sm leading-relaxed text-muted-foreground font-medium">{t.desc}</p>
              </div>
            ))}
          </div>
        </section>

        <section className="space-y-6">
          <h2 className="text-2xl font-bold text-foreground flex items-center gap-2">
            <Sparkles className="w-5 h-5 text-violet-400" />
            Microproteins &amp; uORF Candidates
          </h2>
          <p className="leading-relaxed">
            The 50-nt rule was developed for canonical protein-coding transcripts, not for short
            translated ORFs — applying it to them would produce misleading calls. Such short ORFs
            are increasingly recognised as a genuine layer of the proteome — the expanding catalogue
            of{" "}
            <a href="https://www.nature.com/articles/s41586-026-10459-x" target="_blank" rel="noopener noreferrer" className="text-primary hover:underline font-semibold">microproteins and peptideins</a>.
            SPLISOFORMS flags these candidates separately and routes them to <strong>Unknown</strong>.
          </p>
          <div className="grid gap-3 md:grid-cols-2">
            <div className="rounded-xl border border-violet-500/30 bg-violet-500/5 p-4">
              <div className="flex items-center gap-2 mb-2">
                 <Sparkles className="w-4 h-4 text-violet-400" />
                 <span className="text-sm font-bold text-violet-400">is_candidate_smorf</span>
              </div>
              <p className="text-sm text-muted-foreground leading-relaxed">
                Isoform &lt; 100 aa and canonical &lt; 150 aa. May encode <i>bona fide</i> microproteins or peptideins.
              </p>
            </div>
            <div className="rounded-xl border border-fuchsia-500/30 bg-fuchsia-500/5 p-4">
              <div className="flex items-center gap-2 mb-2">
                 <Microscope className="w-4 h-4 text-fuchsia-400" />
                 <span className="text-sm font-bold text-fuchsia-400">is_candidate_uorf</span>
              </div>
              <p className="text-sm text-muted-foreground leading-relaxed">
                Isoform &lt; 50 aa with a long canonical (≥ 150 aa). Heuristic flag for NMD-resistant upstream ORFs.
              </p>
            </div>
          </div>
        </section>

        {/* ── Tools & sources ────────────────────────────────────────────── */}
        <section className="space-y-3 pt-2">
          <h2 className="text-lg font-bold text-foreground">Tools &amp; sources</h2>
          <ul className="grid gap-2 sm:grid-cols-2 not-prose list-none pl-0">
            {[
              { label: "GENCODE v47", note: "reference CDS coordinates & transcript tags", href: "https://www.gencodegenes.org/human/release_47.html" },
              { label: "TransDecoder", note: "ORF / CDS prediction from transcripts", href: "https://github.com/TransDecoder/TransDecoder" },
              { label: "SQANTI3", note: "long-read isoform & ORF-type classification", href: "https://github.com/ConesaLab/SQANTI3" },
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

        <div className="pt-10 border-t border-border/20 flex items-center justify-between gap-4">
          <Link href="/documentation/functional-domains" className="text-sm font-bold text-muted-foreground hover:text-primary transition-colors flex items-center gap-2">
            <span className="opacity-60">&larr;</span> Prev: Functional Domains
          </Link>
          <div className="flex flex-col items-end gap-1">
            <span className="text-[10px] uppercase font-bold text-muted-foreground tracking-widest text-right">Next Section</span>
            <Link href="/documentation/neoantigens" className="text-lg font-bold text-primary hover:underline flex items-center gap-2">
              Neoantigens <ExternalLink className="w-4 h-4" />
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}
