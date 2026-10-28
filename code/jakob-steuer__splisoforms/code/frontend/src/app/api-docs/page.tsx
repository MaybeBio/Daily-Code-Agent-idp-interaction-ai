"use client";

import { Code2, ArrowLeft, KeyRound, Gauge, ShieldCheck, Terminal } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import Link from "next/link";

const endpoints = [
  { method: "GET", path: "/genes/{symbol}", desc: "Gene summary and its isoforms." },
  { method: "GET", path: "/genes/{symbol}/isoforms", desc: "Paginated isoforms of a gene (limit ≤ 100)." },
  { method: "GET", path: "/isoforms/{transcript_id}", desc: "Core record + summary metrics (Δlength, NMD, domains, neoantigens, PTM, pLDDT)." },
  { method: "GET", path: "/isoforms/{transcript_id}/structure", desc: "Redirect (307) to the predicted mmCIF structure." },
  { method: "GET", path: "/isoforms/{transcript_id}/domains", desc: "Pfam + TED domain annotations." },
  { method: "GET", path: "/isoforms/{transcript_id}/neoantigens", desc: "Top predicted neoantigens (limit ≤ 50)." },
  { method: "GET", path: "/search?query=", desc: "Search genes + isoforms (limit ≤ 25)." },
];

export default function ApiDocsPage() {
  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-3xl mx-auto px-6 py-12 space-y-10">
        <div>
          <Link href="/" className="inline-flex items-center gap-1.5 text-xs text-muted-foreground hover:text-primary mb-6">
            <ArrowLeft className="w-3.5 h-3.5" /> Back to Dashboard
          </Link>
          <div className="flex items-center gap-3">
            <Code2 className="w-7 h-7 text-secondary" />
            <h1 className="text-4xl font-extrabold tracking-tight">REST API</h1>
            <Badge variant="outline" className="ml-1">v1</Badge>
          </div>
          <p className="text-muted-foreground mt-3 leading-relaxed">
            A curated, open, read-only HTTP API for programmatic access to genes, isoforms, summary metrics,
            and predicted structures. Base URL: <code className="font-mono text-foreground">https://splisoforms.org/api/v1</code>.
            Isoform endpoints accept either a transcript id (ENST…/PB…) or a transcript name (e.g. F2-201).
            All responses are JSON, except the structure endpoint, which 307-redirects to an mmCIF file.
          </p>
        </div>

        {/* Auth */}
        <section className="space-y-3">
          <h2 className="flex items-center gap-2 text-lg font-bold"><KeyRound className="w-4 h-4 text-primary" /> Authentication</h2>
          <p className="text-sm text-muted-foreground leading-relaxed">
            <strong className="text-foreground">No key is required</strong> — the API is open and read-only.
            Anonymous requests are rate-limited per IP (see below). For higher sustained or bulk volume, request
            an API key by emailing{" "}
            <a href="mailto:abdullah.kahraman@fhnw.ch?subject=SPLISOFORMS%20API%20key%20request" className="text-primary hover:underline">abdullah.kahraman@fhnw.ch</a>{" "}
            (FHNW · Kahraman Lab) with your name, affiliation, and intended use — full contact details are on the{" "}
            <Link href="/imprint" className="text-primary hover:underline">Imprint</Link> page. Send the key as a
            Bearer token:
          </p>
          <pre className="rounded-xl border border-border/50 bg-card/50 p-4 text-xs font-mono overflow-x-auto">
{`Authorization: Bearer spl_live_xxxxxxxxxxxxxxxxxxxx   # optional — higher-volume tier only`}
          </pre>
        </section>

        {/* Rate limits */}
        <section className="space-y-3">
          <h2 className="flex items-center gap-2 text-lg font-bold"><Gauge className="w-4 h-4 text-amber-500" /> Rate limits</h2>
          <p className="text-sm text-muted-foreground leading-relaxed">There are two access tiers:</p>
          <div className="rounded-xl border border-border/40 divide-y divide-border/30 overflow-hidden">
            <div className="flex flex-col sm:flex-row sm:gap-4 px-4 py-3">
              <span className="font-semibold text-foreground sm:w-28 shrink-0">Anonymous</span>
              <span className="text-sm text-muted-foreground leading-relaxed">
                Open, no key. Rate-limited per IP; requests over the limit return{" "}
                <code className="font-mono">429 Too Many Requests</code>. Fine for interactive and light
                programmatic use.
              </span>
            </div>
            <div className="flex flex-col sm:flex-row sm:gap-4 px-4 py-3">
              <span className="font-semibold text-foreground sm:w-28 shrink-0">API key</span>
              <span className="text-sm text-muted-foreground leading-relaxed">
                A higher per-key quota for sustained or bulk access — request one as described above.
              </span>
            </div>
          </div>
        </section>

        {/* Endpoints */}
        <section className="space-y-3">
          <h2 className="flex items-center gap-2 text-lg font-bold"><Terminal className="w-4 h-4 text-emerald-500" /> Endpoints</h2>
          <div className="rounded-xl border border-border/40 divide-y divide-border/30 overflow-hidden">
            {endpoints.map((e) => (
              <div key={e.path} className="flex flex-col sm:flex-row sm:items-center gap-1 sm:gap-3 px-4 py-3">
                <div className="flex items-center gap-2 shrink-0">
                  <Badge className="bg-emerald-500/15 text-emerald-400 border-emerald-500/30 text-[10px]">{e.method}</Badge>
                  <code className="font-mono text-sm text-foreground">{e.path}</code>
                </div>
                <span className="text-xs text-muted-foreground sm:ml-auto sm:text-right">{e.desc}</span>
              </div>
            ))}
          </div>
        </section>

        {/* Example */}
        <section className="space-y-3">
          <h2 className="text-lg font-bold">Example</h2>
          <pre className="rounded-xl border border-border/50 bg-card/50 p-4 text-xs font-mono overflow-x-auto leading-relaxed">
{`# no key needed for standard use
curl https://splisoforms.org/api/v1/genes/TP53
curl https://splisoforms.org/api/v1/isoforms/TP53-208

# fetch the structure (follow the 307 redirect to the mmCIF)
curl -L https://splisoforms.org/api/v1/isoforms/TP53-208/structure -o model.cif

# only for higher-volume access: add your key
curl -H "Authorization: Bearer $SPLISOFORMS_KEY" \\
     "https://splisoforms.org/api/v1/search?query=kinase"`}
          </pre>
        </section>

        {/* Scope note */}
        <section className="space-y-3">
          <h2 className="flex items-center gap-2 text-lg font-bold"><ShieldCheck className="w-4 h-4 text-primary" /> Scope</h2>
          <p className="text-sm text-muted-foreground leading-relaxed">
            The API is intentionally coarse-grained: it serves gene/isoform records, summary metrics, and
            structures. It does <strong>not</strong> provide bulk database dumps, raw per-residue mass data,
            or write/compute access. For full interactive exploration, use the web interface.
          </p>
        </section>
      </div>
    </div>
  );
}
