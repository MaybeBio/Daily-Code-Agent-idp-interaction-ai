import type { Metadata } from "next";
import Link from "next/link";

export const metadata: Metadata = {
  title: "Imprint · SPLISOFORMS",
  description:
    "Legal information, operator details, and contact for the SPLISOFORMS Splicing Proteomics Knowledgebase.",
};

const SECTION = "space-y-2";
const H2 = "text-lg font-semibold text-foreground";
const P = "text-sm text-muted-foreground leading-relaxed";

export default function ImprintPage() {
  return (
    <div className="min-h-screen bg-background text-foreground">
      <div className="max-w-3xl mx-auto px-6 py-16 space-y-10">
        <header className="space-y-2">
          <h1 className="text-3xl font-bold tracking-tight">Imprint &amp; Legal Notice</h1>
          <p className="text-sm text-muted-foreground">
            Information about the operator of this site, in line with Art. 3
            UWG and Swiss good-practice for academic research portals.
          </p>
        </header>

        <section className={SECTION}>
          <h2 className={H2}>Operator</h2>
          <p className={P}>
            <strong>The AIsoform Lab</strong> &mdash; AI for Structural Oncology and Splicing<br />
            Institute for Chemistry and Bioanalytics (ICB)<br />
            School of Life Sciences (HLS)<br />
            Fachhochschule Nordwestschweiz FHNW<br />
            Hofackerstrasse 30, 4132 Muttenz, Switzerland
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>Responsible for content</h2>
          <p className={P}>
            Prof. Dr. Abdullah Kahraman &mdash; group leader<br />
            Dr. Jakob Steuer &mdash; lead developer
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>Contact</h2>
          <p className={P}>
            Email:{" "}
            <a href="mailto:abdullah.kahraman@fhnw.ch" className="text-primary hover:underline">
              abdullah.kahraman@fhnw.ch
            </a>
            <br />
            Web:{" "}
            <a
              href="https://www.fhnw.ch/plattformen/icb/datascience/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-primary hover:underline"
            >
              fhnw.ch/plattformen/icb/datascience
            </a>
            <br />
            Bug reports &amp; feature requests:{" "}
            <a
              href="https://github.com/jakob-steuer/splisoforms/issues"
              target="_blank"
              rel="noopener noreferrer"
              className="text-primary hover:underline"
            >
              open an issue on GitHub
            </a>
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>Hosting</h2>
          <p className={P}>
            Hosted on FHNW infrastructure. No third-party CDN. Server
            location: Switzerland.
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>Privacy</h2>
          <p className={P}>
            SPLISOFORMS does <strong>not</strong> use cookies, tracking pixels,
            or third-party analytics. Standard server-side access logs (nginx
            <code className="font-mono mx-1">combined</code> format: IP address,
            timestamp, request line, status code, response size, referrer, and
            user agent) are retained for up to <strong>14 days</strong> for
            operational debugging and abuse prevention, then deleted. API keys
            are transmitted in the HTTP <code className="font-mono mx-1">Authorization</code>
            header and are not written to these logs. No personal data is shared
            with third parties.
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>Funding</h2>
          <p className={P}>
            This work is supported by the Swiss National Science Foundation
            under{" "}
            <a
              href="https://data.snf.ch/grants/grant/222878"
              target="_blank"
              rel="noopener noreferrer"
              className="text-primary hover:underline"
            >
              Grant 222878
            </a>
            .
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>License &amp; intended use</h2>
          <p className={P}>
            SPLISOFORMS is freely available with no login required, and its
            analysis layers are queryable through a public REST API (endpoints
            under <code className="font-mono">/api/v1/</code>). Three licenses
            apply:
          </p>
          <ul className="text-sm text-muted-foreground leading-relaxed space-y-1.5 list-disc pl-5">
            <li>
              <strong className="text-foreground">Source code</strong> &mdash; released under the{" "}
              <a
                href="https://opensource.org/license/mit"
                target="_blank"
                rel="noopener noreferrer"
                className="text-primary hover:underline"
              >
                MIT License
              </a>{" "}
              (
              <a
                href="https://github.com/jakob-steuer/splisoforms"
                target="_blank"
                rel="noopener noreferrer"
                className="text-primary hover:underline"
              >
                github.com/jakob-steuer/splisoforms
              </a>
              ).
            </li>
            <li>
              <strong className="text-foreground">Database content &amp; annotations</strong> &mdash; provided under{" "}
              <a
                href="https://creativecommons.org/licenses/by-nc/4.0/"
                target="_blank"
                rel="noopener noreferrer"
                className="text-primary hover:underline"
              >
                CC BY-NC 4.0
              </a>
              .
            </li>
            <li>
              <strong className="text-foreground">Predicted structures</strong> &mdash; all AlphaFold 3-predicted
              structures are additionally subject to the{" "}
              <a
                href="https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md"
                target="_blank"
                rel="noopener noreferrer"
                className="text-primary hover:underline"
              >
                AlphaFold 3 Model Parameters Terms of Use
              </a>
              ; any ongoing use or adaptation of the structural models remains
              governed by those terms.
            </li>
          </ul>
          <p className={P}>
            The long-read transcriptome datasets analysed here are available from
            their original publications; reference annotation and structures
            derive from GENCODE v47, UniProt, and the AlphaFold Database (AFDB).
            Bulk structure files and full database dumps are available from the
            authors upon reasonable request. For non-clinical research use only
            &mdash; not a medical device, no diagnostic guarantee.
          </p>
        </section>

        <section className={SECTION}>
          <h2 className={H2}>Disclaimer of liability</h2>
          <p className={P}>
            Despite careful curation, the operator assumes no liability for
            the topicality, correctness, completeness, or quality of the
            information provided. Liability claims against the operator for
            material or immaterial damages arising from access to, use, or
            non-use of the published information are excluded, unless the
            operator demonstrably acted with wilful intent or gross
            negligence. External links open content that the operator has
            no influence over; responsibility lies with the linked
            providers.
          </p>
        </section>

        <div className="pt-6 border-t border-border/30">
          <Link href="/" className="text-sm text-primary hover:underline">
            &larr; Back to home
          </Link>
        </div>
      </div>
    </div>
  );
}
