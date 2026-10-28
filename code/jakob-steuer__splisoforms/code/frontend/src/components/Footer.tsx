import React from "react";
import Link from "next/link";
import { Github } from "lucide-react";

const GROUP_URL = "https://www.fhnw.ch/plattformen/icb/datascience/";
const FHNW_URL = "https://www.fhnw.ch/en";
const FHNW_LS_URL = "https://www.fhnw.ch/en/life-sciences";
const LICENSE_URL = "https://creativecommons.org/licenses/by-nc/4.0/";
const SNF_GRANT_URL = "https://data.snf.ch/grants/grant/222878";
const REPO_URL = "https://github.com/jakob-steuer/splisoforms";
const PREPRINT_DOI = "10.64898/2026.08.06.743218";
const PREPRINT_URL = `https://doi.org/${PREPRINT_DOI}`;

const linkClass =
  "text-primary hover:text-primary/80 underline decoration-primary/30 underline-offset-4";

export function Footer() {
  return (
    <footer className="w-full border-t border-border/20 bg-background/80 backdrop-blur-xl py-10 px-6 mt-auto relative z-10 transition-colors duration-300">
      <div className="max-w-[1800px] mx-auto grid grid-cols-1 md:grid-cols-3 gap-10 text-sm text-muted-foreground">
        {/* LEFT — three uniformly-styled meta lines, vertically centered to align with the logos (spans 2/3) */}
        <div className="md:col-span-2 flex flex-col justify-center min-h-[140px] space-y-1 text-sm leading-snug text-muted-foreground/90">
          <p>
            © 2026{" "}
            <a
              href={FHNW_LS_URL}
              target="_blank"
              rel="noopener noreferrer"
              className={linkClass}
            >
              FHNW School of Life Sciences
            </a>
          </p>
          <p>
            Licensed under{" "}
            <a
              href={LICENSE_URL}
              target="_blank"
              rel="noopener noreferrer"
              className={linkClass}
            >
              CC BY-NC 4.0
            </a>
            .
          </p>
          <p>
            Supported by{" "}
            <a
              href={SNF_GRANT_URL}
              target="_blank"
              rel="noopener noreferrer"
              className={linkClass}
            >
              SNF Grant 222878
            </a>
            .
          </p>
          <p>
            Cite: Steuer J &amp; Kahraman A,{" "}
            <a
              href={PREPRINT_URL}
              target="_blank"
              rel="noopener noreferrer"
              className={linkClass}
              title="SPLISOFORMS: a Structure-Resolved Knowledge Base of Alternative Splicing Isoforms (bioRxiv preprint)"
            >
              SPLISOFORMS
            </a>
            , bioRxiv 2026, doi:{PREPRINT_DOI}.
          </p>
          <p className="pt-1 flex flex-wrap items-center gap-x-3 gap-y-1">
            <Link href="/imprint" className={linkClass}>
              Imprint &amp; Legal
            </Link>
            <span aria-hidden className="text-muted-foreground/40">·</span>
            <a
              href={REPO_URL}
              target="_blank"
              rel="noopener noreferrer"
              className={`${linkClass} inline-flex items-center gap-1.5`}
            >
              <Github className="w-3.5 h-3.5" /> Source &amp; issues on GitHub
            </a>
          </p>
        </div>

        {/* RIGHT — affiliation logos side-by-side, FHNW left (→ fhnw.ch), group right.
            Intrinsic width/height on each <img> lock the aspect ratio
            (fhnw 1280x295, group 800x392); shrink-0 prevents the flex row
            from squeezing them horizontally on narrow viewports. */}
        <div className="flex flex-row items-center justify-end gap-10 min-h-[140px]">
          <a
            href={FHNW_URL}
            target="_blank"
            rel="noopener noreferrer"
            className="hover:opacity-80 transition-opacity block shrink-0"
            title="FHNW"
          >
            <img
              src="/fhnw-logo.png"
              alt="FHNW — Fachhochschule Nordwestschweiz"
              width={1280}
              height={295}
              className="h-12 w-auto object-contain dark:invert"
            />
          </a>
          <a
            href={GROUP_URL}
            target="_blank"
            rel="noopener noreferrer"
            className="hover:opacity-80 transition-opacity block shrink-0"
            title="The AIsoform Lab"
          >
            <img
              src="/aisoform-lab-logo.png"
              alt="The AIsoform Lab"
              width={800}
              height={392}
              className="h-24 w-auto object-contain"
            />
          </a>
        </div>
      </div>
    </footer>
  );
}
