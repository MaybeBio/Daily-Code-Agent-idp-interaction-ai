# Changelog

All notable changes to this project are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed

- Replace certainty language and residue-independence assumptions with a reproducible protein-
  cluster sensitivity; distinguish DisProt reference annotations from complete biological ground
  truth and label evidence-code mechanisms as exploratory.

## [1.0.0] - 2026-08-13

### Added

- Stable product release around the unchanged frozen v0.1 structural-biology study.
- Release audit plus contribution, conduct and security documentation.

### Fixed

- Hardened the transitive `image-size` parsers against malformed ICNS and ISO-BMFF input,
  added executable security regression probes, and kept the production audit strict for every
  other high-severity advisory.
- Made frozen-registry verification portable across operating systems by comparing floating-point
  values with a narrow numerical tolerance, including an explicitly bounded tolerance for seeded
  classifier-bootstrap interval boundaries, while retaining exact checks for all other data.
- Aligned the Cloudflare compatibility date with the version supported by the pinned runtime.

### Added

- `scripts/generate_threshold_sweep.py` and `site/app/threshold-explorer.tsx`: an interactive
  pLDDT-threshold explorer on the classifier section of the site. The script reads the real
  `data/external/joined_residues.csv` and precomputes precision/recall/F1/MCC and confusion
  counts at every integer threshold from 0 to 100, overall and split by the same protein-level
  conditional-folding-flag breakdown `docs/research-report.md` already reports (verified
  internally consistent with the frozen registry's own point estimates at threshold=70: see
  `site/tests/threshold-explorer.test.mjs`), into `site/app/data/threshold-sweep.json`. The
  slider is a real `<input type="range">` with a visible label, `aria-valuetext` describing the
  selected threshold and its headline stats in words, an `aria-live="polite"` region for
  screen-reader users, a "Reset to preregistered threshold (70)" control, a vertical marker on
  the existing pLDDT ECDF chart showing where the threshold falls relative to both
  distributions, and a static, always-visible reference table (thresholds 50/70/90, aggregate
  and by conditional-folding flag) as a non-interactive fallback. This is a presentation-layer,
  site-scoped addition: it does not touch `src/foldings_edge/`, the frozen registry, or the
  research docs, and does not change any reported headline number.
- `site/app/hero-chain.tsx`: a real-data hero animation — a 60-residue slice (residues
  645–704) of UniProt Q92793's actual per-residue AlphaFold2 pLDDT (`data/external/
  joined_residues.csv`), rendered as a connected chain on canvas whose per-residue jitter
  amplitude is driven by `(100 − pLDDT)` for that residue, so confidently-predicted residues
  stay nearly still and low-pLDDT residues (the curated-disorder region this slice runs into,
  starting at residue 674) wiggle visibly. Respects `prefers-reduced-motion` (a single static
  frame instead of the loop), sized to its own numerically-verified content extent, and paired
  with an accessible caption stating the real accession, residue range, and data source
  (canvas itself is `aria-hidden`). This is a visualization-only, site-scoped change: no new
  statistics, citations, or changes to `src/foldings_edge/`, the frozen registry, or the
  research docs.

## [0.1.0] - 2026-08-13

### Added

- Preregistered research protocol (`docs/research-protocol.md`): data sources, exact sample
  construction, exclusion criteria, ground-truth region definition, and exact statistical tests,
  fixed before any analysis was run, plus an explicit framing of this project's relationship to
  Alderson et al. 2023 (PNAS) as a qualitative, not numerical, comparison.
- `scripts/fetch_data.py`: fetches human DisProt entries (`release=current`) and joins them to
  AlphaFold DB per-residue pLDDT by UniProt accession, with disclosed exclusion logging.
  Discovered and worked around two undocumented API behaviors during the real fetch: DisProt's
  search endpoint ignoring `limit`/`offset` under an `organism` filter, and AlphaFold DB's
  metadata endpoint returning UniProt isoform records rather than the expected long-protein
  sequence-fragment records.
- `data/provenance.json`: full data provenance for both sources, license status, sample
  construction, and the structural verification performed.
- `data/external/joined_residues.csv` and `data/external/exclusions_log.json`: the committed,
  derived 387-protein / 228,662-residue joined dataset (CC BY 4.0 downstream of two CC BY 4.0
  sources) and its exclusion log.
- `src/foldings_edge`: joined-residue loading (`dataset.py`), the Mann–Whitney/KS/bootstrap/
  classifier statistical pipeline (`stats.py`), the result-registry builder (`registry.py`), and
  a small CLI (`cli.py`).
- `reports/v0.1-foldings-edge-registry.json`: the frozen, reproducible v0.1.0 result registry,
  generated from the real joined dataset.
- `docs/research-report.md`: hypothesis-by-hypothesis findings, including where the pLDDT/
  disorder relationship breaks down (by conditional-folding proxy, by DisProt evidence code, and
  by individual protein), and full limitations.
- `site/`: an accessible Next.js (vinext) interactive site covering the pLDDT-vs-disorder
  distribution comparison, classifier performance with uncertainty, full data tables, and a
  provenance page, built for Cloudflare Workers deployment as `foldings-edge-interactive`.
- Full repository hygiene: `pyproject.toml` (ruff, mypy strict, pytest with a 95% coverage gate),
  `CITATION.cff`, `ACCESSIBILITY.md`, and CI workflows (`ci.yml`, `codeql.yml`).
