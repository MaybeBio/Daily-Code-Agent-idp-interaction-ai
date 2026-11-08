# Foldings Edge

<p><a href="https://github.com/lindgreendavid/lindgreendavid/tree/main/brand"><img src="https://raw.githubusercontent.com/lindgreendavid/lindgreendavid/main/brand/lab-notes-mark.svg" width="52" align="right" alt="Lab Notes research-cycle mark"></a></p>

**Part of the [Lab Notes Research Portfolio](https://notes.lindgreendavid.com/)** · Structural biology · Question → evidence → finding → boundary

A reproducible reanalysis testing whether AlphaFold2's per-residue pLDDT confidence score
predicts experimentally-curated intrinsic disorder (DisProt) — and where that relationship
breaks down.

**[Open the live interactive site](https://foldings-edge-interactive.lindgreendavid.workers.dev)** · **[Read the plain-language write-up](https://notes.lindgreendavid.com/posts/foldings-edge-plddt-disorder)**

**Stable product release:** [v1.0.0](https://github.com/lindgreendavid/foldings-edge/releases/tag/v1.0.0) · **Study:** unchanged frozen v0.1 registry.

**Research question:** does AlphaFold2's pLDDT confidence score predict DisProt-curated
intrinsic disorder for human proteins, and where does that relationship break down? This
project is motivated by, but is **not a replication of**, T. Reid Alderson et al.'s 2023 PNAS
paper "Systematic identification of conditionally folded intrinsically disordered regions by
AlphaFold2" (DOI: [10.1073/pnas.2304302120](https://doi.org/10.1073/pnas.2304302120)). It runs
an independent, preregistered test of the same underlying idea — pLDDT vs. curated disorder —
on a fresh 2026 DisProt snapshot joined to real AlphaFold DB predictions, using its own sample,
ground-truth scope, and classifier, and compares its results to that paper only qualitatively.

**Headline finding (v0.1.0):** on 387 human proteins (228,662 residues, 400-of-1,339 fixed-seed
DisProt sample), pLDDT is **substantially lower inside curated disorder regions** (median 45.1
vs. 88.8). The direction survives a protein-level paired sensitivity (median within-protein
difference 36.47, cluster-bootstrap 95% CI 34.29–40.07), so it is not an artifact of treating
all 228,662 correlated residues as independent observations.
A fixed pLDDT<70 classifier catches most curated disorder (76.5% recall) but is a **noisy, only
moderately strong** predictor of it (34.1% precision, F1 0.472, MCC 0.367) — H2 confirmed above
chance, not strongly. Recall is lower on disorder regions that co-occur with a conditional-folding annotation
elsewhere on the same protein (0.681 vs 0.780; protein-clustered 95% interval of the difference 0.0005 to 0.212, AUC
difference interval -0.005 to 0.115) and on disorder supported by hydrogen-deuterium exchange mass spectrometry
(35.1% recall, but from only 4 proteins). (An earlier version of this README said precision and MCC "drop sharply" for
the conditional-folding group; that is a size artifact of scoring each subgroup against a shared non-disorder pool, and
is corrected in the report's erratum.) Those subgroup patterns are exploratory and do not establish why performance differs; the overall
direction is only a qualitative match to Alderson et al., not a numerical replication. Full
reasoning, every hypothesis's disposition, and every limitation: see
[`docs/research-report.md`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/docs/research-report.md).

**Paper:** [When Precision Misleads (PDF)](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/paper/paper.pdf) · [citation and status](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/paper/README.md) · [version history](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/history.md)

## What this contributes

- An independently-run test of the same general relationship a significant published paper
  investigates, using two real, live public data sources joined at reproduction time, rather
  than a restatement of the paper's own claims: does a simple pLDDT<70 classifier actually
  predict curated disorder on a fresh sample, and specifically where does it fail? (Answer: the
  distributional difference is unambiguous; the classifier works better than chance but is far
  from a clean disorder detector, and recall is lower where the motivating
  paper's own theory predicts it should be, though with 61 conditional proteins the evidence for that gap is weak.)
- Two disclosed, real data-access deviations, found and fixed during the actual fetch rather than
  hidden: DisProt's search API ignores `limit`/`offset` under an `organism` filter, and AlphaFold
  DB's metadata endpoint returns UniProt isoform records rather than the long-protein sequence
  fragments this project initially expected. Both are documented in
  [`data/provenance.json`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/data/provenance.json) and fixed in
  [`scripts/fetch_data.py`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/scripts/fetch_data.py) before any result was generated.
- A concrete, per-protein and per-evidence-code failure analysis, not just an aggregate
  precision/recall number — see `docs/research-report.md` for named examples of both failure
  directions (curated disorder that received confident pLDDT, and non-disorder that received low
  pLDDT).
- What it does **not** contribute: a replication of Alderson et al. (2023)'s exact classifier,
  database combination, or ROC/precision numbers; any claim about AlphaFold3; or any claim about
  intrinsically disordered proteins beyond this specific 387-protein sample.

## What's here

| Path | What it is |
| --- | --- |
| [`docs/research-protocol.md`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/docs/research-protocol.md) | The preregistered hypotheses, data sources, sample construction, exclusion criteria, and exact statistical tests — written and committed before any result existed. |
| [`docs/research-report.md`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/docs/research-report.md) | What the frozen registry actually shows, hypothesis by hypothesis, including the breakdown analyses and every limitation. |
| [`data/provenance.json`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/data/provenance.json) | Full data provenance: sources, access date, license status, sample construction, and the two disclosed API-behavior deviations found during the real fetch. |
| [`scripts/fetch_data.py`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/scripts/fetch_data.py) | Fetches and joins DisProt human disorder annotations with AlphaFold DB per-residue pLDDT at reproduction time. |
| [`data/external/joined_residues.csv`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/data/external/joined_residues.csv) | The committed, derived 387-protein / 228,662-residue joined dataset (CC BY 4.0 downstream of two CC BY 4.0 sources). |
| [`src/foldings_edge/`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/src/foldings_edge) | The Python package: residue loading (`dataset.py`), the Mann–Whitney/KS/bootstrap/classifier statistics (`stats.py`), and the registry builder (`registry.py`). |
| [`tests/`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/tests) | Unit tests (small synthetic fixtures, not real network calls) and a byte-comparison test against the frozen registry. |
| [`reports/v0.1-foldings-edge-registry.json`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/reports/v0.1-foldings-edge-registry.json) | The frozen, deterministic analysis output, generated from the real joined dataset. |
| [`site/`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/site) | An accessible Next.js (vinext) interactive site: a hero animation of a real protein's per-residue pLDDT driving a wiggling chain, pLDDT-vs-disorder distribution view, classifier performance with uncertainty, an interactive pLDDT-threshold explorer (drag a slider to see live precision/recall/F1/MCC and confusion counts at any threshold, split by the conditional-folding-flag breakdown, with a one-click reset to the preregistered threshold of 70), full data tables, and a provenance page, built for Cloudflare Workers. |

## Reproduce the analysis locally

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'

# Fetch and join the real data (400-of-1,339 fixed-seed DisProt sample; ~5-10 minutes,
# polite to both live APIs). The committed data/external/joined_residues.csv already
# contains the result of this exact run:
python scripts/fetch_data.py

# Regenerate the registry and confirm it matches the committed one:
python scripts/generate_registry.py --output /tmp/v0.1-foldings-edge-registry.json \
    --residues data/external/joined_residues.csv
cmp reports/v0.1-foldings-edge-registry.json /tmp/v0.1-foldings-edge-registry.json  # should be silent

# Or just summarize the dataset:
foldings-edge summarize data/external/joined_residues.csv
```

This is deterministic given a fixed joined dataset (fixed seed throughout); CI regenerates and
byte-compares the registry from the committed dataset on every push.

## Run the interactive site locally

```bash
cd site
pnpm install
pnpm run dev      # local development server
pnpm run build    # production build (Cloudflare Workers target)
pnpm run test     # build + node --test
pnpm run lint     # eslint
```

`site/wrangler.jsonc` is configured for deployment to Cloudflare Workers as
`foldings-edge-interactive`. This repository does not run `wrangler deploy` — that is a manual
step the maintainer runs after reviewing a build.

## Quality gates

```bash
ruff check .
ruff format --check .
mypy src
pytest                      # includes a 95% coverage gate
python -m build

cd site && pnpm run lint && pnpm run build && pnpm run test
```

## Data sources and citations

- **DisProt** (ground truth): human entries, `release=current`, accessed 2026-08-13.
  CC BY 4.0. https://disprot.org
- **AlphaFold Protein Structure Database** (per-residue pLDDT): accessed 2026-08-13. CC BY 4.0.
  Jumper, J., Evans, R., Pritzel, A. et al. "Highly accurate protein structure prediction with
  AlphaFold." *Nature* 596, 583–589 (2021). DOI:
  [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2).
- **Motivating citation** (not replicated, compared qualitatively): Alderson, T.R., Pritišanac,
  I., Kolarić, Đ., Moses, A.M., Forman-Kay, J.D. "Systematic identification of conditionally
  folded intrinsically disordered regions by AlphaFold2." *PNAS* 120(44):e2304302120 (2023).
  DOI: [10.1073/pnas.2304302120](https://doi.org/10.1073/pnas.2304302120).

See [`data/provenance.json`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/data/provenance.json) for the full provenance record and
[`CITATION.cff`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/CITATION.cff) for citing this software.

## Scope and limitations (short version — full version in the report)

Single DisProt snapshot, human only, 400-of-1,339 fixed-seed sample, 387 proteins after
disclosed exclusions. No claim about AlphaFold3. No claim about intrinsically disordered
proteins beyond this dataset. The "conditional folding" breakdown uses a protein-level proxy,
not Alderson et al.'s own region-level five-database definition. See
[`docs/research-protocol.md`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/docs/research-protocol.md#scope-boundaries-declared-before-results)
and [`docs/research-report.md`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/docs/research-report.md#limitations) for the complete, disclosed
list.

## License

MIT. See [`LICENSE`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/LICENSE). DisProt and AlphaFold DB data are not redistributed under this
license — see [`data/provenance.json`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/data/provenance.json) for their own CC BY 4.0 terms.

## Citation

See [`CITATION.cff`](https://github.com/lindgreendavid/foldings-edge/blob/HEAD/CITATION.cff).

## Deploying the interactive site

The interactive site in `site/` is a Cloudflare Worker. Its home page is large, and on the free Workers plan rendering it on every visit can exceed the CPU limit: in November 2026 the live site returned empty 503 errors for that reason. The fix is to save the home page as a static file that Cloudflare serves without running the Worker.

```bash
cd site
pnpm install
pnpm run deploy:static                  # build, save the home page as a static file, then deploy with wrangler
node scripts/deploy-static.mjs --no-deploy   # same, but stop before deploying, to inspect dist/client/index.html
```

Do not deploy with `pnpm run deploy:vinext` or `vinext-cloudflare deploy`: they build again first and delete the saved page, and the 503 errors can come back. If the site moves to the paid Workers plan, the normal deploy works and this script is no longer needed.
