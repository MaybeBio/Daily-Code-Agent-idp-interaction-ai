# IDR Explorer

The disordered regions of the reviewed human proteome, read two ways, with the modifications that tune their charge:

- **IDRs**: every MobiDB disordered region of 20 or more residues (33,542 on 12,297 proteins), with its charge, its
  charge patterning, its acidic stretches, RG/RGG clusters and basic patches, and its sites.
- **Domain flanks**: the residues beside every folded-domain edge (121,444; the window a domain-tethered IDR, or
  dtIDR, would occupy), scored from sequence for basic and for acidic character. The top tenth of each score are the
  **basic flanks** (arginine and lysine writers) and **acidic flanks** (CK2, glutamylases).
- **Charge, in words**: every IDR and flank is acidic, basic, mixed or has few charges (net charge per residue, then
  the share charged; `idrkit.charge_type`), and in the Das–Pappu diagram of states with plain names and the shape each
  class predicts (localCIDER; Holehouse et al. 2017), with κ and Ω.

The site is in [`docs/`](https://github.com/Shechterlab/idr-explorer/blob/HEAD/docs) (GitHub Pages: `main`, folder `/docs`). Its Overview first asks whether a domain flank is a unit at
all (read without a window, family by family for every InterPro entry, and by charge class), then where each
modification that changes charge sits (acidic and basic flanks, blocky IDRs, acidic stretches, basic patches, RG/RGG
clusters, kinase motifs, charge blocks, occupancy, protein classes); every test runs again in the browser on a
selection of proteins and evidence. Tables of every protein, flank, IDR and site, each under a plot that selects (hover
a point, click to open it, drag to pick rows), and a page for every protein: its function (UniProt), what is unusual
about it, a map with a read-out of every residue (sites, charge, IDRs, flanks, UniProt regions, RNA and DNA binding in
cells, AlphaFold confidence), each charged segment with its sequence and what its modifications would do to its charge,
its sites with their sources, and its AlphaFold model. Methylarginine evidence at residue level is in the
[Methylarginine Site Explorer](https://shechterlab.github.io/PTM_Pipeline/).

## Code shared with the Methylarginine Site Explorer

The Overview's charts, scorecard tiles, tooltips, SVG/PNG/TSV exports and statistics (Mantel-Haenszel, Fisher,
hypergeometric, Benjamini-Hochberg, Mann-Whitney), its tables (sortable headers, row clicks) and the marks, facts and
folds of its protein pages are that site's own code. `tools/sync-ptm-core.mjs` copies the
declarations named in `tools/ptm-core.json`, and everything they use, verbatim from `docs/assets/explorer.js` of
Shechterlab/PTM_Pipeline at a pinned commit into `docs/assets/ptm-core.js` (and the matching rules into
`ptm-core.css`, with the rules for the classes the IDR Explorer's own markup uses, listed there too); the IDR Explorer supplies its data and selection with `PTMCore.use({...})`. Never edit the generated
files: change the explorer, update the pinned commit, and run

```bash
npm ci --prefix tools
node tools/sync-ptm-core.mjs --from ../PTM_Pipeline     # or without --from: fetch the pinned commit from GitHub
node tools/sync-ptm-core.mjs --check                   # CI: fail if the files differ from what the explorer gives
node --test tests/overview.test.mjs                    # the browser's whole-dataset tests equal the Python tables
```

## idrkit

The sequence definitions the lab's IDR work shares, so each is written once:

| | |
| --- | --- |
| `charge_vector`, `ncpr`, `fcr`, `HIS` | residue charge (histidine under the value each tool uses) |
| `delta`, `delta_z` | Das & Pappu's delta (blobs of 5 and 6, as localCIDER computes it), and its z-score against shuffles of the same residues: local charge blocks |
| `scd`, `scd_z` | sequence charge decoration, and its shuffle z-score: long-range segregation |
| `acidic_stretches`, `acidic_runs`, `basic_patches`, `rg_clusters` | the charged sequence elements |
| `e_tract_count`, `ck2_motif_count`, `rg_motif_count`, `longest_acidic_run` | the grammars the dtIDR scores count |
| `idr_intervals`, `merge_intervals` | IDRs from disorder intervals |
| `modified_charge`, `charge_shift`, `MODIFIED_CHARGE` | the charge of modified residues (phospho −2, acetyl-K 0, citrulline 0, ...) and the NCPR and delta a sequence has with its sites modified |

```bash
pip install "idrkit @ git+https://github.com/Shechterlab/idr-explorer"
python -c "import idrkit as K; print(K.delta_z('KKKKKEEEEE' * 6))"
pytest            # in a clone: the tests include localCIDER's delta and a literature SCD value
```

glu_proteomics' `glu.sequence` imports its acidic stretches, RG clusters and charge from here.

## How the data are made

`docs/data/` is generated, never edited, by `Analysis/scripts/build_idr_explorer.py` in the lab's MicroDomains
repository, with [PTM_Pipeline](https://github.com/Shechterlab/PTM_Pipeline), glu_proteomics and this repository cloned
beside it. The regions, contexts and modification kinds are defined once (MicroDomains `idr_universe.py`), the
whole-dataset statistics come from `idr_overview_stats.py` with PTM_Pipeline's implementations, and
`docs/data/build.json` records the commit of each repository and the MD5 of every input. MicroDomains and
glu_proteomics are private until publication; the inputs they read are public and listed below.

| Data | Source | Licence |
| --- | --- | --- |
| reference proteome | UniProt 2026_03 reviewed human | CC BY 4.0 |
| disorder | MobiDB (Piovesan et al. 2025) | CC BY 4.0 |
| domains | InterPro 108.0 | CC0 |
| methylarginines | PTM_Pipeline release v1.4.1 (DeAngelo, Silverstein & Shechter, Cell Chem Biol 2026) | CC BY-NC-SA 4.0 |
| phosphorylation, lysine acylation, ubiquitin/SUMO, lysine methylation | iPTMnet export (Huang et al. 2018) staged in PTM_Pipeline | CC BY-NC-SA 4.0 |
| modified residues with enzymes (CK2), ADP-ribosylation, citrullination, polyglutamylation | UniProt 2026_03 MOD_RES (PTM_Pipeline scripts/53) | CC BY 4.0 |
| acetyl, ubiquitin, SUMO and methyl sites, Gold/Silver/Bronze | PTMeXchange human builds, PRIDE PXD079218, PXD068989, PXD068489, PXD076372 (PTM_Pipeline scripts/54) | CC0 |
| citrullination | Rebak et al., Nat Struct Mol Biol 2024, doi:10.1038/s41594-024-01214-9 (PTM_Pipeline scripts/55) | CC BY 4.0 |
| ADP-ribosylation | Hendriks et al., Nat Commun 2021, doi:10.1038/s41467-021-26172-4 (scripts/56) | CC BY 4.0 |
| phosphosites and occupancy in 33 cell lines | Koenig et al., Nat Struct Mol Biol 2026, doi:10.1038/s41594-026-01877-6 (scripts/57) | CC BY 4.0 |
| acetylation stoichiometry | Hansen et al., Nat Commun 2019, doi:10.1038/s41467-019-09024-0 (scripts/58) | CC BY 4.0 |
| kinase motifs | Johnson et al., Nature 2023, doi:10.1038/s41586-022-05575-3, scored with kinase-library 1.8.0 | CC BY 4.0; CC BY-NC-SA 3.0 |
| κ, Ω, diagram-of-states classes | computed with localCIDER 0.1.21 (Holehouse et al., Biophys J 2017, doi:10.1016/j.bpj.2016.11.3200), run as a program | GPL v2 (the program) |
| glutamylation | TTLL4 re-analysis of PRIDE PXD065251 (Dey et al., J Cell Sci 2026) | lab analysis |

## Licence

Code (idrkit, the site): MIT. Data in `docs/data/`: CC BY-NC-SA 4.0, because the methylarginine release and iPTMnet
are; the sources keep their own terms (table above). See [LICENSE.md](https://github.com/Shechterlab/idr-explorer/blob/HEAD/LICENSE.md).
