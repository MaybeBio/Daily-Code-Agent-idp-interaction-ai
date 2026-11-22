# Paper drafts

**These drafts describe the pre-rewrite system and are kept as history, not as
current claims.** Read them that way or not at all.

`preprint.md` (April 2026) and `preprint_v050.md` (April 2026) were written
against the Python service that has since been removed. Specifically:

- They describe structure prediction, inverse folding, multimer prediction,
  molecular dynamics and "~50 sequence- and structure-based analyses". The
  Rust application does none of that. Most of what they describe was
  implemented as heuristics behind a fallback — `predictors.py` reimplemented
  TMHMM, SignalP, TANGO and Camsol as hydropathy arithmetic, and `docking.py`
  returned an affinity in kcal/mol from a pocket estimator — which is the
  specific practice the rewrite exists to remove.
- They state the licence as "Apache 2.0 with Commons Clause". The project is
  Apache-2.0; that grant appears in no `LICENSE` file and in no crate's
  metadata. See the [licence section of the README](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/README.md#licence).
- `preprint_v050.md` leads on post-quantum cryptography, which is not a
  property of the current application.
- The benchmarking claims rest on a TM-score *proxy* computed from RMSD.
  `BENCHMARKS.md` in this directory is careful about that and says "TM-score
  proxy" throughout; a `docs/BENCHMARKS.md` page that dropped the qualifier
  and tabled it against real CASP15 medians has been deleted.

Nothing here has been edited to match the current system. Rewriting a dated
draft in place would destroy the record of what was claimed when, and these
files are the record.

A paper describing what actually exists would be a different document: three
heads on frozen ESM-2 embeddings, a licence-audited dataset pipeline, and
external benchmark numbers with the protocol stated — including a negative
result on a larger embedding model. The honest version of that is in
[training/benchmarks](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/training/benchmarks).
