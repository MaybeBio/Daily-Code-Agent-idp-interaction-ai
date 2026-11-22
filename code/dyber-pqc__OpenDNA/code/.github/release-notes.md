## Read this first

**This is a rewrite, not an upgrade.** OpenDNA is now pure Rust. The Python
backend, the React interface, ESMFold structure prediction, ESM-IF1 inverse
folding, the molecular dynamics tools, the antibody and plasmid tools and the
model hub have all been removed. If you depend on any of those, stay on
`v0.5.0-rc1`, the last release that has them.
[The changelog](https://github.com/dyber-pqc/OpenDNA/blob/main/docs/CHANGELOG.md)
lists the removals in full.

**This is a release candidate.** The tag that produced these files is the first
thing the release workflow has ever built, no installer has previously been
produced on macOS or Linux, and none of them are signed, so your operating
system will warn you. The command-line tool is the part that has been measured
and tested.

## What is in the files below

| File | What it is |
|---|---|
| `opendna-*-x86_64-unknown-linux-gnu.tar.gz` | command-line tool, Linux |
| `opendna-*-aarch64-apple-darwin.tar.gz` | command-line tool, Apple Silicon |
| `opendna-*-x86_64-apple-darwin.tar.gz` | command-line tool, Intel Mac |
| `opendna-*-x86_64-pc-windows-msvc.zip` | command-line tool, Windows |
| `*-setup.exe`, `*.msi` | desktop application, Windows |
| `*.dmg` | desktop application, macOS |
| `*.AppImage`, `*.deb`, `*.rpm` | desktop application, Linux |
| `SHA256SUMS` | checksums for everything above |

Check what you downloaded before running it:

```bash
sha256sum --check --ignore-missing SHA256SUMS
```

A file missing from the list is a platform whose build failed. That is visible
on purpose rather than hidden by a release that refuses to publish at all.

## Predictions need the trained heads

Nothing ships with model weights. Download the heads you want and point
`OPENDNA_MODELS` at the directory:

```bash
hf download dyber-pqc/opendna-disorder --local-dir models
hf download dyber-pqc/opendna-burial --local-dir models
hf download dyber-pqc/opendna-stability --local-dir models
hf download dyber-pqc/opendna-solubility --local-dir models
export OPENDNA_MODELS=$PWD/models
```

If one of those reports that the repository does not exist, that head has not
been published yet.

Then prove the installation reproduces the numbers it is supposed to:

```bash
opendna selftest
```

Eighteen checks, nine of them values recorded against a specific build. This
is worth running rather than skipping: `opendna models` passes an install
whose weights have been swapped for another head's, because the files are
intact and verify — it is only running them that catches it.

The first prediction also downloads the ESM-2 embedding model, a few hundred
megabytes, after which it works offline. For a machine with no network, or for
running several predictions at once, stage a bundle first:

```bash
opendna stage --to /mnt/shared/opendna
opendna --models /mnt/shared/opendna --offline predict --task disorder \
    -i candidates.fasta -o disorder.csv
```

Staging is not only for air-gapped machines: two runs cannot read the Hugging
Face cache at the same time, so anything running parts in parallel needs a
bundle.

Without a head installed, OpenDNA reports the prediction as unavailable. It
does not estimate the value some other way.

## How good the models are

| Head | In-domain test | External |
|---|---|---|
| Disorder | AUROC 0.883 | **AUROC 0.941** on CAID3 Disorder-PDB |
| Accessibility | Spearman 0.737 | not benchmarked externally |
| Stability (ΔΔG) | Spearman 0.674 | Spearman 0.442 over 4 ProteinGym assays |
| Solubility | Spearman 0.701 | none exists under usable terms |

Disorder is the strong one. Stability is a ranking signal and its in-domain
figure overstates generalisation by about 0.23 Spearman — treat an individual
ΔΔG as weak, and read the interval printed beside it, which for a prediction
near zero will span zero.

Each head's accuracy, baseline and limitations are on its model card, and
every result file carries those limitations as comments. Read them before
acting on a number. [The validation page](https://github.com/dyber-pqc/OpenDNA/blob/main/docs/VALIDATION.md)
is written for whoever has to sign this off, and says what has not been
checked.

## Two licence obligations worth knowing

The solubility head is the one to watch. It is trained on eSOL, which is
CC BY-SA 2.1 JP, and because that licence is share-alike the weights carry the
same terms rather than the Apache-2.0 of the other three. `opendna stage`
writes a `LICENSES.md` that groups the heads by licence and flags that one,
because a result passed on without it is an obligation somebody does not know
they have.

The binaries contain two MPL-2.0 components, `colored` and `option-ext`, which
arrive through the Hugging Face client used to fetch the embedding model.
Neither is modified; their source is the published crate at the version in
`Cargo.lock`. If you redistribute these binaries, that obligation travels with
them. The command-line README has the detail.
