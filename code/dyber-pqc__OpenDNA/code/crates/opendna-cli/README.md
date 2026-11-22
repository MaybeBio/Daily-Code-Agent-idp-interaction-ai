# opendna-cli

The `opendna` command: runs the trained heads over a whole FASTA file.

The desktop application works one protein at a time, which is the wrong shape
for the job people actually have — a few hundred candidates to triage. This
runs the same prediction code, so the two cannot disagree about a number.

```bash
opendna models
opendna predict --task disorder    -i candidates.fasta -o disorder.csv
opendna predict --task solubility  -i candidates.fasta -o solubility.csv
opendna predict --task stability   -i candidates.fasta --mutations muts.csv -o ddg.csv
opendna predict --task burial      -i candidates.fasta -o accessibility.csv
```

`burial` gives relative solvent accessibility per residue — roughly 0 for a
residue in the core, roughly 1 for one fully exposed — which is what answers
where a mutation is tolerable and where a linker can go. Read it beside
`disorder`: accessibility describes a folded structure, so for a region with
no stable fold the number is still produced and still means nothing.

`--device cpu` forces the CPU on a build with GPU support compiled in, which
is what you want when the GPU is busy with someone else's job and what makes
it possible to check the two agree.

`--format` takes `csv`, `tsv`, `json`, or `caid` for a CAID submission — see
[CAID.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-cli/CAID.md). With no `-o`, the table goes to standard output and
progress goes to standard error, so a pipe stays clean.

Sequences are read as written rather than validated against the twenty
standard amino acids. Real target files contain selenocysteine, pyrrolysine
and the ambiguity codes — CAID's own `DP03745` has a `U` at position 48 — and
ESM-2 has distinct tokens for all of them, so collapsing them to `X` would
feed the model a residue the file does not specify.

Heads are found via `--models`, then `$OPENDNA_MODELS`, then `./models`:

```bash
hf download dyber-pqc/opendna-disorder --local-dir models
export OPENDNA_MODELS=$PWD/models
```

## Knowing what you ran

A result is only defensible if the model behind it can be named and checked.

```bash
opendna models                       # what is installed, and does it verify
opendna models --lock models.lock    # pin it
opendna models --check models.lock   # fail unless the install is still that
```

Each head records a semantic version, a build id digesting everything that
determined the weights (the dataset manifest, the hyperparameters, the
embedding model, the training code), and the SHA-256 of the weights file
beside it. `models` rehashes the weights and reports one of three states:

| | |
|---|---|
| `verified` | the weights are the file the metadata names |
| `NOT VERIFIABLE` | metadata predating hashes; usable, but unprovable |
| `CORRUPT` | the weights are not that file |

Corrupt weights make `models` exit non-zero, and `predict` refuses before
loading anything. A half-finished download otherwise produces numbers from a
file nobody can identify, which is worse than no numbers.

`--check` exits non-zero when the installed heads differ from the lockfile,
including a retrain that kept the same version number — the build id catches
that. Commit the lockfile beside the analysis it produced.

## Running on a machine with no network

The heads are about a megabyte each and copy anywhere. The embedding model does
not: it is 568 MB, it lives in the Hugging Face cache rather than beside the
heads, and a machine that has never run a prediction downloads it on its first
one. On an instrument computer or anything behind an air gap that is not a slow
first run, it is a failure at the moment somebody needs a result.

```bash
opendna stage --to /mnt/usb/opendna-models        # on a connected machine
opendna --models ./opendna-models --offline selftest
```

`stage` writes one directory holding the four heads, their cards, the embedding
model's own files, `models.lock`, and a `sha256sums.txt` in the format
`sha256sum -c` reads. It copies what this machine already resolved rather than
re-fetching on the target, so the bundle is exactly what was tested here. It
refuses to stage an installation that does not verify — a bundle nobody can
trust, carried to a machine where nobody can re-download, is the worst possible
place to discover a bad file.

`--offline` is a hard promise, not a preference: nothing will reach the
network, and if the embedding model is not staged the command says so and names
the two files it wanted. Without the flag a missing model is downloaded, which
is right on a laptop and wrong on anything validated.

`opendna models` says which of the two situations a directory is in.

Verified by staging a bundle, pointing `HF_HOME` at an empty directory so the
cache held nothing, and running `--offline`: the self-test passed all
seventeen checks and a disorder run gave values identical to the online one.

## Proving the install works here

Verifying the files is not the same as checking the machine. A head can hash
correctly and still produce nonsense because this computer resolved a different
BLAS, because the GPU path is broken, or because the embedding model in the
cache is not the one the head was trained against.

```bash
opendna selftest
```

It runs all four heads over ubiquitin and checks seventeen things, going
through the same code `predict` does rather than a parallel harness. Eight are
invariants that hold for any build of these heads — the output has one value
per residue, probabilities stay in `[0, 1]`, the hydrophobic core reads as less
accessible than the chain-forming lysines, `L8P` is more destabilising than
`L8M`, and a second identical run gives identical values. Nine are numbers
measured from named builds, with a tolerance.

A pinned number is only checked against the build it was measured from. Against
any other build it reports `not checked`, naming both builds, because nobody
has measured what that head does on this fixture. A retrain is therefore
expected to leave those unchecked until someone records new values; the
invariants still run.

Exit status is 0 when nothing failed and non-zero when anything did, so it goes
straight into an installation script. It catches what `models` cannot: swapping
the burial weights into the disorder head and rewriting the metadata to match
passes `opendna models` completely clean and fails three of the pinned checks.

The tolerances are a judgement, not a measurement — the recorded values come
from one machine, on CPU. A failure that lands just outside one is worth
reporting rather than assuming the bound is too tight.

## Disorder probabilities are calibrated

The disorder head's probabilities are corrected before they are reported, and
the threshold for a call comes from the head rather than from a constant.

Measured on the held-out split: expected calibration error 0.217 before the
correction, 0.032 after. Where the raw head said 0.55, 13% of those residues
were actually annotated disordered. The correction is monotonic, so AUROC is
unchanged; what changes is the meaning of the number.

That moves the threshold. The conventional 0.5 was never right for this head —
on its raw scale the best threshold is 0.73 — and on the calibrated scale it is
0.34, chosen on the validation split. Both the map and the threshold live in
`disorder.metadata.json`, so `predict`, the desktop application and the CAID
output all use the same ones, and recalibrating the head changes the call
without changing any code.

A head carrying no recorded calibration falls back to 0.5 and its card says
calibration was not measured. "Measured and left alone" and "never measured"
are different states and the card distinguishes them.

## Regressor outputs carry an interval

`stability`, `solubility` and `burial` report `interval_low` and
`interval_high` beside their point estimate, and a `#` comment saying what
confidence they carry and what coverage that confidence actually achieved.

```
# interval_low/high: 90% prediction interval, +/-1.327 kcal/mol, which covered
# 93.7% of held-out cases. ...
name,mutation,ddg_kcal_per_mol,log_odds,interval_low,interval_high
ubiquitin,L8P,-0.779,-4.878,-2.1061,0.5481
ubiquitin,I44D,-1.589,-12.891,-2.9161,-0.2619
```

L8P's interval spans zero: at 90% confidence the head cannot say whether that
mutation destabilises. I44D's does not. That distinction is the useful
content, and it is invisible in the point estimate alone.

The intervals come from split conformal prediction, fitted on the validation
split and measured on the held-out one. The coverage actually achieved is
reported beside the level asked for, because they differ — the splits are
clustered by sequence identity, which is the right way to build them and
breaks the exchangeability a conformal guarantee assumes. The burial head
undercovers at every level. See each model card, and
[docs/VALIDATION.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/VALIDATION.md).

A head with no recorded intervals gets no interval columns rather than
invented ones.

## Not every residue deserves equal trust

Every disorder output carries a `confidence` column: how far the calibrated
probability sits from the decision threshold. Measured on held-out data, it
orders residues by how likely the call is to be wrong.

```bash
opendna predict --task disorder -i candidates.fasta --min-confidence 0.25 -o out.csv
```

That adds a `low_confidence` column marking the rows below the cut-off, and
prints what the cut-off was measured to buy:

```
Marking calls below confidence 0.25. Measured at the nearest recorded cut-off
(0.25) on held-out data: 57% of residues clear it and 5.3% of those calls were
wrong, against 15.7% over all residues.
```

Rows are marked, never removed. A table quietly missing rows looks exactly like
a complete one. The cut-off is yours because the right one depends on what a
wrong call costs you; what the head supplies is the measured consequence of
each choice, recorded in its metadata.

The signal captures about half the error reduction an oracle would achieve, so
it is worth having and is not a guarantee. The one out-of-domain proxy measured
alongside it — whether a residue's 5-mer appeared in training — carried no
information about error, and the card records that so nobody re-tests it.

## Stability needs to be told which mutations

The head scores named substitutions rather than scanning everything, so it
takes a CSV of `name,mutation` where the name matches the FASTA header and the
mutation is written the usual way, 1-based:

```csv
name,mutation
ubiquitin,L8P
ubiquitin,K6A
```

The wild-type letter in each label is checked against the sequence at that
position, and a disagreement stops the run. That check exists because a
mutation list numbered against a slightly different construct — a tagged
version, a different isoform — otherwise scores a neighbouring residue and
produces a table that looks completely normal.

## Every output says what produced it

Delimited output carries `#` comment lines with the head version, the
embedding model, the device, and the limitations taken from the model card.
JSON carries the same as fields. Readers skip `#` by convention: pandas with
`comment='#'`, R with `comment.char='#'`.

The limitations are read out of the card rather than written here, so the CLI
cannot state a different set than the card does. This matters more than it
looks: a prediction that travels without them becomes a number in a
spreadsheet, then a number in a slide, and by then nobody can tell which model
made it or what it does not cover.

## Every run writes down what made it

The `#` comments say which head produced a file, which is what a person
reading it needs. It is not enough to reproduce the run: which input, which
bytes of that input, which build of the head, which embedding model, on what,
for how long. So every run that writes an output file writes a manifest beside
it.

```
out.csv
out.manifest.json
```

It is JSON because something else has to read it — a pipeline collecting
provenance, a notebook, a script checking two runs used the same model. One run
writes one file and never appends, so there is no race to worry about.

| | |
|---|---|
| `input` | path, SHA-256 of the bytes read, sequence and residue counts |
| `output` | path, SHA-256 of the bytes written, row count |
| `head` | task, version, build id, weights SHA-256, path |
| `embedding_model` | repo, where it came from, size, digest if known |
| `run` | UTC start, wall-clock seconds, device, the arguments as given |
| `limitations` | the model card's, so the manifest alone states the caveats |

The input is hashed from the bytes in memory rather than by re-reading the
file, so the digest describes what was actually used even when the input came
through a pipe and has no path at all.

The embedding model is the one thing not hashed by default: it is 568 MB and
digesting it takes about as long as a small run. A staged bundle already wrote
the digest into its `sha256sums.txt`, so for a staged install the manifest
carries it for free; otherwise the field is null and `sha256_source` says why.
An absent field with a reason is honest — a missing field somebody assumes was
checked is not.

`--no-manifest` turns it off, for the case where something downstream cannot
tolerate a second file appearing.

## Speed, and what it costs

Sequences are embedded in batches rather than one at a time. Every member of a
batch is padded to the longest, so the grouping sorts by length first and then
fills each batch to `--batch-tokens` padded tokens — which is also what bounds
the memory a batch needs.

Measured on 200 sequences and 35,033 residues, interleaving the settings across
rounds on an otherwise idle machine:

| device | `--batch-tokens 0` | `--batch-tokens 8192` |
|---|---|---|
| CPU | 60.8 s (59.1–63.8) | **54.5 s** (51.7–55.5) |
| RTX 2000 Ada | **6.4 s** (6.3–6.4) | 6.8 s (6.8–6.8) |

Batching helps the CPU by about a tenth and costs the GPU about a sixteenth,
which was the opposite of what we expected on the GPU: its arithmetic is
already wide enough that the padded work costs more than the per-call overhead
it saves. So one number cannot be right for both, and left unset
`--batch-tokens` follows the device. Your GPU may differ — the flag is there so
you can measure yours rather than inherit ours.

If you benchmark this, interleave the settings and report the spread. Two
earlier attempts here ran one setting after the other while a compile competed
for the cores, and produced first an overstated gain and then no detectable
one. A third, cleaner run of the same comparison put the CPU gain nearer a
sixth where this one says a tenth, so treat the size as approximate and the
direction as settled.

Sequences longer than the embedding model's limit cannot be batched — a
windowed embedding is assembled from several passes over the same sequence — so
they run one at a time and the progress line says so.

## How reproducible the numbers are

Exactly as reproducible as this, and no more:

| | |
|---|---|
| Same input, same device, same `--batch-tokens` | **byte-identical**, verified on both CPU and GPU |
| Same input, same device, different `--batch-tokens` | differs in the last digit of ~0.2% of values |
| Same input, interrupted and resumed at the same `--chunk` | **byte-identical** |
| Same input, different `--chunk`, above `--batch-tokens 0` | differs in the last digit of ~0.1% of values |
| Same input, CPU vs GPU | differs in the last digit of ~0.3% of values |

Across devices, 102 of 35,033 disorder probabilities differed, by at most
0.0001 and by 0.0000003 on average — the f32 rounding floor, not a disagreement
about biology. Nothing in those runs changes which residues are called
disordered at any sensible threshold. But it does mean a result file from the
GPU and one from the CPU are not the same bytes, so do not diff them and do not
compare their checksums.

Batch size and device are both recorded in every run manifest, which is what
makes this checkable rather than a claim: two manifests that agree on input
digest, head build, device and batch size should have produced the same output
digest, and if they did not, something is wrong that is worth finding.

`opendna selftest` is the tool for this. Its pinned values were recorded on
CPU; running it on a GPU checks that the GPU reproduces them inside the
fixture's tolerances, and on the hardware here all seventeen checks pass on
both.

The one thing to avoid: driving one model from several threads at once on
CUDA changes the low bits, because cuBLAS picks differently shaped reductions
under concurrent load. The CLI is serial, so this does not arise here, but it
is why it is serial.

## Running it across a cluster

`split` cuts a FASTA into parts balanced by residue count, because runtime
scales with residues and a scheduler waits for the slowest part. `merge` joins
the results back, reading each part's run manifest and refusing to join parts
produced with a different head build, on a different device, or at a different
batch size — which `cat` cannot notice and a cluster upgraded one node at a
time will eventually produce.

```bash
opendna stage --to /mnt/shared/opendna
opendna split -i candidates.fasta --out-dir parts --parts 20
# ... one task per part, each with --models /mnt/shared/opendna --offline ...
opendna merge --output disorder.csv parts/part-*.csv
```

**Stage the bundle first if parts will run at the same time.** This is not an
optimisation. Two runs that both load the embedding model through the Hugging
Face cache cannot do it at once: on Windows the second is refused the weights
file by the operating system and exits 3 saying the installation needs
fixing, which is unhelpful but not wrong — what needs fixing is that the
model is being read out of a cache instead of a bundle. Measured here: two
concurrent runs through the cache, one failed; three concurrent runs against a
staged directory, all three fine. A staged bundle is opened read-only and
shares without trouble, which is also why `--offline` is the right flag for a
cluster task.

Nextflow and Snakemake modules are in [pipelines/](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/pipelines), along with
the measurement showing that splitting and merging reproduces a whole-file run
exactly at `--batch-tokens 0`, and does not above it.

## Resuming a run that was interrupted

A run over fifty thousand sequences takes hours, and the things that stop it
are ordinary: a scheduler's wall-clock limit, a power cut, somebody else
needing the GPU. `--checkpoint` makes that survivable.

```bash
opendna predict --task disorder -i candidates.fasta -o disorder.csv \
  --checkpoint disorder.ckpt --chunk 64
```

Results are recorded in `disorder.ckpt` as they are produced and flushed to
disk after every chunk, so a process that is killed loses at most one chunk.
Give the same command again and it continues from where it stopped. The
results file is written only when the run finishes: the checkpoint is a record
of work done, not a results file.

**Being interrupted does not change the answer.** That is the point, and it is
the reason `--chunk` exists rather than the work being resumed at whatever
sequence it happened to reach. Sequences are batched for embedding within the
chunk they fall in, and batch composition changes the order of floating-point
sums, so a run cut at an arbitrary point would not match an uninterrupted one.
Cutting at the same boundaries every time removes the question. Verified on a
twelve-sequence file: a run resumed from a checkpoint whose last record had
been truncated mid-write produced the same 854 rows, byte for byte, as a run
that was never interrupted — the same output digest in both manifests.

What it costs is that chunking is itself a setting that changes the numbers,
in the same way and to the same degree as `--batch-tokens` does. Measured on
the same file, against a run with no checkpoint at all:

| | |
|---|---|
| `--chunk` anything, with `--batch-tokens 0` | **identical** |
| `--chunk 1` or `4`, at the CPU default of 8192 | 1 of 853 accessibility values differs, by 0.0001 |

At `--batch-tokens 0` every sequence is its own batch, so how the file is
divided cannot matter and chunking is free. Above it, the chunk size joins the
list of things a manifest has to record, and it does: a run that used a
checkpoint records its chunk size, a run that did not records nothing there,
and `opendna merge` treats those as different runs and refuses to join them.

Resuming is refused, naming the field that differs, if the input bytes, the
head build, the head's weights digest, the embedding model, the device, the
batch size, the mutation file or the confidence cut-off is not the one the
checkpoint was written with — or if the chunk size differs. A table that is
half one run and half another is the mixture `merge` exists to refuse between
files, and creating it inside one file would be worse, not better.

Two more things it does deliberately:

- A file that is not an OpenDNA checkpoint is refused and left alone. A
  mistyped `--checkpoint` must not destroy what it was pointed at.
- A truncated last record is dropped with a message saying so, and the file is
  cut back to its last complete record before anything is appended. An
  interrupted write costs one chunk, not the run.

## What the exit status means

A pipeline cannot read an error message. It reads a number, and what it needs
from that number is what to do next.

| | | |
|---|---|---|
| 0 | it worked | — |
| 1 | an unanticipated failure | report it |
| 2 | the command asked for something that does not exist | fix the command |
| 3 | the installation is wrong | fix the machine; retrying will not help |
| 4 | the data is wrong | fix the input; retrying will not help |
| 5 | the checks ran and did not hold | do not use results from this machine |

2 is deliberately the same code `clap` uses for a usage error, since a mistyped
flag and an impossible combination of valid ones want the same response.
3 covers a missing or corrupt head, a lockfile that does not match, and
`--offline` with nothing staged. 4 covers unparseable FASTA, a mutation naming
a residue the sequence does not have, a sequence longer than the embedding
model takes, and an antibody with no 9-mer `humanness` can score. 5 comes from
`selftest`; from `fold` when an engine exits cleanly having returned a
structure that is not the sequence it was sent; and from `humanness` when the
threshold asked for is finer than the table can express.

A folding engine that fails outright is 1, deliberately. The cause could be the
installation (a missing CUDA library), the input (a sequence longer than the
engine's memory allows) or the engine itself, and nothing here can tell those
apart; filing it under a category a caller would act on would be guessing at
the action. What the engine said is printed instead.

Anything not classified is 1. That is the honest default: a failure nobody
anticipated is not one a caller can be told how to handle, and filing it under
a known category would be worse than admitting it.

Every failure prints to standard error and names which of the above it is.
`predict` streams results to standard output, so a failing run leaves standard
output empty — a pipe collecting results never receives an error as though it
were data.

Checks that need only the arguments and the input happen before any model
loads, so a typo in a mutation file is refused in about 50 ms rather than after
568 MB has loaded for a run that cannot work.

The contract is exercised against the compiled binary in
[tests/exit_codes.rs](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/crates/opendna-cli/tests/exit_codes.rs). It needs the heads installed and
skips with a message when they are not, because a test that fails on every
machine without them gets disabled and then checks nothing at all.

## Folding a sequence

`opendna fold` predicts structures by driving an installed ESMFold or
ColabFold. OpenDNA does not fold proteins and has nothing to fall back on: with
no engine installed this exits 3 naming what it looked for, because coordinates
from something that cannot fold proteins would be worse than no coordinates.

```bash
opendna fold -i candidates.fasta -o structures/
```

```text
3 sequences from candidates.fasta through esmfold (esmfold 2.0.0)
Running: /opt/conda/bin/esm-fold -i /work/structures/engine/input.fasta -o /work/structures/engine/out
name                      residues       mean B  structure
ubiquitin                       76         91.4  structures/ubiquitin.pdb
p53_dbd                        198         88.0  structures/p53_dbd.pdb
designed_07                     64         54.2  structures/designed_07.pdb
```

One `.pdb` per sequence, named after its FASTA header, plus a
`.manifest.json` beside each naming the engine that produced it, the version
it reported, the digest of the file that ran and the command line it was
given. Everything the engine wrote is kept under `structures/engine/`, so the
things OpenDNA does not read — ColabFold's per-model scores, its alignments,
its log — are still there.

**A structure is matched to its sequence by its residues, not by its
filename.** Every file in the engine's output directory is parsed, and the one
whose chain is the sequence that was sent is the one written out. This is the
part that is easy to skip and expensive to skip: a folding run writes outputs
under names of its own choosing, several per sequence in ColabFold's case, and
the ordinary way a batch goes quietly wrong is a structure filed against the
wrong sequence — a stale file from a previous run, a name mangled on the way
through, a failure that produced fewer outputs than inputs without saying so.
An engine that returns a structure of something else exits 5 rather than
writing it.

Where that leaves ties, they are broken in the open: ColabFold's ranked models
all have the right sequence, so `rank_001` wins; two inputs that are the same
sequence under different names are separated by name; and anything still
ambiguous is refused with the candidates listed, rather than resolved by
directory order.

The `mean B` column is the B-factor column averaged over each chain's Cα
atoms, which is where both engines write their per-residue confidence (pLDDT).
It is reported because it is useful and named for what it is: nothing here
checked what a given engine meant by that column, and the number is the
engine's estimate of its own output, not a measurement.

Anything the engine should be told beyond the input and the output goes
through `--engine-arg`, verbatim and repeatable:

```bash
opendna fold -i candidates.fasta -o structures/ --engine-arg --num-recycles --engine-arg 4
```

Nothing is translated between engines. A flag whose meaning OpenDNA invented
would be a flag that quietly means something else after the engine's next
release, and the full command line is in every manifest the run writes.

With both engines installed, `--engine esmfold` or `--engine colabfold` is
required rather than optional: they do not produce the same structures, and
choosing for you would leave the result depending on something nobody wrote
down. `--check engines.lock` refuses to start unless the installed engines are
the pinned ones.

What is not here: multimers, which need the engines' own chain-separator
conventions and a different check than "one chain, this sequence"; and any use
of the predicted coordinates beyond writing them out. The desktop application
opens a PDB file, so a folded model can be viewed; nothing yet feeds one back
into a head.

## How human an antibody looks

`opendna humanness` cuts a heavy-chain sequence into overlapping 9-residue
peptides and reports, for each one, the fraction of human donors whose antibody
repertoire contains that exact peptide. The measure follows OASis (Prihoda et
al. 2022).

```bash
opendna humanness -i candidate.fasta -t humanness/table.bin
```

No model runs. Given the table, every number is an exact count, so there is no
accuracy to quote and no interval to report: the table decides what the numbers
mean, and its card says what it was built from. Without a table this exits 3 —
a score with no denominator is not a weaker score, it is not a score.

The per-window detail is the output worth reading, and it is what this prints
by default. Three real heavy chains from the antibody cohorts — one human, one
C57BL/6 mouse, one from a Kymouse carrying human immunoglobulin loci — scored
against the 141-donor table:

```text
$ opendna humanness -i example.fasta -t datasets/out/humanness/table.bin
table: 141 donors, 1096088 peptides, from datasets/out/humanness/table.bin
name                     identity  windows   skipped  rarest window
human-IGHM                  79.8%      109         0  YCARVPPFP at 79 (0.00% of donors)
mouse-C57BL6                24.3%      111         0  KESGPGLVA at 4 (0.00% of donors)
kymouse-transgenic          90.3%      113         0  YYCASSPLE at 93 (0.00% of donors)

human-IGHM: 1 stretch below the medium (>=10% of donors)
       residues  windows    lowest  rarest peptide in the stretch
         78-108       22     0.00%  YCARVPPFP

mouse-C57BL6: 6 stretches below the medium (>=10% of donors)
       residues  windows    lowest  rarest peptide in the stretch
           4-25       13     0.00%  KESGPGLVA
          18-28        2     0.00%  SITCTVSGF
          22-37        7     0.00%  TVSGFSLTS
          31-41        2     0.00%  YGVSWVRQP
          43-95       44     0.00%  GLEWLGVIW
         89-113       16     0.00%  ATYYCAKTD

kymouse-transgenic: 1 stretch below the medium (>=10% of donors)
       residues  windows    lowest  rarest peptide in the stretch
         92-111       11     0.00%  YYCASSPLE
```

That is the shape of the answer. The human antibody and the humanised one each
have exactly one unusual stretch, their CDR3, which is unique to an individual
by construction and is *supposed* to be rare. The mouse antibody is unusual in
six places spread through its framework, which is where a humanisation campaign
would have to work. A single percentage cannot tell those two situations apart;
the stretch list does, and it names the residues.

Stretches are reported rather than the *n* rarest windows because rare windows
cluster: nine overlapping windows over one odd CDR3 are one finding, and a list
of the rarest came back as the same stretch nine times over, every entry at the
same floor.

A single percentage also invites comparison with a published figure computed
from other repertoires, and that comparison is wrong — prevalence depends
entirely on which donors are in the table and how deeply each was sequenced.

### What it was measured to do

Scored against the four antibody evaluation cohorts, whose donors are disjoint
from the table's. Median per-sequence identity at the 10% threshold:

| Cohort | median identity | the right answer |
|---|---:|---|
| human | 70.7% | high — donors the table has never seen |
| mouse | 25.9% | lower — mouse germline |
| diagnostic (mice with human Ig loci) | 90.0% | high, and it is |
| out of distribution | 36.4% | lower; rhesus 62.5%, camel 23.7% |

The diagnostic cohort is the one worth the trouble. Those are mice, and the
measure calls their antibodies human, so it tracks the immunoglobulin rather
than the animal. That survives the obvious confound: the diagnostic cohort is
mostly IGHM while the human cohort spans isotypes, and a class-switched
sequence has been mutated away from germline, but matched isotype by isotype
the diagnostic cohort still scores higher on all three.

It scores *higher* than human because transgenic repertoires are less
affinity-matured and so sit nearer germline, which is where prevalence is
highest — not because transgenic mice make better human antibodies. The same
effect limits how any single score can be read: within the human cohort alone,
median identity at the 50% threshold runs from 41.1% (IGHA) to 75.8% (IGHM), a
35-point spread between isotypes of one species. An absolute identity cannot be
read as "how human" without knowing how mutated the sequence is, which is the
third reason the per-stretch output is the one to use.

`datasets/out/humanness/EVALUATION.md` holds the full tables, regenerated by
`opendna-data evaluate humanness`.

Three thresholds are available for the summary identity, the ones OASis
defines: `--stringency loose` (1% of donors), `medium` (10%, the default) and
`strict` (50%). A threshold finer than the table can express is **refused**
rather than reported, with exit 5: from *n* donors the smallest non-zero
prevalence is 1/*n*, so "1% of donors" against a 20-donor table really means
"at least one donor", and printing it as the former would state a number
nothing measured.

Two more refusals rather than plausible answers. A sequence with no scorable
window — every 9-mer containing an `X` or a gap — exits 4 instead of scoring
zero, because zero asserts the sequence looks unlike a human antibody and
nothing here measured that. And a table whose entries are unsorted, or which
claims a peptide occurs in more donors than it has, is rejected on reading:
the first would make lookups silently miss and read low, the second would
produce prevalences above one.

What a high prevalence does not mean: that an antibody carrying the peptide is
safe, non-immunogenic or developable. Immunogenicity depends on MHC
presentation, dose, route and the recipient, none of which this measures. What
a low prevalence does not mean: that no human has the peptide — only that it
was not in these donors' sampled sequences.

The table is built by the dataset pipeline, which is where the donor choice and
its consequences are documented:

```bash
cd datasets && cargo run --release -- build humanness
```

## External engines, when there are any

Some of what this field wants is decades of other people's software: folding a
sequence, docking a ligand, running a simulation. OpenDNA drives those rather
than reimplementing them badly, and `opendna engines` is how an installation
is checked before anything depends on it.

```bash
opendna engines
```

```text
External engines OpenDNA can drive:

  esmfold     not installed
                structure prediction from sequence -- esmfold is not installed. ...
  vina        AutoDock Vina v1.2.5  3f1a9c2e5b7d0411  /usr/local/bin/vina
                ligand docking -- executable
```

Each one found is identified by the version string it prints, verbatim, and
the SHA-256 of the file that would be run. The version is not parsed into a
number: third-party tools word it however they like and change the wording
between releases, so a manifest recording what the tool actually said can be
checked by eye, while a parsed number can be wrong without anyone noticing.

`--lock` and `--check` work as they do for `models`:

```bash
opendna engines --lock engines.lock
opendna engines --check engines.lock
```

An engine that is absent is simply not written down, rather than pinned to its
absence — a lockfile saying "vina was not installed" would start failing the
day somebody installs it.

**What a digest here covers is worth knowing.** For a self-contained program
it identifies what runs. For an entry point — a script, a launcher, a Python
interpreter that imports a library — it identifies that file and not the
environment behind it, and the version string carries most of the identity.
The listing and the lockfile both say which case each engine is in, rather
than letting a reader assume the stronger one.

**Two of the four are driven.** `opendna fold` uses esmfold or colabfold, as
described above. Vina and openmm are listed so an installation can be checked
and pinned before the commands that use them exist; the order those arrive in
is in [docs/ROADMAP.md](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/docs/ROADMAP.md).

## Licences travel with the weights

Three heads are Apache-2.0. The solubility head is not: it is trained on eSOL,
which is CC BY-SA 2.1 JP, and because that licence is share-alike — and
because whether trained weights are a derivative work of their training data
is an open legal question this project will not answer in its own favour — the
weights carry the same terms.

That matters when you pass a staged bundle on, which is a redistribution.
`opendna stage` writes a `LICENSES.md` grouping the heads by licence and saying
plainly when redistribution carries obligations, and each head's model card
travels beside its weights with the attribution its licence requires. The run
manifest records the licence of the head behind each result, so a number does
not outlive the knowledge that passing it on has conditions.

All of it is read out of the cards rather than written into the program, so the
notice cannot claim terms the file beside the weights does not. A head whose
card declares no licence is reported as `unstated` and treated as possibly
obliging — defaulting it to the permissive licence the others use would be a
false statement about a file somebody is about to redistribute.

## Third-party code in the binary

Most of what `opendna` is built from is MIT or Apache-2.0. Two components are
not, and both reach the compiled binary:

| Crate | Licence | Why it is here |
|---|---|---|
| [`colored`](https://crates.io/crates/colored) | MPL-2.0 | via `hf-hub` -> `hf-xet` -> `xet-runtime` |
| [`option-ext`](https://crates.io/crates/option-ext) | MPL-2.0 | via `hf-hub` -> `hf-xet` -> `xet-runtime` -> `shellexpand` -> `dirs` |

The MPL is file-level copyleft: distributing a binary containing those files
obliges whoever distributes it to say so and to make their source available.
Neither crate is modified here, so their source is the published crate at the
version in `Cargo.lock`, and this table is the notice. OpenDNA's own code stays
Apache-2.0.

They are not optional. `hf-hub` 1.0 has no feature that turns off Xet support,
and `hf-hub` is how the embedding model is fetched. A build that never needs to
fetch it still links them.

`cargo deny check licenses` enforces the whole list, and `deny.toml` records why
each non-permissive licence was accepted rather than quietly widening the
allowed set.

## What it will not do

There is no fallback. A missing head is an error naming the files it looked
for and the command to fetch one, and a missing folding engine is an error
naming what it looked for rather than coordinates from something else. A sequence longer than the embedding model's
limit is refused for solubility, because that head pools over a single pass
and pooling across windows is a different quantity than the one it was trained
on — disorder windows instead, and says so in the output.

Read each head's model card before acting on a number. The stability head's
honest out-of-domain correlation is about 0.44 Spearman, which makes it a
ranking signal and not a measurement.
