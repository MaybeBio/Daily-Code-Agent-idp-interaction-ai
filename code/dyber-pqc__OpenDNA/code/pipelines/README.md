# Running OpenDNA under a workflow engine

Modules for [Nextflow](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/pipelines/nextflow) and [Snakemake](https://github.com/dyber-pqc/OpenDNA/blob/HEAD/pipelines/snakemake), for the job
OpenDNA is actually for: a few hundred candidate sequences to triage, spread
across a cluster.

They are deliberately thin. The CLI already answers the questions a scheduler
cannot — which models are installed, whether they verify, what each failure
means, whether two parts were produced the same way — and these wire that up
rather than reimplementing it.

```bash
nextflow run pipelines/nextflow/main.nf \
    --fasta candidates.fasta \
    --models /mnt/shared/opendna-models \
    --lockfile /mnt/shared/opendna-models/models.lock \
    --task disorder --parts 20
```

```bash
snakemake --snakefile pipelines/snakemake/Snakefile --cores 8 \
    --config fasta=candidates.fasta models=/mnt/shared/opendna-models \
             task=disorder parts=20
```

The models directory is one produced by `opendna stage`: heads, embedding
model, lockfile and checksums in one place. Put it on shared storage the
compute nodes can read. Nothing in the run touches the network.

## What they do that a shell loop does not

**Fail before the queue fills.** `models --check` and `selftest` run once,
before anything else is scheduled, and the whole run stops if the installed
heads are not the ones in the lockfile or if the machine does not reproduce
their recorded results. Finding that out after four hundred tasks is finding
it out too late.

**Retry only what retrying can fix.** Exit 3 means the installation is wrong,
4 means the data is wrong, 5 means the install does not reproduce itself. None
of those will be different the second time, so only exit 1 — an unanticipated
failure, which might be a full disk or an evicted node — is retried.

**Refuse a mixture.** The merge step is `opendna merge`, not `cat`. It reads
each part's run manifest and stops if the parts were produced with a different
head build, on a different device, or at a different batch size. On a cluster
where nodes are upgraded one at a time, half a run landing on machines with a
GPU is a real thing that happens, and the resulting file is a mixture nobody
can describe.

**Balance the parts by work.** `opendna split --parts N` divides by residue
count rather than sequence count, because runtime scales with residues and a
scheduler waits for the slowest part. On a 200-sequence test file, four parts
came out at 8852, 8735, 8740 and 8706 residues — holding 50, 52, 50 and 48
sequences.

## Why `--batch-tokens 0`

Both pipelines set it, and it costs about 10% of the CPU runtime.

Sequences are sorted by length and batched within whatever file they arrive
in, so dividing the input differently gives different batch-mates, different
padding, and a different order of floating-point sums. The part count is a
scheduling decision — somebody will change it from 20 to 50 because the queue
is busy — and it should not be able to change a number in a results table.

Measured on 200 sequences and 35,033 residues, comparing a whole-file run
against a four-way split:

| | |
|---|---|
| `--batch-tokens 0` | **identical**, all 35,034 rows |
| `--batch-tokens 8192` | 23 of 35,033 probabilities differ, no call flipped |

The differences are in the fourth decimal. None of them flipped a disorder
call in that run, but one could for a residue sitting within 1e-4 of the
threshold, so the guarantee at 0 is worth the 10%. Raise it if you would
rather have the time and your results do not need to be independent of how the
work was divided — the batch size is recorded in every run manifest either
way.

## Stage the bundle before running parts in parallel

Both pipelines take a models directory. Point it at a directory produced by
`opendna stage`, not at a bare heads directory, and pass `--offline`.

The reason is measured rather than stylistic: a run with no staged bundle
loads the embedding model through the Hugging Face cache, and two processes
cannot read that cache at once — on Windows the second is refused the weights
file and exits 3. A pipeline whose whole purpose is to run parts side by side
will hit this on the first node that gets two tasks. Three concurrent runs
against a staged bundle were fine.

## Checking it yourself

```bash
pipelines/verify.sh /mnt/shared/opendna-models candidates.fasta
```

Runs what the pipelines run — check, split, predict each part, merge — and
compares the merged result against predicting the whole file in one go. It
also relabels one part's manifest as though it had run on another device and
confirms the merge refuses it. Exit 0 means chunking does not change the
answers on your hardware.

`BATCH_TOKENS=8192 pipelines/verify.sh ...` shows the contrary case, and
explains itself when it fails.

## What has not been checked

**Neither pipeline has been run against its engine.** Nextflow and Snakemake
are not installed on the machine these were written on, so the `.nf` and
`Snakefile` are unverified as syntax and unverified as workflows. What *is*
verified is the behaviour they depend on, by `verify.sh` above: that splitting
and merging reproduces a whole-file run exactly, and that a mixture is
refused. Treat the engine files as a starting point to adapt, not as something
known to run.

Neither pipeline uses `opendna predict --checkpoint`, which would let a part
that failed halfway continue instead of rerunning whole. The CLI supports it
and it is verified to produce byte-identical results across an interruption,
but wiring it in here is not a one-line change: both engines give a retried
task a fresh working directory, so the checkpoint has to live somewhere that
outlives the attempt, and that choice depends on the cluster's filesystem
rather than on anything testable here. What the pipelines get for free is
Snakemake rebuilding only what is missing, which works because `opendna split`
is deterministic — the same input and part count give the same chunks byte for
byte.

Chunking would at least be numerically free in these pipelines: they set
`--batch-tokens 0`, where every sequence is its own batch and how the work is
divided cannot change a value. Measured, not assumed.
