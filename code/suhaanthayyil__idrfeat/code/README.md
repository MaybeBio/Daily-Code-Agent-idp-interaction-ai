# idrfeat

`idrfeat` annotates intrinsically disordered regions (IDRs) in proteins and computes a wide
feature row for every IDR segment. It is the feature-extraction layer for a project on whether
missense variants perturb IDR function, using labels from a per-residue protein fitness screen.
The features cover amino acid composition, charge and its patterning, hydropathy,
k-mer content, low complexity, disorder and fold-on-binding propensity, short linear motifs,
PTM sites, and phase-separation signals.

The library is built so that the sequence-based features always work with no downloads. Disorder
predictors, annotation databases, and conservation tracks are optional: when one is missing, its
columns are filled with missing values or zero counts and a flag records that, rather than
failing the run.

## Install

The project targets Python 3.11. Using [uv](https://docs.astral.sh/uv/):

```
uv venv --python 3.11
uv pip install -e ".[dev]"
```

That installs the core library, which computes every sequence-based feature. The heavier
scientific tools are optional extras:

```
uv pip install -e ".[predictors]"     # metapredict, localCIDER
uv pip install -e ".[binding]"        # AIUPred (disorder and fold-on-binding)
uv pip install -e ".[conservation]"   # pyBigWig (phyloP, phastCons)
```

metapredict is the primary disorder predictor. When it is not installed, IDR segments are called
with a transparent built-in heuristic so the pipeline still runs. The heuristic is a smoke-test
stand-in, not a replacement for metapredict. Every output row records which backend defined its
segment in the `disorder_backend` column.

## Quickstart

Run the feature extractor on the bundled example, which needs no downloads:

```
idrfeat features --fasta src/idrfeat/data/example.fasta --out features.parquet
```

This writes `features.parquet` and `features.csv`, one row per IDR segment. Point `--fasta` at
your own protein FASTA to run it on real data.

## Annotation databases

Annotation features (ELM motifs, dbPTM sites, PhaSePro and DrLLPS phase-separation membership)
turn on when you pass the matching file. To fetch the databases:

```
python scripts/fetch_databases.py --out data/raw
```

Then pass them on the command line, for example:

```
idrfeat features --fasta proteins.fasta --out features.parquet \
    --elm data/raw/elm_instances.tsv \
    --dbptm data/raw/dbPTM_Phosphorylation.gz data/raw/dbPTM_Acetylation.gz \
    --phasepro data/raw/phasepro_full.json \
    --drllps data/raw/drllps_LLPS.txt
```

URLs, versions, and licenses for every source are recorded in `data/raw/SOURCES.txt` at fetch
time.

## Develop

```
uv run pytest          # tests for the feature logic
uv run ruff check .    # lint
uv run black --check . # format check
```

## License

MIT, see `LICENSE`. Each external database and predictor keeps its own upstream license; ELM is
academic and non-commercial.
