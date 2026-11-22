# OpenDNA from Python

For the notebook where this work usually happens.

```python
import opendna

run = opendna.predict("disorder", {"ubiquitin": "MQIFVKTLTGK..."})
run.frame.head()
```

```text
      name  position residue  probability  disordered  confidence
 ubiquitin         1       M       0.1775           0      0.1625
 ubiquitin         2       Q       0.1699           0      0.1699
```

## It brings the provenance with it

That is the reason this exists rather than a `subprocess.run` in your own
notebook. A data frame that has lost which model made it is how a prediction
becomes a fact in a spreadsheet, and then in a slide, with nobody able to say
later what produced it.

```python
run.build                      # 'f7d99273e4b9b961' -- the exact weights
run.manifest["input"]["sha256"]  # what went in
run.limitations                # the model card's caveats, in full
print(run.cite())
```

```text
disorder predicted with OpenDNA head version 1.1.0, build f7d99273e4b9b961,
on frozen facebook/esm2_t30_150M_UR50D embeddings (CPU)

Limitations stated by the model card:
  - DisProt records disorder that has been experimentally demonstrated. ...
```

## Installing

It needs the `opendna` binary, which is where the predictions actually happen.
Take a build from the [releases
page](https://github.com/dyber-pqc/OpenDNA/releases) and put it on your path,
or point `OPENDNA_BIN` at it.

```bash
pip install -e bindings/python          # from a checkout
pip install -e 'bindings/python[frames]'  # with pandas, for .frame
```

Predictions also need the trained heads, which are not bundled:

```bash
hf download dyber-pqc/opendna-disorder --local-dir models
export OPENDNA_MODELS=$PWD/models
```

## The four heads

```python
opendna.predict("disorder",   {"name": sequence})   # per residue
opendna.predict("burial",     {"name": sequence})   # per residue
opendna.predict("solubility", {"name": sequence})   # per protein
opendna.predict("stability",  {"name": sequence},
                mutations=[("name", "L8P"), ("name", "I44D")])
```

`stability` scores named substitutions and needs to be told which, so it takes
`mutations`. Positions are 1-based, as they are written. Every regressor
returns `interval_low` and `interval_high` beside its point estimate; for a
prediction near zero the interval will span zero, which is the head saying it
cannot tell you the direction.

## Nothing here estimates anything

A head that is not installed raises `HeadMissing`. It does not fall back to a
heuristic, because a plausible number from no model is worse than an error.
The exceptions mirror the command-line tool's exit codes:

| | |
|---|---|
| `ToolMissing` | no `opendna` binary |
| `HeadMissing` | the head is absent or does not verify (exit 3) |
| `InputRejected` | the sequences or mutations are not usable (exit 2 or 4) |
| `OpenDnaError` | anything else, including a failed self-check (exit 5) |

## Why it shells out

The command-line tool is the path with the exit-code contract, the weight
verification and the run manifest already tested around it. Shelling out
inherits all three. A compiled extension would mean a wheel per platform and a
second code path to keep honest, for the same four predictions — worth doing
only if the per-call overhead ever matters, and for a run that spends minutes
embedding, it does not.

The cost is one process per call, so batch your sequences into one `predict`
rather than calling it in a loop. That is faster anyway: the heads sort
sequences by length and embed them together.

## Tests

```bash
pip install -e 'bindings/python[test]'
pytest bindings/python
```

The refusal tests need only the binary. The prediction tests need the heads and
skip with a message when they are absent.
