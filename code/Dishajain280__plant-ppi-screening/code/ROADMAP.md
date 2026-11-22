# Status and roadmap

## Where this stands

**Built and verified:** the full harness. Symbol resolution against UniProt, token
budgeting and batch generation within the server's quota, parsing of summary
confidence files, ranking metrics with bootstrap intervals, a dependency-free
figure, and 35 tests that run in CI on every push.

**Not done:** the predictions. No AlphaFold job has been submitted, so this
repository contains **no biological result**. That is stated in the README and it
is deliberate — a number produced from unverified labels would be worse than no
number.

**Why it is paused:** the two remaining steps both need a human. Labels require
reading papers and signing off on each one; the AlphaFold Server has no API, so
jobs are uploaded by hand at 30 per day.

## Next steps, with honest effort estimates

| # | Step | Time | Blocked by |
| --- | --- | --- | --- |
| 1 | Verify 10 pairs following [docs/curating-labels.md](https://github.com/Dishajain280/plant-ppi-screening/blob/HEAD/docs/curating-labels.md) | ~2 h | reading the cited papers |
| 2 | Fill `range_a`/`range_b` for the BRI1–BAK1 calibration pair from InterPro | 20 min | choosing defensible domain boundaries |
| 3 | Run the calibration pair, check the full-length versus domain gap | 20 min | 1 day of server quota |
| 4 | Submit the benchmark pairs and download the results | ~1 day, mostly waiting | 30 jobs/day quota |
| 5 | `ingest` → `evaluate` → update this file with the measured result | 1 h | steps 1–4 |

Total focused effort: roughly one evening plus one day of waiting.

## Definition of done

- At least 10 verified benchmark pairs, committed **before** any prediction is seen.
- Calibration reproduced: full-length BRI1–BAK1 scoring clearly below its truncated
  domains.
- A reported AUC **with its interval**, recall in the top decile, and decoy scores
  kept separate from curated negatives.
- The README figure replaced with real scores, and the demo figure labelled as demo.

## Non-goals

- Automating the web interface to bypass the daily quota. It breaches the server's
  terms and would cost access mid-project.
- Reporting metrics on labels nobody verified. The evaluator refuses this, and that
  guard stays.
- Predicting before the label set is frozen in a commit.

## Success looks like this

A small, honest number — even an unimpressive one. On plant proteins, where
alignments are sparse and disorder is common, "the ranking was no better than
chance on 12 verified pairs" is a legitimate finding and is the exact question the
project asks. A suspiciously clean AUC on a set of easy decoys is not.
