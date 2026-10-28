# One-command regeneration. `make all` runs every stage after the download; `make page` rebuilds
# the project page from the tracked results/. BSC_ROOT moves data and results elsewhere;
# BSC_DENSS_ENV points at the DENSS environment (numpy < 2); ARCHIVE at Predictions.tar.gz for coordrg.
# `coverage` reads results/coverage/, written by analysis/coverage/run_all.sh (network, extra requirements).
UV ?= uv run
ARCHIVE ?=
export BSC_N_JOBS ?= 6

.PHONY: all fetch score coordrg analyse select_eval predict envelopes coverage report page test lint licenses

all: score coordrg analyse select_eval predict envelopes coverage report page

fetch:
	$(UV) bsc fetch

score:
	$(UV) bsc score

coordrg:
	$(UV) bsc coordrg $(if $(ARCHIVE),-- --archive $(ARCHIVE),)

analyse:
	$(UV) bsc analyse

select_eval:
	$(UV) bsc select_eval

predict:
	$(UV) bsc predict

envelopes:
	$(UV) bsc envelopes

coverage:
	$(UV) bsc coverage

report:
	$(UV) bsc report

page:
	$(UV) python docs/site/build.py

test:
	$(UV) pytest -q

lint:
	$(UV) ruff check src tests docs scripts

# licence register: THIRD_PARTY_LICENSES.md and docs/licenses/ from curated.json and the environments
licenses:
	python3 scripts/licenses.py
