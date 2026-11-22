# metadata

Small text files that describe the data. Committed to Git.

- `sources.tsv`: where each raw file comes from (URL, database version, license).
- `domains.tsv`: the CATH domain set (Phase 1), one row per domain.
- `idp_proteins.tsv`: proteins for the disorder extension (Phase 5), one row per protein.
- `checksums.sha256`: written by the download scripts; check with
  `shasum -a 256 -c metadata/checksums.sha256`.

## Data dictionary
(Add a column table for each file when it is created.)
