#!/usr/bin/env python3
"""
backfill_gene_registry.py — STEP 1 of the ground-truth gene re-attribution.

Using each long-read dataset's OWN gene assignment (BC = GTF gene_id ENSG;
ccRCC = SQANTI associated_gene) and the FULL GENCODE v47 annotation as the gene
registry, this script:

  1. Populates genomic coords (chrom/strand/start/end) on EXISTING genes from the
     annotation (many were empty — also fixes the bloated-locus display).
  2. Creates gene rows for every ENSG a long-read transcript is assigned to that is
     not yet in the DB (pseudogenes / lncRNA / etc.), consolidating any symbol-only
     stub gene onto its real ENSG.
  3. Writes data/regene/targets.tsv = (isoform_id, current_gene_id, target_gene_id)
     — the exact re-attribution map that backfill_regene_from_ground_truth.py applies.

Ground truth join is 100% for both sources. READ-ONLY unless --apply.

Usage:
    python backend/scripts/backfill_gene_registry.py            # dry-run
    python backend/scripts/backfill_gene_registry.py --apply
"""
import os
import sys
import csv
import gzip
import re
import argparse
from collections import defaultdict
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.models.bio import Gene, Isoform, Source

csv.field_size_limit(1 << 24)
import os
# Project root that holds references/ and datasets/ (the parent of db/<repo>);
# override with SPLISOFORMS_PROJECT_ROOT when the checkout lives elsewhere.
PROJECT_ROOT = Path(os.environ.get("SPLISOFORMS_PROJECT_ROOT") or Path(__file__).resolve().parents[4])
SP = Path(__file__).resolve().parents[2]
SOURCES = SP / "data" / "sources"
# Raw GTF ground-truth inputs. Paths were reorganised in the 2026-06 data-layout
# cleanup (gencode → references/, BC → datasets/); a couple of candidate locations
# are probed so a future move fails loudly (assertion below) rather than silently.
def _first_existing(*cands: Path, what: str) -> Path:
    for c in cands:
        if c.exists():
            return c
    raise SystemExit(
        f"backfill_gene_registry: could not locate {what}. Tried:\n  "
        + "\n  ".join(str(c) for c in cands)
        + "\nUpdate the path in backend/scripts/backfill_gene_registry.py."
    )

ANNOT_GTF = _first_existing(
    PROJECT_ROOT /"references" / "gencode_v47" / "gencode.v47.annotation.gtf.gz",
    PROJECT_ROOT /"gencode_v47" / "gencode.v47.annotation.gtf.gz",
    what="GENCODE v47 annotation GTF",
)
BC_GTF_DIR = _first_existing(
    PROJECT_ROOT /"datasets" / "bc" / "19666157",
    PROJECT_ROOT /"lr_breast_cancer_data" / "19666157",
    what="Breast-cancer assembly GTF dir (gene_id ENSG ground truth)",
)
TARGETS_OUT = SP / "data" / "regene" / "targets.tsv"

_ensg = re.compile(r'gene_id "([^"]+)"')
_gname = re.compile(r'gene_name "([^"]+)"')
_gtype = re.compile(r'gene_type "([^"]+)"')
_txid = re.compile(r'transcript_id "([^"]+)"')
strip = lambda e: e.split(".")[0] if e else e
_parent_ensg = re.compile(r'(ENSG\d+)')


def build_registry():
    reg, name_to_ensg = {}, {}
    with gzip.open(ANNOT_GTF, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = _ensg.search(f[8])
            if not m:
                continue
            e = strip(m.group(1))
            nm = (_gname.search(f[8]).group(1) if _gname.search(f[8]) else "")
            bt = (_gtype.search(f[8]).group(1) if _gtype.search(f[8]) else "")
            reg[e] = {"name": nm, "biotype": bt, "chrom": f[0],
                      "strand": f[6], "start": int(f[3]), "end": int(f[4])}
            if nm:
                name_to_ensg.setdefault(nm, e)
    return reg, name_to_ensg


def bc_orig_to_ensg():
    o2e = {}
    for g in BC_GTF_DIR.glob("*.gtf"):
        with open(g) as fh:
            for ln in fh:
                if "\ttranscript\t" not in ln:
                    continue
                t, e = _txid.search(ln), _ensg.search(ln)
                if t and e:
                    o2e[t.group(1)] = strip(e.group(1))
    return o2e


def ccrcc_assoc():
    a = {}
    for f in (SOURCES / "ccrcc" / "sqanti3").glob("*/scisoseq_transcripts.sorted.filtered_lite_classification.txt"):
        with open(f) as fh:
            gi = fh.readline().rstrip("\n").split("\t").index("associated_gene")
            for ln in fh:
                c = ln.rstrip("\n").split("\t")
                if len(c) > gi:
                    a[c[0]] = c[gi].strip()
    return a


def main():
    ap = argparse.ArgumentParser(description="Gene registry + target map for re-attribution")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    db = SessionLocal()
    print("building registry / loading ground truth …", flush=True)
    reg, name_to_ensg = build_registry()
    o2e = bc_orig_to_ensg()
    bc_tx_to_orig = {r["transcript_id"].strip(): r["orig_transcript_id"].strip()
                     for r in csv.DictReader(open(SOURCES / "breast_cancer" / "metadata.tsv"), delimiter="\t")}
    ccrcc = ccrcc_assoc()
    print(f"  registry={len(reg)}  bc_gtf={len(o2e)}  ccrcc_assoc={len(ccrcc)}", flush=True)

    # DB gene lookups (mutated as we create)
    gene_rows = db.query(Gene.id, Gene.symbol, Gene.ensembl_gene_id).all()
    gene_by_ensg = {strip(e): gid for gid, s, e in gene_rows if e}
    gene_by_sym = {s: gid for gid, s, e in gene_rows}
    created_syms = set()
    stat = defaultdict(int)

    # ── 1. coords on existing genes ──────────────────────────────────────────
    if args.apply:
        n = 0
        for gid, sym, e in gene_rows:
            if e and strip(e) in reg:
                m = reg[strip(e)]
                db.query(Gene).filter(Gene.id == gid).update(
                    {"chromosome": m["chrom"], "strand": m["strand"],
                     "genomic_start": m["start"], "genomic_end": m["end"]},
                    synchronize_session=False)
                n += 1
        db.commit()
        stat["coords_set"] = n
    else:
        stat["coords_set"] = sum(1 for _, _, e in gene_rows if e and strip(e) in reg)

    # ── gene get-or-create helpers ───────────────────────────────────────────
    def gene_for_ensg(e):
        if e in gene_by_ensg:
            return gene_by_ensg[e]
        meta = reg.get(e)
        name = meta["name"] if (meta and meta["name"]) else e
        # consolidate a symbol-only stub onto this ENSG
        if name in gene_by_sym:
            gid = gene_by_sym[name]
            if args.apply:
                upd = {"ensembl_gene_id": e}
                if meta:
                    upd.update(chromosome=meta["chrom"], strand=meta["strand"],
                               genomic_start=meta["start"], genomic_end=meta["end"],
                               description=meta["biotype"])
                db.query(Gene).filter(Gene.id == gid).update(upd, synchronize_session=False)
            gene_by_ensg[e] = gid
            stat["consolidated"] += 1
            return gid
        sym = name if (name not in gene_by_sym and name not in created_syms) else e
        gid = None
        if args.apply:
            g = Gene(symbol=sym, ensembl_gene_id=e,
                     chromosome=meta["chrom"] if meta else None,
                     strand=meta["strand"] if meta else None,
                     genomic_start=meta["start"] if meta else None,
                     genomic_end=meta["end"] if meta else None,
                     description=meta["biotype"] if meta else None)
            db.add(g); db.flush(); gid = g.id
        gene_by_ensg[e] = gid
        gene_by_sym[sym] = gid; created_syms.add(sym)
        stat["created_ensg"] += 1
        return gid

    def gene_for_novel(ng):
        if ng in gene_by_sym:
            return gene_by_sym[ng]
        pe = _parent_ensg.search(ng)
        chrom = reg.get(strip(pe.group(1)), {}).get("chrom") if pe else None
        gid = None
        if args.apply:
            g = Gene(symbol=ng, chromosome=chrom, description="novel_locus (SQANTI)")
            db.add(g); db.flush(); gid = g.id
        gene_by_sym[ng] = gid; created_syms.add(ng)
        stat["created_novelgene"] += 1
        return gid

    def resolve(is_bc, tx):
        """→ ('gid', gid) | ('unresolved', reason)"""
        if is_bc:
            e = o2e.get(bc_tx_to_orig.get(tx, ""))
            if not e:
                return ("unresolved", "bc_no_ensg")
            if e not in reg:
                return ("unresolved", "ensg_not_in_v47")
            return ("gid", gene_for_ensg(e))
        ag = ccrcc.get(tx)
        if not ag:
            return ("unresolved", "ccrcc_no_assoc")
        if ag.startswith("ENSG"):
            e = strip(ag)
            if e not in reg:
                return ("unresolved", "ensg_not_in_v47")
            return ("gid", gene_for_ensg(e))
        if ag.startswith("novelGene"):
            return ("gid", gene_for_novel(ag))
        e = name_to_ensg.get(ag)
        if e and e in reg:
            return ("gid", gene_for_ensg(e))
        if ag in gene_by_sym:
            return ("gid", gene_by_sym[ag])
        return ("unresolved", "symbol_not_found")

    # ── 2/3. resolve every novel isoform, build target map ───────────────────
    novel = db.query(Isoform.id, Isoform.transcript_id, Isoform.gene_id, Source.name)\
              .join(Source, Source.id == Isoform.source_id)\
              .filter(Isoform.ensembl_transcript_id.is_(None)).all()
    targets = []  # (iso_id, cur_gid, target_gid)
    for iso_id, tx, cur_gid, src in novel:
        is_bc = "Breast" in (src or "")
        kind, val = resolve(is_bc, tx)
        if kind == "unresolved":
            stat[f"unresolved:{val}"] += 1
            continue
        targets.append((iso_id, cur_gid, val))
        stat["changed" if val != cur_gid else "unchanged"] += 1

    if args.apply:
        db.commit()
        TARGETS_OUT.parent.mkdir(parents=True, exist_ok=True)
        with open(TARGETS_OUT, "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["isoform_id", "current_gene_id", "target_gene_id"])
            for iso_id, cur, tgt in targets:
                if tgt is not None:
                    w.writerow([iso_id, cur, tgt])

    print("\n═══ GENE REGISTRY / TARGET MAP ═══════════════════════════")
    print(f"existing-gene coords set     : {stat['coords_set']}")
    print(f"new genes created (ENSG)     : {stat['created_ensg']}")
    print(f"  consolidated symbol-stubs  : {stat['consolidated']}")
    print(f"new novelGene loci created   : {stat['created_novelgene']}")
    print(f"isoforms resolved            : {stat['unchanged'] + stat['changed']}")
    print(f"  changed (will re-home)     : {stat['changed']}")
    print(f"  unchanged                  : {stat['unchanged']}")
    for k in sorted(stat):
        if k.startswith("unresolved:"):
            print(f"  {k:26s}: {stat[k]}")
    print("\n" + ("wrote " + str(TARGETS_OUT) if args.apply else "(dry-run — re-run with --apply)"))
    db.close()


if __name__ == "__main__":
    main()
