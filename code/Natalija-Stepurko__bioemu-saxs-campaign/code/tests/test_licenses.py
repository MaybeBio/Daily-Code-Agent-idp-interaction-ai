"""scripts/licenses.py: the curated register is complete and in sync with THIRD_PARTY_LICENSES.md, and the
licence of a package is read from its metadata in the documented order."""

import datetime
import importlib.util
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("licenses", REPO / "scripts" / "licenses.py")
lic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lic)
CUR = json.loads(lic.CURATED.read_text())


def dist(**kw):
    base = {"name": "pkg", "version": "1.0", "license_expression": "", "license": "", "classifiers": []}
    return {**base, **kw}


def test_every_curated_entry_names_its_source_and_date():
    entries = CUR["tools"] + CUR["data"] + list(CUR["package_corrections"].values())
    for e in entries:
        assert e["licence"].strip()
        assert e["sources"], e
        for s in e["sources"]:
            assert s["url"].startswith("https://") or (REPO / s["url"]).exists(), s
            datetime.date.fromisoformat(s["checked"])
    for e in CUR["tools"] + CUR["data"]:
        assert e["status"].split(" ")[0] in {"confirmed", "not"} or e["status"] == "not applicable", e[
            "status"
        ]


def test_register_blocks_match_curated_file():
    text = lic.REGISTER.read_text()
    for name, body in lic.curated_blocks(CUR).items():
        m = re.search(rf"<!-- generated:{name}:start -->\n(.*?)\n<!-- generated:{name}:end -->", text, re.S)
        assert m, name
        assert m.group(1) == body, f"{name} block is stale: run python3 scripts/licenses.py"


def test_licence_order_expression_field_classifier():
    assert lic.licence_of(dist(license_expression="MIT", license="BSD"), {}) == ("MIT", "License-Expression")
    assert lic.licence_of(dist(license="BSD 3-Clause License\n\nCopyright ..."), {}) == (
        "BSD-3-Clause",
        "License",
    )
    full_text = dist(
        license="Copyright (c) 2005, X\nAll rights reserved.",
        classifiers=["License :: OSI Approved :: BSD License"],
    )
    assert lic.licence_of(full_text, {}) == ("BSD (variant not stated)", "classifier")
    assert lic.licence_of(dist(license="UNKNOWN"), {}) == ("unknown", "none")
    fix = {"pkg": {"licence": "Apache-2.0", "versions": ["1.0"]}}
    assert lic.licence_of(dist(license="MIT"), fix) == ("Apache-2.0", "curated")
    assert lic.licence_of(dist(version="2.0", license="MIT"), fix) == ("MIT", "License")


def test_review_flags():
    assert lic.review_flag("LGPL-2.1-or-later") == "copyleft"
    assert lic.review_flag("MPL-2.0 AND MIT") == "copyleft"
    assert lic.review_flag("GPL-3.0-only") == "copyleft"
    assert lic.review_flag("unknown") == "unknown"
    assert lic.review_flag("BSD-3-Clause") == ""
    assert lic.review_flag("Apache-2.0 AND Apache-2.0 WITH LLVM-exception") == ""


def test_direct_dependencies_from_project_files():
    envs = {e["id"]: e for e in CUR["environments"]}
    main = lic.direct_names(envs["main"])
    assert {"numpy", "mdtraj", "pytest", "ruff"} <= set(main)
    sim = lic.direct_names(envs["simulation"])
    assert "calvados" in sim and "openmm" in sim
    assert lic.requirement_name("calvados @ git+https://github.com/x/y.git@abc") == "calvados"
    assert lic.requirement_name("Jinja2==3.1.6") == "jinja2"
    assert lic.requirement_name('torch>=2.6.0; extra == "md"') == "torch"
