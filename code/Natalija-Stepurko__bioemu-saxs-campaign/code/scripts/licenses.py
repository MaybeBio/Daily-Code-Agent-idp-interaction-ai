"""Licence register of the project's software and data.

    python3 scripts/licenses.py [--venv ID=PATH ...] [--allow-missing]

Reads the curated entries (tools, data sources, services, licence corrections and the list of
environments) from docs/licenses/curated.json and the package metadata of every environment listed
there. The metadata are read by running a short probe with that environment's own interpreter
(importlib.metadata), so nothing is imported from the packages themselves. Writes:

- docs/licenses/python-packages.csv   every installed distribution of every environment
- docs/licenses/python-packages.md    the same as a readable appendix, one section per environment
- THIRD_PARTY_LICENSES.md             the blocks between `<!-- generated:NAME:start/end -->` markers
                                      (tools, data, direct Python dependencies, packages for review);
                                      the text outside the markers is written by hand

Standard library only; any Python >= 3.11 runs it. The licence of a package is taken, in this order,
from a curated correction (with the URL and date it was checked), the License-Expression field, the
License field when it is a short identifier, the trove classifiers, and otherwise left "unknown".
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CURATED = REPO / "docs" / "licenses" / "curated.json"
OUT_CSV = REPO / "docs" / "licenses" / "python-packages.csv"
OUT_MD = REPO / "docs" / "licenses" / "python-packages.md"
REGISTER = REPO / "THIRD_PARTY_LICENSES.md"

# Run inside each environment's interpreter; prints one JSON list of distributions.
PROBE = r"""
import importlib.metadata as m, json, sys
out = []
for d in m.distributions():
    md = d.metadata
    urls = {}
    for u in md.get_all("Project-URL") or []:
        k, _, v = u.partition(",")
        urls[k.strip()] = v.strip()
    out.append({
        "name": md["Name"], "version": md["Version"],
        "license_expression": md.get("License-Expression") or "",
        "license": md.get("License") or "",
        "classifiers": [c for c in (md.get_all("Classifier") or []) if c.startswith("License ::")],
        "home_page": md.get("Home-page") or "", "project_urls": urls,
        "requires": md.get_all("Requires-Dist") or [],
        "license_files": md.get_all("License-File") or [],
    })
json.dump({"python": sys.version.split()[0], "dists": out}, sys.stdout)
"""

# Short License-field values and trove classifiers mapped to SPDX identifiers. Only mappings that are
# unambiguous; anything else is reported as found. A classifier that names no version or variant
# ("BSD License", "Apache Software License") is reported as a family with "(variant not stated)".
FIELD_TO_SPDX = {
    "mit": "MIT",
    "mit license": "MIT",
    "bsd-3-clause": "BSD-3-Clause",
    "bsd 3-clause": "BSD-3-Clause",
    "bsd 3-clause license": "BSD-3-Clause",
    "3-clause bsd license": "BSD-3-Clause",
    "bsd-2-clause": "BSD-2-Clause",
    "2-clause bsd": "BSD-2-Clause",
    "apache 2.0": "Apache-2.0",
    "apache-2.0": "Apache-2.0",
    "apache license 2.0": "Apache-2.0",
    "apache license, version 2.0": "Apache-2.0",
    "isc license": "ISC",
    "psf-2.0": "PSF-2.0",
    "mpl-2.0": "MPL-2.0",
    "mpl-2.0 and mit": "MPL-2.0 AND MIT",
    "gplv3": "GPL-3.0",
    "lgplv2.1+": "LGPL-2.1-or-later",
}
CLASSIFIER_TO_SPDX = {
    "MIT License": "MIT",
    "ISC License (ISCL)": "ISC",
    "Mozilla Public License 2.0 (MPL 2.0)": "MPL-2.0",
    "GNU General Public License v3 (GPLv3)": "GPL-3.0",
    "GNU General Public License v3 or later (GPLv3+)": "GPL-3.0-or-later",
    "GNU Lesser General Public License v3 or later (LGPLv3+)": "LGPL-3.0-or-later",
    "GNU Lesser General Public License v2 or later (LGPLv2+)": "LGPL-2.0-or-later",
    "GNU Lesser General Public License v2 (LGPLv2)": "LGPL-2.0",
    "BSD License": "BSD (variant not stated)",
    "Apache Software License": "Apache (version not stated)",
    "Python Software Foundation License": "PSF-style (classifier)",
}
COPYLEFT = re.compile(r"\b(A?GPL|LGPL|MPL)", re.I)


def norm_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def licence_of(dist: dict, corrections: dict) -> tuple[str, str]:
    """(licence, basis) for one distribution's metadata."""
    fix = corrections.get(norm_name(dist["name"]))
    if fix and dist["version"] in fix.get("versions", [dist["version"]]):
        return fix["licence"], "curated"
    if dist["license_expression"].strip():
        return dist["license_expression"].strip(), "License-Expression"
    first = dist["license"].strip().splitlines()[0].strip() if dist["license"].strip() else ""
    if first.lower() in FIELD_TO_SPDX:
        return FIELD_TO_SPDX[first.lower()], "License"
    names = [c.split(" :: ")[-1] for c in dist["classifiers"]]
    mapped = [CLASSIFIER_TO_SPDX.get(n, n) for n in names]
    if mapped:
        return " / ".join(dict.fromkeys(mapped)), "classifier"
    if first and len(first) <= 40 and first.upper() != "UNKNOWN":
        return first, "License (as written)"
    return "unknown", "none"


def review_flag(licence: str) -> str:
    if licence == "unknown":
        return "unknown"
    if COPYLEFT.search(licence):
        return "copyleft"
    return ""


def homepage(dist: dict) -> str:
    urls = {k.lower(): v for k, v in dist["project_urls"].items()}
    for key in ("homepage", "home", "source", "source code", "repository", "github", "documentation"):
        if urls.get(key):
            return urls[key]
    return dist["home_page"] or next(iter(dist["project_urls"].values()), "")


def requirement_name(spec: str) -> str:
    spec = spec.split(";")[0].strip()
    if "@" in spec:
        spec = spec.split("@")[0]
    return norm_name(re.split(r"[\s<>=!~\[(]", spec, maxsplit=1)[0])


def direct_names(env: dict) -> list[str]:
    """Direct dependencies of an environment: from pyproject.toml, a requirements file, or a list."""
    src = env["direct"]
    if src["kind"] == "pyproject":
        p = tomllib.loads((REPO / src["path"]).read_text())
        specs = list(p["project"]["dependencies"])
        for group in p.get("dependency-groups", {}).values():
            specs += [s for s in group if isinstance(s, str)]
    elif src["kind"] == "requirements":
        specs = []
        for line in (REPO / src["path"]).read_text().splitlines():
            line = line.split("#")[0].strip()
            if line and not line.startswith("-"):
                specs.append(line)
    else:
        specs = src["names"]
    return sorted({requirement_name(s) for s in specs})


def probe(python: Path) -> dict:
    res = subprocess.run([str(python), "-c", PROBE], capture_output=True, text=True, check=True)
    return json.loads(res.stdout)


def roots(dists: list[dict]) -> list[str]:
    """Installed distributions that no other installed distribution names in Requires-Dist (requirements
    that apply only with an extra are left out, since which extras were installed is not recorded)."""
    required = {
        requirement_name(r) for d in dists for r in d["requires"] if "extra" not in r.partition(";")[2]
    }
    return sorted(norm_name(d["name"]) for d in dists if norm_name(d["name"]) not in required)


def md_cell(text: str) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def table(header: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(md_cell(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def sources_cell(sources: list[dict]) -> str:
    return "; ".join(f"[{s['what']}]({s['url']}) ({s['checked']})" for s in sources)


def render_tools(cur: dict) -> str:
    rows = [
        [
            f"**{t['name']}**" + (f" {t['version']}" if t.get("version") else ""),
            t["use"],
            t["obtained"],
            t["redistributed"],
            t["licence"],
            t.get("restrictions", ""),
            t.get("attribution", ""),
            t["status"],
            sources_cell(t["sources"]),
        ]
        for t in cur["tools"]
    ]
    head = [
        "Tool",
        "How the project uses it",
        "Obtained",
        "In this repo",
        "Licence or terms",
        "Restrictions in practice",
        "Attribution",
        "Status",
        "Checked at",
    ]
    return table(head, rows)


def render_data(cur: dict) -> str:
    rows = [
        [
            f"**{d['name']}**" + (f" ({d['version']})" if d.get("version") else ""),
            d["use"],
            d["licence"],
            d.get("attribution", ""),
            d["tracked"],
            d["raw_redistributed"],
            d["status"],
            sources_cell(d["sources"]),
        ]
        for d in cur["data"]
    ]
    head = [
        "Source",
        "How the project uses it",
        "Licence or terms",
        "Attribution",
        "What the repo tracks from it",
        "Raw data redistributed",
        "Status",
        "Checked at",
    ]
    return table(head, rows)


def collect(cur: dict, venv_override: dict, allow_missing: bool) -> list[dict]:
    corrections = {norm_name(k): v for k, v in cur["package_corrections"].items()}
    envs = []
    for env in cur["environments"]:
        path = Path(venv_override.get(env["id"], env["venv"]))
        python = path / "bin" / "python"
        if not python.exists():
            if allow_missing:
                print(f"skipped {env['id']}: no interpreter at {python}", file=sys.stderr)
                continue
            sys.exit(
                f"{env['id']}: no interpreter at {python} (pass --venv {env['id']}=PATH or --allow-missing)"
            )
        info = probe(python)
        direct = direct_names(env)
        rows = []
        for d in sorted(info["dists"], key=lambda d: norm_name(d["name"])):
            lic, basis = licence_of(d, corrections)
            rows.append(
                {
                    "environment": env["id"],
                    "package": d["name"],
                    "version": d["version"],
                    "direct": "yes" if norm_name(d["name"]) in direct else "",
                    "licence": lic,
                    "basis": basis,
                    "review": review_flag(lic),
                    "homepage": homepage(d),
                    "license_expression": d["license_expression"],
                    "license_field": (d["license"].strip().splitlines() or [""])[0][:100],
                    "classifiers": "; ".join(c.replace("License :: ", "") for c in d["classifiers"]),
                }
            )
        installed = {norm_name(r["package"]) for r in rows}
        envs.append(
            {
                **env,
                "venv": str(path),
                "python": info["python"],
                "rows": rows,
                "direct_names": direct,
                "direct_missing": [n for n in direct if n not in installed],
                "roots": roots(info["dists"]),
            }
        )
    return envs


def write_csv(envs: list[dict]) -> None:
    cols = [
        "environment",
        "package",
        "version",
        "direct",
        "licence",
        "basis",
        "review",
        "homepage",
        "license_expression",
        "license_field",
        "classifiers",
    ]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, lineterminator="\n")
    w.writeheader()
    for e in envs:
        w.writerows(e["rows"])
    OUT_CSV.write_text(buf.getvalue())


def write_appendix(envs: list[dict], cur: dict) -> None:
    out = [
        "# Python packages by environment",
        "",
        "Generated by `scripts/licenses.py` from the installed package metadata of each environment;",
        "do not edit by hand. The same rows are in `python-packages.csv`. Licence basis:",
        "`License-Expression` and `License` are the package's own metadata fields, `classifier` its trove",
        "classifiers, `curated` a correction recorded with its source in `curated.json` (see the table at",
        "the end). Review: `copyleft` marks GPL, LGPL, AGPL or MPL terms, `unknown` a package whose metadata",
        "state no licence.",
        "",
        table(
            ["Environment", "Purpose", "Interpreter", "Python", "Installed", "Direct", "For review"],
            [
                [
                    e["id"],
                    e["purpose"],
                    f"`{e['venv']}`",
                    e["python"],
                    str(len(e["rows"])),
                    str(len(e["direct_names"])),
                    str(sum(1 for r in e["rows"] if r["review"])),
                ]
                for e in envs
            ],
        ),
        "",
    ]
    for e in envs:
        src = e["direct"]
        where = f"`{src['path']}`" if src["kind"] in ("pyproject", "requirements") else src["source"]
        out += [
            f"## {e['id']}",
            "",
            f"{e['purpose']}. Interpreter `{e['venv']}`, Python {e['python']}. "
            f"Direct dependencies from {where}: " + ", ".join(f"`{n}`" for n in e["direct_names"]) + ".",
            "",
        ]
        if e["direct_missing"]:
            out += [
                "Named as direct but not installed: "
                + ", ".join(f"`{n}`" for n in e["direct_missing"])
                + ".",
                "",
            ]
        out += [
            "Installed packages that no other installed package requires (a check on the direct list): "
            + ", ".join(f"`{n}`" for n in e["roots"])
            + ".",
            "",
            table(
                ["Package", "Version", "Direct", "Licence", "Basis", "Review", "Homepage"],
                [
                    [
                        r["package"],
                        r["version"],
                        r["direct"],
                        r["licence"],
                        r["basis"],
                        r["review"],
                        r["homepage"],
                    ]
                    for r in e["rows"]
                ],
            ),
            "",
        ]
    out += [
        "## Curated licence corrections",
        "",
        "Applied where a package's metadata state no licence, state it without its variant, or state one",
        "that its upstream licence file contradicts.",
        "",
        table(
            ["Package", "Licence", "Why", "Checked at"],
            [
                [
                    k + (" " + ", ".join(v["versions"]) if v.get("versions") else ""),
                    v["licence"],
                    v["note"],
                    sources_cell(v["sources"]),
                ]
                for k, v in cur["package_corrections"].items()
            ],
        ),
        "",
    ]
    OUT_MD.write_text("\n".join(out))


def render_direct(envs: list[dict]) -> str:
    rows = [
        [r["package"], e["id"], r["version"], r["licence"], r["review"], r["homepage"]]
        for e in envs
        for r in e["rows"]
        if r["direct"]
    ]
    return table(["Package", "Environment", "Version", "Licence", "Review", "Homepage"], rows)


def render_review(envs: list[dict]) -> str:
    seen: dict[tuple, list[str]] = {}
    for e in envs:
        for r in e["rows"]:
            if r["review"]:
                key = (r["package"], r["licence"], r["review"])
                seen.setdefault(key, []).append(
                    f"{e['id']} {r['version']}" + (" (direct)" if r["direct"] else "")
                )
    rows = [
        [p, lic, flag, ", ".join(where)]
        for (p, lic, flag), where in sorted(seen.items(), key=lambda x: x[0][0].lower())
    ]
    return table(["Package", "Licence", "Review", "Environments (version)"], rows)


def replace_block(text: str, name: str, body: str) -> str:
    pat = re.compile(rf"(<!-- generated:{name}:start -->\n)(?:.*?\n)?(<!-- generated:{name}:end -->)", re.S)
    if not pat.search(text):
        sys.exit(f"THIRD_PARTY_LICENSES.md has no generated:{name} block")
    return pat.sub(lambda m: m.group(1) + body + "\n" + m.group(2), text)


def curated_blocks(cur: dict) -> dict[str, str]:
    return {"tools": render_tools(cur), "data": render_data(cur)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--venv",
        action="append",
        default=[],
        metavar="ID=PATH",
        help="use PATH for the environment ID listed in curated.json",
    )
    ap.add_argument(
        "--allow-missing", action="store_true", help="skip environments whose interpreter is absent"
    )
    a = ap.parse_args()
    cur = json.loads(CURATED.read_text())
    envs = collect(cur, dict(v.split("=", 1) for v in a.venv), a.allow_missing)
    write_csv(envs)
    write_appendix(envs, cur)
    text = REGISTER.read_text()
    blocks = {
        **curated_blocks(cur),
        "python-direct": render_direct(envs),
        "python-review": render_review(envs),
    }
    for name, body in blocks.items():
        text = replace_block(text, name, body)
    REGISTER.write_text(text)
    n = sum(len(e["rows"]) for e in envs)
    print(
        f"{len(envs)} environments, {n} package rows, "
        f"{sum(1 for e in envs for r in e['rows'] if r['review'])} for review"
    )


if __name__ == "__main__":
    main()
