#!/usr/bin/env python3
"""
organize_papers.py — PDF Article Organizer

Files journal articles from a downloads folder into a flat library named
`firstauthor-year.pdf`, and applies Finder tags.

Identity (author + year) comes from the article's DOI, resolved against
Crossref (or arXiv), not from guessing at the PDF text. Every resolution is
then verified by checking that the resolved *title actually appears in the
document* — without that gate, any PDF that merely cites a DOI (a grant
application, a referee report, a set of research notes) resolves to a real
paper and gets filed under its name.

Usage:
    # Always dry run first
    python organize_papers.py --dry-run --limit 20 --verbose

    # Full run, 8 workers
    python organize_papers.py --jobs 8

    # No network: fall back to text heuristics (much less accurate)
    python organize_papers.py --offline
"""

import argparse
import csv
import importlib.util
import json
import plistlib
import re
import shutil
import subprocess
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Lock

# pypdf is used by the extraction subprocess, not by this process directly,
# so check availability here rather than importing it.
if importlib.util.find_spec("pypdf") is None:
    sys.exit("pypdf not found. Run: pip install pypdf")


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

KNOWN_TAGS = {
    "IDP": "2",
    "LLPS": "2",
    "Peptides": "2",
    "Metals": "2",
    "Cyclic": "2",
    "NMR": "3",
    "FRET": "3",
    "Evolution": "3",
    "Mechanics": "3",
    "Allostery": "4",
    "Binding": "4",
    "Bioinformatics": "4",
    "Ensembles": "4",
    "Entropy": "4",
    "Enzymes": "4",
    "Heme": "4",
    "Hydrogenase": "4",
    "Interactions": "4",
    "Polymer": "4",
    "QM": "4",
    "Thermodynamics": "4",
    "Water": "4",
    "CoarseGrained": "5",
    "Design": "5",
    "ForceFields": "5",
    "FreeEnergy": "5",
    "Kinetics": "5",
    "ML": "5",
    "MSM": "5",
    "Simulation": "5",
    "COVID19": "6",
    "DNA/RNA": "6",
    "Aggregation": "7",
    "Folding": "7",
    "Mutation": "7",
    "Structure": "7",
    "TP/Friction": "7",
    "Downhill": "1",
    "Smfs": "0",
}

# Keywords for tag assignment. Checked against title + abstract text (case-insensitive).
TAG_KEYWORDS = {
    "IDP": [
        "intrinsically disordered", "intrinsically unstructured", "idp", "idr",
        "disordered protein", "disordered region", "unfolded protein", "natively unfolded",
    ],
    "LLPS": [
        "liquid-liquid phase separation", "llps", "phase separation", "condensate",
        "biomolecular condensate", "membraneless organelle", "coacervate",
    ],
    "Peptides": [
        "peptide", "dipeptide", "tripeptide", "oligopeptide", "cyclic peptide",
        "antimicrobial peptide", "amp",
    ],
    "Metals": [
        "metal ion", "metalloprotein", "zinc", "copper", "iron", "calcium",
        "magnesium", "manganese", "cobalt", "nickel", "metalloenzyme",
    ],
    "Cyclic": [
        "cyclic peptide", "cyclization", "cyclic dinucleotide",
    ],
    "NMR": [
        " nmr ", "nuclear magnetic resonance", "chemical shift", "relaxation dispersion",
        "noe ", "noesy", "hsqc", "residual dipolar coupling", "rdc",
    ],
    "FRET": [
        "fret", "förster resonance", "forster resonance", "single-molecule fret",
        "smfret", "fluorescence resonance",
    ],
    "Evolution": [
        "evolution", "evolutionary", "phylogenetic", "natural selection",
        "conservation", "sequence conservation", "homolog",
    ],
    "Mechanics": [
        "atomic force microscopy", "afm", "mechanical unfolding", "force spectroscopy",
        "stretching", "elasticity", "mechanical stability", "optical tweezer",
    ],
    "Allostery": [
        "allostery", "allosteric", "conformational change", "effector",
    ],
    "Binding": [
        "binding affinity", "dissociation constant", "kd ", "kon ", "koff",
        "protein-protein interaction", "ligand binding", "docking",
    ],
    "Bioinformatics": [
        "bioinformatics", "sequence alignment", "database", "genome", "proteome",
        "blast", "hmm", "sequence analysis", "deep learning", "neural network",
    ],
    "Ensembles": [
        "conformational ensemble", "structural ensemble", "ensemble",
        "population", "heterogeneous", "conformational sampling",
    ],
    "Entropy": [
        "entropy", "entropic", "free energy landscape", "thermodynamic",
        "configurational entropy",
    ],
    "Enzymes": [
        "enzyme", "catalysis", "catalytic", "substrate", "active site",
        "michaelis", "turnover", "enzymatic",
    ],
    "Heme": [
        "heme", "haeme", "cytochrome", "hemoglobin", "myoglobin", "porphyrin",
    ],
    "Hydrogenase": [
        "hydrogenase", "hydrogen evolution", "h2 production",
    ],
    "Interactions": [
        "protein interaction", "molecular interaction", "electrostatic",
        "van der waals", "hydrophobic interaction",
    ],
    "Polymer": [
        "polymer", "polyelectrolyte", "polymer physics", "worm-like chain",
        "freely jointed chain", "persistence length",
    ],
    "QM": [
        "quantum mechanics", "qm/mm", "density functional", "dft",
        "ab initio", "quantum chemical", "hartree-fock",
    ],
    "Thermodynamics": [
        "thermodynamics", "enthalpy", "heat capacity", "calorimetry",
        "differential scanning", "isothermal titration",
    ],
    "Water": [
        "water molecule", "hydration", "solvation", "water dynamics",
        "hydrophobic effect", "water network",
    ],
    "CoarseGrained": [
        "coarse-grained", "coarse grained", "cg model", "martini",
        "go model", "go-like", "elastic network",
    ],
    "Design": [
        "protein design", "de novo design", "computational design",
        "directed evolution", "engineering",
    ],
    "ForceFields": [
        "force field", "amber ff", "charmm36", "opls", "gromos",
        "force-field", "parameterization", "forcefield",
    ],
    "FreeEnergy": [
        "free energy", "potential of mean force", "pmf", "umbrella sampling",
        "free energy perturbation", "fep", "thermodynamic integration",
        "metadynamics", "replica exchange",
    ],
    "Kinetics": [
        "kinetics", "rate constant", "rate coefficient", "folding rate",
        "unfolding rate", "transition state", "kramers", "diffusion",
    ],
    "ML": [
        "machine learning", "deep learning", "neural network", "alphafold",
        "transformer", "graph neural", "random forest", "support vector",
    ],
    "MSM": [
        "markov state model", "msm", "transition network", "kinetic network",
        "committor", "metzner",
    ],
    "Simulation": [
        "molecular dynamics", "monte carlo simulation", "md simulation",
        "amber", "gromacs", "namd", "openmm", "charmm",
        "trajectory", "force field",
    ],
    "COVID19": [
        "covid", "sars-cov", "coronavirus", "spike protein", "ace2",
        "pandemic", "viral",
    ],
    "DNA/RNA": [
        "dna", "rna", "nucleic acid", "nucleotide", "base pair",
        "double helix", "ribosome", "transcript",
    ],
    "Aggregation": [
        "aggregation", "amyloid", "fibril", "prion", "misfolding",
        "inclusion body", "oligomer",
    ],
    "Folding": [
        "protein folding", "unfolding", "refolding", "folding pathway",
        "chaperone", "denaturation", "two-state", "folding kinetics",
    ],
    "Mutation": [
        "mutation", "mutant", "point mutation", "substitution",
        "single nucleotide polymorphism", "snp", "variant",
    ],
    "Structure": [
        "crystal structure", "x-ray crystallography", "cryo-em", "cryoem",
        "structure determination", "pdb", "secondary structure",
        "tertiary structure", "quaternary",
    ],
    "TP/Friction": [
        "transition path", "friction", "recrossing", "diffusion coefficient",
        "memory kernel", "generalized langevin",
    ],
    "Downhill": [
        "downhill folding", "barrierless folding", "one-state folding",
        "type 0 folding",
    ],
    "Smfs": [
        "single molecule", "single-molecule", "smfs", "optical trap",
        "magnetic tweezer",
    ],
}

# Patterns that strongly suggest a non-article document
NON_ARTICLE_PATTERNS = [
    r"\binvoice\b",
    r"\breceipt\b",
    r"\border\s+confirmation\b",
    r"\bpurchase\s+order\b",
    r"\btax\s+invoice\b",
    r"\bstatement\s+of\s+account\b",
    r"\bpayment\s+due\b",
    r"\bamount\s+due\b",
    r"\bvat\s+number\b",
    r"\bcurriculum\s+vitae\b",
    r"\bresume\b",
    r"\bcover\s+letter\b",
    r"\bboard\s+of\s+directors\b",
    r"\bslide\s+\d+\b",           # presentation slides
    r"\blecture\s+notes?\b",
    r"\bchapter\s+\d+\b",         # thesis/book chapters
    r"\bthesis\b",
    r"\bdissertation\b",
    r"\bindex\b.*\bpage\b",
]

# Known journal name fragments (not exhaustive but covers common ones)
JOURNAL_PATTERNS = [
    r"nature\s+(communications|chemistry|physics|methods|structural|chemical)",
    r"journal\s+of\s+(the\s+)?",
    r"physical\s+review\s+letters?",
    r"proc\.?\s+natl\.?\s+acad\.?\s+sci",
    r"pnas\b",
    r"j\.?\s+chem\.?\s+phys",
    r"j\.?\s+am\.?\s+chem\.?\s+soc",
    r"jacs\b",
    r"biophys\.?\s+j",
    r"proteins?:",
    r"biochemistry\b",
    r"angew\.?\s+chem",
    r"chemphyschem",
    r"plos\s+(one|comput|biol)",
    r"elife\b",
    r"science\b",
    r"cell\b",
    r"\bpeerj\b",
    r"scientific\s+reports",
    r"nucleic\s+acids",
    r"structure\b",
    r"curr\.?\s+opin",
    r"annu\.?\s+rev",
    r"chem\.?\s+rev",
    r"accounts\s+of\s+chem",
    r"biopolymers\b",
    r"protein\s+sci",
    r"j\.?\s+mol\.?\s+biol",
    r"j\.?\s+phys\.?\s+chem",
    r"langmuir\b",
    r"soft\s+matter",
    r"macromolecules\b",
    r"acs\s+(nano|cent|catal)",
    r"nanoscale\b",
]
DEFAULT_DOWNLOADS = Path.home() / "Downloads"
DEFAULT_ARTICLES = Path.home() / "Documents" / "Work" / "Articles"
DEFAULT_LOG = Path.home() / "Downloads" / "organize_papers.csv"
DEFAULT_CACHE = Path.home() / ".cache" / "pdf-organizer" / "doi-cache.json"

CSV_FIELDS = [
    "source_path",
    "dest_path",
    "original_filename",
    "new_filename",
    "author",
    "year",
    "tags",
    "is_article",
    "confidence",
    "action",
    "error",
    "doi",
    "title",
    "pages",
    "resolved_via",
    "title_match",
]

# Crossref asks for a contact address in the User-Agent for the "polite pool".
# We deliberately do NOT send one: the user's email is personal data and this
# tool has no business leaking it to a third party. The anonymous pool is fine.
USER_AGENT = "pdf-organizer/2.0"

# Article types that are commentary rather than primary research. These are
# filed only with --include-commentary.
FRONT_MATTER = re.compile(
    r"news\s*&\s*views|news and views|research briefing|correction to|"
    r"retraction note|matters arising|editorial\b|book review|milestone",
    re.IGNORECASE,
)
# Matched without word boundaries: PDF text extraction routinely glues words
# together ("SciencePERSPECTIVES", "1Matters arising").
COMMENTARY = re.compile(r"comment\b|commentary|perspectives", re.IGNORECASE)

# Nature's news/careers/comment DOI prefix — never primary research.
NEWS_FILENAME = re.compile(r"^d41586[-_]", re.IGNORECASE)

# Supplementary information, not the article itself.
SUPPLEMENTARY = re.compile(r"moesm|_esm\b|_si_\d|_si_001|supporting[-_]info", re.IGNORECASE)

CROSSREF_TYPES_OK = {
    "journal-article", "posted-content", "preprint",
    "book-chapter", "proceedings-article", "report",
}

_print_lock = Lock()
_cache_lock = Lock()


def say(msg: str) -> None:
    with _print_lock:
        print(msg, flush=True)


# ---------------------------------------------------------------------------
# Text extraction (in a subprocess, so a malformed PDF cannot hang the run)
# ---------------------------------------------------------------------------

_EXTRACT_WORKER = r"""
import sys, json, re
import pypdf

path, npages = sys.argv[1], int(sys.argv[2])
out = {"text": "", "metadata": {}, "pages": 0, "error": ""}
try:
    reader = pypdf.PdfReader(path)
    out["pages"] = len(reader.pages)
    try:
        meta = reader.metadata or {}
        out["metadata"] = {str(k).lstrip("/"): str(v) for k, v in meta.items() if v}
    except Exception:
        pass
    n = min(len(reader.pages), npages)
    text = "".join((reader.pages[i].extract_text() or "") for i in range(n))
    out["text"] = re.sub(r"\s+", " ", text)[:12000]
except Exception as exc:
    out["error"] = str(exc)[:150]
sys.stdout.write(json.dumps(out))
"""


def extract_pdf_info(pdf_path: Path, pages: int = 3, timeout: int = 45) -> dict:
    """
    Extract text + metadata from the first `pages` pages.

    pypdf has no internal timeout and will hang forever on some malformed
    files (a single scanned invoice once stalled a 608-file run at file 55),
    so it runs in a subprocess we can kill.
    """
    info = {"text": "", "metadata": {}, "pages": 0, "error": ""}
    try:
        proc = subprocess.run(
            [sys.executable, "-c", _EXTRACT_WORKER, str(pdf_path), str(pages)],
            capture_output=True, text=True, timeout=timeout,
        )
        if proc.stdout.strip():
            info.update(json.loads(proc.stdout))
    except subprocess.TimeoutExpired:
        info["error"] = f"extraction timed out after {timeout}s"
    except Exception as exc:
        info["error"] = str(exc)[:150]

    if len(info["text"].strip()) < 50:
        try:
            result = subprocess.run(
                ["gs", "-sDEVICE=txtwrite", "-dNOPAUSE", "-dBATCH", "-dQUIET",
                 "-sOutputFile=-", "-dLastPage=2", str(pdf_path)],
                capture_output=True, text=True, timeout=30,
            )
            if len(result.stdout.strip()) > len(info["text"].strip()):
                info["text"] = re.sub(r"\s+", " ", result.stdout)[:12000]
        except Exception:
            pass

    return info


# ---------------------------------------------------------------------------
# DOI extraction and resolution
# ---------------------------------------------------------------------------

def find_doi(text: str, metadata: dict) -> str:
    blob = text + "\n" + json.dumps(metadata)
    match = re.search(r"\b10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", blob)
    return match.group(0).rstrip(").,;:") if match else ""


def find_arxiv(text: str, metadata: dict, doi: str) -> str:
    m = re.search(r"(?i)arxiv[.:/ ]\s*(\d{4}\.\d{4,5})", doi or "")
    if m:
        return m.group(1)
    m = re.search(r"(?i)arxiv:\s*(\d{4}\.\d{4,5})", text + json.dumps(metadata))
    return m.group(1) if m else ""


def doi_variants(doi: str) -> list:
    """
    DOIs scraped from PDF text usually have running text glued to the end
    ("10.1101/2025.11.10.687719doi", "10.1002/advs.770861of15"). Generate
    progressively-cleaned candidates and try each.
    """
    out, seen = [], set()

    def add(cand: str) -> None:
        cand = cand.strip().rstrip(".,;:)/-")
        if cand and cand not in seen and re.match(r"^10\.\d{4,9}/\S+$", cand):
            seen.add(cand)
            out.append(cand)

    add(doi)
    add(re.sub(r"/-/DCSupplemental.*$", "", doi))
    add(re.sub(r"(?i)doi$", "", doi))
    add(re.sub(r"(?i)\d+of\d+$", "", doi))
    add(re.sub(r"(?i)of$", "", doi))
    add(re.sub(r"www\..*$", "", doi))
    add(re.sub(r"(?<=[a-z0-9])(?=[A-Z])[A-Za-z].*$", "", doi))
    m = re.match(r"^(10\.\d{4,9}/[-._;()/:A-Za-z0-9]*?\d)(?:[A-Za-z]{2,})$", doi)
    if m:
        add(m.group(1))
    add(re.sub(r"/v\d+.*$", "", doi))
    return out


def _http_json(url: str, timeout: int = 25):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.load(resp)


def _year_from(msg: dict) -> str:
    for key in ("published-print", "published-online", "issued", "created"):
        val = msg.get(key)
        if val and val.get("date-parts") and val["date-parts"][0] and val["date-parts"][0][0]:
            return str(val["date-parts"][0][0])
    return ""


def _pack_crossref(msg: dict, via: str) -> dict:
    authors = msg.get("author") or []
    family = ""
    for a in authors:
        if a.get("sequence") == "first" and a.get("family"):
            family = a["family"]
            break
    if not family and authors:
        family = authors[0].get("family") or authors[0].get("name") or ""
    return {
        "family": family,
        "year": _year_from(msg),
        "title": (msg.get("title") or [""])[0],
        "type": msg.get("type", ""),
        "container": (msg.get("container-title") or [""])[0],
        "doi": msg.get("DOI", ""),
        "via": via,
    }


def crossref_by_doi(doi: str) -> dict:
    for cand in doi_variants(doi):
        try:
            msg = _http_json(
                "https://api.crossref.org/works/" + urllib.parse.quote(cand, safe="")
            )["message"]
            if msg.get("author"):
                return _pack_crossref(msg, "crossref")
        except Exception:
            pass
        time.sleep(0.1)
    return {}


def crossref_by_title(text: str) -> dict:
    query = " ".join(text.split()[:22])
    if len(query) < 40:
        return {}
    try:
        res = _http_json(
            "https://api.crossref.org/works?rows=1&query.bibliographic="
            + urllib.parse.quote(query)
        )
        items = res["message"]["items"]
        if items and items[0].get("author"):
            return _pack_crossref(items[0], "crossref-title")
    except Exception:
        pass
    return {}


def arxiv_by_id(arxiv_id: str) -> dict:
    try:
        req = urllib.request.Request(
            "http://export.arxiv.org/api/query?id_list=" + arxiv_id,
            headers={"User-Agent": USER_AGENT},
        )
        with urllib.request.urlopen(req, timeout=25) as resp:
            root = ET.fromstring(resp.read())
        ns = {"a": "http://www.w3.org/2005/Atom"}
        entry = root.find("a:entry", ns)
        if entry is None:
            return {}
        name = entry.find("a:author/a:name", ns).text
        return {
            "family": name.split()[-1],
            "year": entry.find("a:published", ns).text[:4],
            "title": " ".join(entry.find("a:title", ns).text.split()),
            "type": "preprint",
            "container": "arXiv",
            "doi": "arXiv:" + arxiv_id,
            "via": "arxiv",
        }
    except Exception:
        return {}


class Cache:
    """Persistent DOI -> resolution cache, so re-runs cost no network."""

    def __init__(self, path: Path):
        self.path = path
        self.data = {}
        self.dirty = False
        try:
            self.data = json.loads(path.read_text())
        except Exception:
            pass

    def get(self, key: str):
        with _cache_lock:
            return self.data.get(key)

    def put(self, key: str, value: dict) -> None:
        with _cache_lock:
            self.data[key] = value
            self.dirty = True

    def save(self) -> None:
        if not self.dirty:
            return
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(json.dumps(self.data, indent=1))
        except Exception:
            pass


def resolve_identity(info: dict, cache: Cache, allow_title_search: bool) -> dict:
    """Return authoritative {family, year, title, type, ...} or {}."""
    doi = find_doi(info["text"], info["metadata"])
    arxiv_id = find_arxiv(info["text"], info["metadata"], doi)

    key = ("arxiv:" + arxiv_id) if arxiv_id else ("doi:" + doi) if doi else ""
    if key:
        hit = cache.get(key)
        if hit is not None:
            return hit

    result = {}
    if arxiv_id:
        result = arxiv_by_id(arxiv_id)
    if not result and doi:
        result = crossref_by_doi(doi)
    if not result and allow_title_search:
        result = crossref_by_title(info["text"])

    if key:
        cache.put(key, result)
    return result


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def _norm(s: str) -> str:
    """Lowercase, fold accents, reduce to [a-z0-9 ]."""
    s = s.replace("‐", "-").replace("‑", "-").replace("–", "-").replace("’", "'")
    # Fold accents first, or "Juárez" becomes "ju rez" and never matches.
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9 ]", " ", s.lower())


# Function words that carry no identifying power in a title.
TITLE_STOPWORDS = {
    "with", "from", "that", "this", "then", "than", "they", "them", "these",
    "those", "into", "onto", "over", "under", "between", "among", "during",
    "after", "before", "which", "while", "where", "when", "what", "whom",
    "whose", "been", "being", "have", "having", "using", "used", "uses",
    "also", "such", "some", "more", "most", "other", "based", "toward",
    "towards", "within", "without", "upon", "about", "across", "their",
    "there", "here", "does", "done", "both", "each", "will", "would",
}

# A title with fewer distinctive words than this cannot be verified by word
# overlap — too few words means near-certain coincidental matches. Such a
# title must instead appear contiguously in the document.
MIN_TITLE_WORDS = 5


def title_signal(title: str, text: str) -> dict:
    """
    Evidence that `text` is the document `title` belongs to.

    Returns {score, words, contiguous}: the fraction of the title's distinctive
    words present, how many such words the title has, and whether the title
    appears contiguously in the document.
    """
    out = {"score": 0.0, "words": 0, "contiguous": False}
    if not title:
        return out
    page, want = _norm(text), _norm(title)
    compact_title = want.replace(" ", "")
    if len(compact_title) > 25 and compact_title[:60] in page.replace(" ", ""):
        out["contiguous"] = True
        out["score"] = 1.0
    words = [w for w in want.split() if len(w) > 3 and w not in TITLE_STOPWORDS]
    out["words"] = len(words)
    if words and not out["contiguous"]:
        out["score"] = sum(1 for w in words if w in page) / len(words)
    return out


def author_in_text(family: str, text: str) -> bool:
    """Does the resolved first author's surname appear in the document?"""
    fam = _norm(family).strip()
    if len(fam) < 3:
        return False
    return re.search(r"(?<![a-z])" + re.escape(fam) + r"(?![a-z])", _norm(text)) is not None


def verify_identity(ident: dict, text: str, floor: float) -> tuple:
    """
    The safety gate. Returns (ok, score, reason).

    A DOI inside a PDF proves nothing about what the PDF *is* — grant
    applications, referee reports and research notes all cite DOIs and would
    otherwise resolve to, and be filed as, someone else's paper.

    Long titles are checked by word overlap. Short ones cannot be ("Provisions
    For The Use Of DNA In Accordance With Law" reduces to three usable words,
    all of which occur in unrelated legal boilerplate, and scored a perfect
    1.00 against a Zoom consent form), so they must match completely *and* be
    corroborated by the author's surname appearing in the document.
    """
    sig = title_signal(ident.get("title", ""), text)
    score = sig["score"]
    family = ident.get("family", "")
    author_ok = author_in_text(family, text)

    # The title appearing verbatim is conclusive on its own.
    if sig["contiguous"]:
        return True, score, ""

    if sig["words"] >= MIN_TITLE_WORDS:
        if score >= floor and author_ok:
            return True, score, ""
        # Overwhelming title overlap stands even if the author line did not
        # survive text extraction.
        if score >= 0.90:
            return True, score, ""
        if score < floor:
            return False, score, (
                f"resolved title absent from document ({score:.2f} < {floor:.2f})"
            )
        # Partial title overlap with no author is how documents *about* a
        # field match papers *in* it — a grant referee report scored 0.75
        # against a paper on the same topic.
        return False, score, (
            f"title only partly present ({score:.2f}) and author "
            f"'{family or '?'}' not found in document"
        )

    # Short title: demand every word plus author corroboration.
    if score >= 0.999 and author_ok:
        return True, score, ""
    return False, score, (
        f"title too short to verify ({sig['words']} words) and author "
        f"'{family or '?'}' not found in document"
    )


def commentary_kind(text: str, pages: int) -> str:
    """Return a label if this looks like commentary rather than research."""
    head = text[:220]
    m = FRONT_MATTER.search(head)
    if m and pages <= 8:
        return m.group(0).strip().lower()
    m = COMMENTARY.search(head[:170])
    if m and pages <= 6:
        return m.group(0).strip().lower()
    return ""


# ---------------------------------------------------------------------------
# Tagging
# ---------------------------------------------------------------------------

def assign_tags(text_lower: str) -> list:
    """
    Return up to 2 matching tags.

    Keywords are matched on word boundaries. Plain substring matching (the
    original behaviour) is badly wrong here: 'amp' (antimicrobial peptide)
    fires on "example" and "sample", 'rna' on "journal", 'iron' on
    "environment". Measured over 148 filed papers, 73% carried at least one
    tag produced this way.

    Custom boundaries rather than \\b because several keywords legitimately
    contain '/' or '-' (e.g. "qm/mm", "coarse-grained").
    """
    scores = {}
    for tag, keywords in TAG_KEYWORDS.items():
        hits = 0
        for kw in keywords:
            kw = kw.strip()
            if not kw:
                continue
            if re.search(r"(?<![a-z0-9])" + re.escape(kw) + r"(?![a-z0-9])", text_lower):
                hits += 1
        if hits:
            scores[tag] = hits
    return [tag for tag, _ in sorted(scores.items(), key=lambda kv: -kv[1])[:2]]
# ---------------------------------------------------------------------------
# Heuristic classification
# ---------------------------------------------------------------------------

def classify(text: str, metadata: dict) -> dict:
    """
    Returns dict with keys: is_article, author, year, tags, confidence, reason.
    Uses only heuristics — no external API.
    """
    text_lower = text.lower()
    combined = (metadata.get("Title", "") + " " + text).lower()

    # --- Non-article veto ---
    for pattern in NON_ARTICLE_PATTERNS:
        if re.search(pattern, text_lower, re.IGNORECASE):
            return {
                "is_article": False,
                "reason": f"non-article pattern: {pattern}",
                "author": "", "year": "", "tags": [], "confidence": "high",
            }

    # --- Positive signals ---
    score = 0
    signals = []

    # DOI is the strongest single signal
    doi_match = re.search(r"\b10\.\d{4,}/\S+", text)
    if doi_match:
        score += 4
        signals.append("doi")

    # Abstract section
    if re.search(r"\babstract\b", text_lower):
        score += 2
        signals.append("abstract")

    # Standard paper section headers
    section_hits = sum(
        1 for pat in [r"\bintroduction\b", r"\bmethod", r"\bresult", r"\bconclusion",
                      r"\bdiscussion\b", r"\breferences\b", r"\backnowledg"]
        if re.search(pat, text_lower)
    )
    if section_hits >= 3:
        score += 2
        signals.append(f"sections({section_hits})")
    elif section_hits >= 1:
        score += 1

    # Known journal name
    for jpat in JOURNAL_PATTERNS:
        if re.search(jpat, text_lower):
            score += 2
            signals.append("journal_name")
            break

    # Author list patterns (e.g. "Smith J," or "Smith, J." or "Smith et al.")
    if re.search(r"[A-Z][a-z]+ [A-Z]\.,?\s+[A-Z][a-z]+|et al\.", text):
        score += 1
        signals.append("author_pattern")

    # Volume/issue/page patterns
    if re.search(r"\bvol\.?\s*\d+|issue\s+\d+|\bpp\.?\s*\d+|\bpages?\s+\d+", text_lower):
        score += 1
        signals.append("vol_issue")

    # Received/accepted/published dates
    if re.search(r"received:?|accepted:?|published:?", text_lower):
        score += 1
        signals.append("pub_dates")

    # --- Decision ---
    if score >= 6:
        confidence = "high"
    elif score >= 3:
        confidence = "medium"
    elif score >= 1:
        confidence = "low"
    else:
        return {
            "is_article": False,
            "reason": "no article signals",
            "author": "", "year": "", "tags": [], "confidence": "low",
        }

    is_article = score >= 3  # require at least medium confidence

    # --- Author extraction ---
    author = extract_author(text, metadata)

    # --- Year extraction ---
    year = extract_year(text, metadata, doi_match)

    # --- Tag assignment ---
    tags = assign_tags(combined)

    return {
        "is_article": is_article,
        "author": author,
        "year": year,
        "tags": tags,
        "confidence": confidence,
        "reason": "|".join(signals),
    }


def extract_author(text: str, metadata: dict) -> str:
    """Extract first author last name, returning ASCII lowercase."""
    # 1. PDF metadata Author field
    meta_author = metadata.get("Author", "").strip()
    if meta_author:
        # Take first author if semicolon/comma separated list
        first = re.split(r"[;,]", meta_author)[0].strip()
        # Last word is typically last name
        parts = first.split()
        if parts:
            return to_ascii_lower(parts[-1])

    # 2. Look for "Firstname Lastname" before affiliations
    #    Common patterns: "John Smith1," or "Smith, John" or "J. Smith"
    patterns = [
        # "Smith J," or "Smith, J" style (last name first)
        r"^([A-Z][a-zÀ-ÿ\-']{2,}),?\s+[A-Z]",
        # "John Smith" at start of line
        r"^[A-Z][a-z]+\s+([A-Z][a-zÀ-ÿ\-']{2,})\d*[,\s]",
        # "J. Smith" style
        r"[A-Z]\.\s+([A-Z][a-zÀ-ÿ\-']{2,})\d*[,\s]",
    ]
    for line in text.splitlines()[:60]:
        line = line.strip()
        if not line or len(line) > 120:
            continue
        for pat in patterns:
            m = re.match(pat, line)
            if m:
                candidate = m.group(1)
                # Reject common non-name words
                if candidate.lower() not in {
                    "abstract", "introduction", "methods", "results",
                    "discussion", "journal", "received", "accepted",
                    "published", "email", "copyright", "figure",
                }:
                    return to_ascii_lower(candidate)

    return ""


def extract_year(text: str, metadata: dict, doi_match) -> str:
    """Extract 4-digit publication year."""
    # 1. From DOI — many DOIs embed year: 10.1021/jacs.2023.xxxxx or similar
    if doi_match:
        m = re.search(r"\b(19|20)\d{2}\b", doi_match.group(0))
        if m:
            return m.group(0)

    # 2. PDF metadata CreationDate or ModDate
    for key in ("CreationDate", "ModDate"):
        val = metadata.get(key, "")
        m = re.search(r"(19|20)\d{2}", val)
        if m:
            yr = int(m.group(0))
            if 1990 <= yr <= 2030:
                return m.group(0)

    # 3. Received/accepted/published dates in text
    m = re.search(
        r"(?:received|accepted|published)[^\n]{0,30}(20\d{2}|19\d{2})",
        text, re.IGNORECASE
    )
    if m:
        return m.group(1)

    # 4. Copyright year
    m = re.search(r"©\s*(20\d{2}|19\d{2})", text)
    if m:
        return m.group(1)

    # 5. Any year in range appearing in first 1000 chars
    years = re.findall(r"\b(20[0-2]\d|19[89]\d)\b", text[:1000])
    if years:
        # Prefer most common year
        from collections import Counter
        return Counter(years).most_common(1)[0][0]

    return ""

# ---------------------------------------------------------------------------
# Filename helpers
# ---------------------------------------------------------------------------

def to_ascii_lower(name: str) -> str:
    """Normalize unicode name to ASCII lowercase alphanumeric."""
    normalized = unicodedata.normalize("NFKD", name)
    ascii_str = "".join(c for c in normalized if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9]", "", ascii_str.lower())

# ---------------------------------------------------------------------------
# Finder tag application
# ---------------------------------------------------------------------------

def apply_tag(filepath: Path, tag_names: list) -> None:
    """Apply Finder tags to a file using xattr."""
    if not tag_names:
        return
    tag_entries = [f"{name}\n{KNOWN_TAGS.get(name, '0')}" for name in tag_names]
    plist_data = plistlib.dumps(tag_entries, fmt=plistlib.FMT_BINARY)
    subprocess.run(
        ["xattr", "-wx", "com.apple.metadata:_kMDItemUserTags",
         plist_data.hex(), str(filepath)],
        check=True, capture_output=True,
    )
# ---------------------------------------------------------------------------
# Deduplication against the destination library
# ---------------------------------------------------------------------------

def normalise_doi(doi: str) -> str:
    return re.sub(r"[^a-z0-9./]", "", (doi or "").lower())


def library_doi_index(dest_dir: Path, stems: set, cache: Cache, jobs: int,
                      timeout: int) -> dict:
    """
    Map normalised DOI -> existing filename, for library files that could
    collide with an incoming paper.

    Only files whose `author-year` stem matches an incoming candidate are
    opened — scanning a 6000-file library in full would be pointless.
    """
    targets = []
    for path in sorted(dest_dir.glob("*.pdf")):
        m = re.match(r"([a-z\-']+)-(\d{4})[a-z]?\.pdf$", path.name)
        if m and (m.group(1), m.group(2)) in stems:
            targets.append(path)
    if not targets:
        return {}

    say(f"Checking {len(targets)} existing library file(s) for duplicate DOIs...")
    index = {}

    def scan(path: Path):
        try:
            key = f"lib:{path.name}:{path.stat().st_size}"
        except OSError:
            return None
        hit = cache.get(key)
        if hit is None:
            info = extract_pdf_info(path, pages=2, timeout=timeout)
            hit = {"doi": normalise_doi(find_doi(info["text"], info["metadata"]))}
            cache.put(key, hit)
        return path.name, hit.get("doi", "")

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for res in pool.map(scan, targets):
            if not res:
                continue
            name, doi = res
            # Short DOIs are truncation artefacts ("10.1073/pnas") and would
            # match every paper from that publisher. Require a full one.
            if len(doi) >= 16:
                index.setdefault(doi, name)
    return index


def resolve_filename(author: str, year: str, taken: set) -> str:
    base = f"{author}-{year}"
    if f"{base}.pdf" not in taken:
        return f"{base}.pdf"
    for suffix in "abcdefghijklmnopqrstuvwxyz":
        if f"{base}{suffix}.pdf" not in taken:
            return f"{base}{suffix}.pdf"
    n = 2
    while f"{base}-{n}.pdf" in taken:
        n += 1
    return f"{base}-{n}.pdf"


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log_action(csv_writer, record: dict) -> None:
    row = {field: record.get(field, "") for field in CSV_FIELDS}
    if isinstance(row["tags"], list):
        row["tags"] = "|".join(row["tags"])
    csv_writer.writerow(row)


# ---------------------------------------------------------------------------
# Per-file analysis (no filesystem changes)
# ---------------------------------------------------------------------------

def analyse(pdf_path: Path, args, cache: Cache) -> dict:
    rec = {
        "source_path": str(pdf_path),
        "original_filename": pdf_path.name,
        "dest_path": "", "new_filename": "",
        "author": "", "year": "", "tags": [],
        "is_article": False, "confidence": "", "action": "", "error": "",
        "doi": "", "title": "", "pages": 0, "resolved_via": "", "title_match": "",
    }

    info = extract_pdf_info(pdf_path, pages=3, timeout=args.timeout)
    text, meta, rec["pages"] = info["text"], info["metadata"], info["pages"]
    if info["error"]:
        rec["error"] = info["error"]
    lower = text.lower()

    # --- cheap structural vetoes -------------------------------------------
    if SUPPLEMENTARY.search(pdf_path.name) and not args.include_supplementary:
        rec.update(action="skipped", error="supplementary material")
        return rec
    if NEWS_FILENAME.match(pdf_path.name) and not args.include_commentary:
        rec.update(action="skipped", error="Nature news/careers (d41586)")
        return rec
    for pattern in NON_ARTICLE_PATTERNS:
        if re.search(pattern, lower, re.IGNORECASE):
            rec.update(action="skipped", error=f"non-article pattern: {pattern}")
            return rec

    # --- identity ----------------------------------------------------------
    ident = {}
    if not args.offline:
        ident = resolve_identity(info, cache, allow_title_search=True)

    if ident and ident.get("family") and ident.get("year"):
        rec["doi"] = ident.get("doi", "")
        rec["title"] = ident.get("title", "")
        rec["resolved_via"] = ident.get("via", "")
        scanned = len(text.strip()) < 200

        # A fuzzy bibliographic search is far weaker evidence than a DOI
        # lookup, so it has to clear a higher bar. Measured on one corpus:
        # correct matches score 1.00, wrong ones 0.56-0.67.
        searched = ident.get("via") == "crossref-title"
        floor = args.min_title_match_search if searched else args.min_title_match

        if scanned and not searched:
            # No usable text (a scan): the title gate cannot run. The DOI came
            # from the file itself, so accept it but say so in the log.
            rec["title_match"] = "n/a (no text)"
            rec["confidence"] = "medium"
        else:
            ok, score, why = verify_identity(ident, text, floor)
            rec["title_match"] = f"{score:.2f}"
            if not ok:
                rec.update(action="skipped", error=why)
                return rec
            rec["confidence"] = "high"

        if ident.get("type") and ident["type"] not in CROSSREF_TYPES_OK \
                and not args.include_commentary:
            rec.update(action="skipped", error=f"not a research item ({ident['type']})")
            return rec

        kind = commentary_kind(text, rec["pages"] or 99)
        if kind and not args.include_commentary:
            rec.update(action="skipped",
                       error=f"{kind} piece, not research ({rec['pages']}p)")
            return rec

        rec["author"] = to_ascii_lower(ident["family"])
        rec["year"] = ident["year"]
        rec["is_article"] = True
        rec["tags"] = assign_tags(lower) if args.tags else []
        if not rec["author"] or not re.match(r"^\d{4}$", rec["year"]):
            rec.update(action="error",
                       error=f"bad metadata: author='{rec['author']}' year='{rec['year']}'")
        return rec

    # --- fallback: text heuristics (offline, or no DOI found) --------------
    # Off by default. This path cannot be verified — there is no resolved
    # title to check against the document — and it is the code that produced
    # `the-2013.pdf`, `administrator-2026.pdf` and `pythondocx-2026.pdf`.
    # Measured 33% correct over 608 files. A file with no resolvable DOI is
    # left in place rather than filed under a guess.
    if not args.allow_heuristic:
        rec.update(action="skipped",
                   error="no DOI resolved (use --allow-heuristic to guess)")
        return rec

    if not text.strip():
        rec.update(action="skipped", error="no text extracted and no DOI")
        return rec

    result = classify(text, meta)
    rec["confidence"] = result["confidence"]
    rec["resolved_via"] = "heuristic"
    rec["tags"] = result.get("tags", []) if args.tags else []
    order = {"high": 3, "medium": 2, "low": 1}
    if not result["is_article"] or order.get(result["confidence"], 0) < order.get(args.min_confidence, 2):
        rec.update(action="skipped",
                   error=result.get("reason", "not_article")[:70])
        return rec

    rec["author"] = to_ascii_lower(result.get("author", ""))
    rec["year"] = result.get("year", "")
    rec["is_article"] = True
    if not rec["author"] or not re.match(r"^\d{4}$", rec["year"]):
        rec.update(action="error",
                   error=f"bad metadata: author='{rec['author']}' year='{rec['year']}'")
    return rec


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="File journal articles into a firstauthor-year.pdf library.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Preview only, no files moved")
    parser.add_argument("--limit", type=int, default=0, help="Process at most N files (0=all)")
    parser.add_argument("--downloads", type=Path, default=DEFAULT_DOWNLOADS)
    parser.add_argument("--articles", type=Path, default=DEFAULT_ARTICLES)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--jobs", type=int, default=8, help="Parallel workers (default 8)")
    parser.add_argument("--timeout", type=int, default=45,
                        help="Per-file text extraction timeout in seconds")
    parser.add_argument("--offline", action="store_true",
                        help="Skip Crossref/arXiv; use text heuristics (much less accurate)")
    parser.add_argument("--include-commentary", action="store_true",
                        help="Also file News & Views, Comments, Perspectives, Nature news")
    parser.add_argument("--include-supplementary", action="store_true",
                        help="Also file supplementary-information PDFs")
    parser.add_argument("--no-dedup", action="store_true",
                        help="Do not check whether a paper is already in the library")
    parser.add_argument("--tags", action="store_true",
                        help="Apply Finder tags from the keyword vocabulary. Off by "
                             "default: keyword matching is too coarse to be trusted, "
                             "and files are better left unlabelled than mislabelled")
    parser.add_argument("--min-title-match", type=float, default=0.70,
                        help="For DOI-resolved papers: fraction of resolved-title words "
                             "that must appear in the document (default 0.70). This is "
                             "the safety gate.")
    parser.add_argument("--min-title-match-search", type=float, default=0.85,
                        help="Same gate for results from fuzzy title search, which is "
                             "weaker evidence than a DOI (default 0.85)")
    parser.add_argument("--allow-heuristic", action="store_true",
                        help="Fall back to guessing author/year from text when no DOI "
                             "resolves. Unverifiable and ~33%% accurate; off by default")
    parser.add_argument(
        "--min-confidence", choices=["high", "medium", "low"], default="medium",
        help="Threshold for the heuristic fallback path only",
    )
    args = parser.parse_args()

    # Offline means there is nothing to resolve against, so the heuristic
    # path is the only one available.
    if args.offline:
        args.allow_heuristic = True

    if not args.downloads.is_dir():
        sys.exit(f"Downloads directory not found: {args.downloads}")
    if not args.dry_run:
        args.articles.mkdir(parents=True, exist_ok=True)

    pdfs = sorted({p for p in args.downloads.glob("*.[Pp][Dd][Ff]")})
    if args.limit:
        pdfs = pdfs[: args.limit]
    if not pdfs:
        sys.exit(f"No PDFs found in {args.downloads}")

    print(f"Found {len(pdfs)} PDF(s) in {args.downloads}")
    if args.dry_run:
        print("DRY RUN — no files will be moved")
    if args.offline:
        print("OFFLINE — falling back to text heuristics; expect poor naming accuracy")
    print(f"Log: {args.log}\n")

    cache = Cache(args.cache)

    # --- phase 1+2: analyse everything in parallel -------------------------
    done = [0]

    def work(path: Path) -> dict:
        rec = analyse(path, args, cache)
        done[0] += 1
        if done[0] % 25 == 0:
            say(f"  ...analysed {done[0]}/{len(pdfs)}")
        return rec

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        records = list(pool.map(work, pdfs))
    cache.save()

    keep = [r for r in records if r["is_article"] and not r["action"]]

    # --- phase 3: deduplicate ---------------------------------------------
    if not args.no_dedup:
        # within this batch: same DOI -> keep the longest document (a 2-page
        # Research Briefing is not a substitute for the 33-page article)
        by_doi = {}
        for r in keep:
            d = normalise_doi(r["doi"])
            if len(d) >= 14:
                by_doi.setdefault(d, []).append(r)
        for group in by_doi.values():
            if len(group) > 1:
                group.sort(key=lambda r: (-(r["pages"] or 0), len(r["original_filename"])))
                for dup in group[1:]:
                    dup.update(action="skipped",
                               error=f"duplicate of {group[0]['original_filename'][:44]}")
        keep = [r for r in keep if not r["action"]]

        stems = {(r["author"], r["year"]) for r in keep}
        index = library_doi_index(args.articles, stems, cache, max(1, args.jobs), args.timeout)
        cache.save()
        for r in keep:
            d = normalise_doi(r["doi"])
            if len(d) >= 16 and d in index:
                r.update(action="skipped", error=f"already in library as {index[d]}")
        keep = [r for r in keep if not r["action"]]

    # --- phase 4: assign names and act ------------------------------------
    taken = {p.name for p in args.articles.glob("*.pdf")} if args.articles.is_dir() else set()
    for r in sorted(keep, key=lambda r: (r["author"], r["year"], r["original_filename"])):
        name = resolve_filename(r["author"], r["year"], taken)
        taken.add(name)
        r["new_filename"] = name
        r["dest_path"] = str(args.articles / name)

    counts = {"moved": 0, "dry-run": 0, "skipped": 0, "error": 0}
    with open(args.log, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()

        for r in records:
            if r["action"] in ("skipped", "error"):
                counts[r["action"]] += 1
                if args.verbose or r["action"] == "error":
                    print(f"  [{r['action']:7}] {r['original_filename'][:52]}  ({r['error'][:56]})")
                log_action(writer, r)

        for r in sorted(keep, key=lambda r: r["new_filename"]):
            tag_str = "|".join(r["tags"]) if r["tags"] else "none"
            label = "dry-run" if args.dry_run else "move"
            print(f"  [{label}] {r['original_filename'][:44]} → {r['new_filename']}  tags={tag_str}")
            if args.verbose:
                print(f"           via={r['resolved_via']} title_match={r['title_match']} "
                      f"pages={r['pages']} doi={r['doi'][:40]}")
            if args.dry_run:
                r["action"] = "dry-run"
                counts["dry-run"] += 1
            else:
                try:
                    shutil.move(r["source_path"], r["dest_path"])
                    if r["tags"]:
                        apply_tag(Path(r["dest_path"]), r["tags"])
                    r["action"] = "moved"
                    counts["moved"] += 1
                except Exception as exc:
                    r["action"] = "error"
                    r["error"] = str(exc)[:120]
                    counts["error"] += 1
                    print(f"  [error] {r['original_filename']}: {exc}")
            log_action(writer, r)

    print("\n" + "=" * 50)
    verb = "would move" if args.dry_run else "moved"
    print(f"{counts['dry-run'] or counts['moved']} {verb}, "
          f"{counts['skipped']} skipped, {counts['error']} error")
    print(f"Log written to: {args.log}")


if __name__ == "__main__":
    main()
