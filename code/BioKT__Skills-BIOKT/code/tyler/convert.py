#!/usr/bin/env python3
"""
PDF-to-Wiki converter for academic literature.

Converts a directory of academic PDFs into a structured markdown wiki with:
- YAML frontmatter (title, authors, year, journal, DOI, abstract, keywords)
- Cleaned body text (stripped headers/footers/page numbers/publisher boilerplate)
- Auto-generated index.md from extracted metadata
- Incremental mode (skips already-converted files)
- Optional BibTeX matching
- Optional recursive directory scanning

Adapted from the tyler skill by Johan Fourie (github.com/johanfourieza/econtools).
"""

import sys
import os
import re
import json
import hashlib
import argparse
from pathlib import Path


# ---------------------------------------------------------------------------
# Filename helpers
# ---------------------------------------------------------------------------

def sanitise_filename(name):
    """Convert a PDF filename to a safe markdown filename."""
    name = os.path.splitext(name)[0]
    name = re.sub(r'[^\w\s-]', '', name)
    name = re.sub(r'[\s]+', '_', name)
    name = name.strip('_')
    return name + '.md'


def file_hash(path):
    """Return SHA-256 hex digest of a file (for incremental mode)."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Metadata extraction heuristics
# ---------------------------------------------------------------------------

def clean_field(text):
    """Remove markdown formatting artifacts from extracted metadata fields."""
    text = re.sub(r'\*\*?', '', text)
    text = re.sub(r'_([^_]+)_', r'\1', text)
    text = re.sub(r'\[(\d+)\]', '', text)
    text = re.sub(r'==>.*?<==', '', text)
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def is_junk_line(line):
    """Return True if a line is publisher boilerplate."""
    line_lower = line.strip().lower()
    junk_patterns = [
        r'==>.*?<==',
        r'downloaded from',
        r'this article is licensed under',
        r'creative commons',
        r'all rights reserved',
        r'published by elsevier',
        r'elsevier b\.v\.',
        r'springer nature',
        r'© 20\d\d',
        r'copyright 20\d\d',
        r'biorxiv preprint',
        r'medrxiv preprint',
        r'posted on biorxiv',
        r'this preprint reports',
        r'pmc funded article',
        r'author manuscript',
        r'nih public access',
        r'europe pmc',
        r'accepted manuscript',
    ]
    for pat in junk_patterns:
        if re.search(pat, line_lower):
            return True
    return False


def extract_title(md_text):
    """Extract the paper title from markdown text."""
    for m in re.finditer(r'^#{1,3}\s+(.+)$', md_text, re.MULTILINE):
        title = m.group(1).strip()
        if len(title) > 10 and len(title) < 300 and '==>' not in title:
            return clean_field(title)

    for m in re.finditer(r'\*\*(.{10,200}?)\*\*', md_text):
        candidate = m.group(1).strip()
        if '==>' not in candidate and not is_junk_line(candidate):
            return clean_field(candidate)

    for line in md_text.split('\n'):
        line = line.strip()
        if len(line) > 15 and not is_junk_line(line) and '==>' not in line:
            return clean_field(line[:200])

    return "Unknown Title"


def extract_authors(md_text, title):
    """Extract authors — typically the lines just after the title."""
    lines = md_text.split('\n')
    title_idx = None

    title_clean = re.sub(r'[#*_\s]+', ' ', title).strip().lower()
    for i, line in enumerate(lines):
        line_clean = re.sub(r'[#*_\s]+', ' ', line).strip().lower()
        if title_clean and title_clean in line_clean:
            title_idx = i
            break

    if title_idx is None:
        return "Unknown"

    candidates = []
    for line in lines[title_idx + 1: title_idx + 8]:
        line = line.strip().strip('*').strip('_').strip()
        if not line:
            continue
        if is_junk_line(line) or '==>' in line:
            continue
        if re.match(r'^(abstract|introduction|#{1,3}\s)', line, re.IGNORECASE):
            break
        if re.match(r'^\d{4}', line):
            break
        if '@' in line and ',' not in line:
            continue
        if len(line) > 200:
            continue
        if len(line) > 5:
            candidates.append(clean_field(line))
        if len(candidates) >= 3:
            break

    if candidates:
        return '; '.join(candidates)
    return "Unknown"


def extract_year(md_text, filename):
    """Extract publication year from text or filename."""
    m = re.search(r'(19|20)\d{2}', filename)
    if m:
        return m.group(0)

    header = md_text[:2000]
    years = re.findall(r'((?:19|20)\d{2})', header)
    if years:
        valid = [y for y in years if 1900 <= int(y) <= 2030]
        if valid:
            return max(valid)

    return "Unknown"


def extract_abstract(md_text):
    """Extract the abstract section from an academic paper."""
    patterns = [
        r'(?:^|\n)\s*#{0,3}\s*\**\s*Abstract\s*\**\s*\n+(.*?)(?=\n\s*#{1,3}\s|\n\s*\**\s*(?:Introduction|1[\.\s]|Keywords)\s*\**)',
        r'(?:^|\n)\s*\**Abstract[\.:]\**\s*(.*?)(?=\n\s*\**\s*(?:Introduction|1[\.\s]|Keywords))',
        r'(?:^|\n)\s*\**Abstract\**\s*[-—:.]?\s*\n*(.*?)(?=\n\n\n|\n\s*\**(?:Introduction|1[\.\s]))',
    ]

    for pattern in patterns:
        m = re.search(pattern, md_text, re.IGNORECASE | re.DOTALL)
        if m:
            abstract = m.group(1).strip()
            abstract = re.sub(r'\s+', ' ', abstract)
            if len(abstract) > 50:
                return abstract[:3000]

    return ""


def extract_keywords(md_text):
    """Extract keywords if present."""
    m = re.search(
        r'(?:^|\n)\s*\**Keywords?\**[\s:]+(.+?)(?:\n\n|\n\s*\**(?:Introduction|1[\.\s]))',
        md_text[:5000], re.IGNORECASE | re.DOTALL
    )
    if m:
        kw = m.group(1).strip()
        kw = re.sub(r'\s+', ' ', kw)
        return kw[:500]
    return ""


def extract_doi(md_text):
    """Extract DOI from the first ~3000 characters of the paper."""
    header = md_text[:3000]
    patterns = [
        r'https?://doi\.org/(10\.\S+)',
        r'https?://dx\.doi\.org/(10\.\S+)',
        r'doi[:\s]+(10\.\d{4,}/\S+)',
        r'DOI[:\s]+(10\.\d{4,}/\S+)',
    ]
    for pat in patterns:
        m = re.search(pat, header, re.IGNORECASE)
        if m:
            doi = m.group(1).rstrip('.,;)')
            return doi
    return ""


def extract_journal(md_text):
    """Extract journal name from the first ~2000 characters of the paper."""
    header = md_text[:2000]

    known_journals = [
        r'Journal of Chemical Physics',
        r'Journal of Chemical Theory and Computation',
        r'Journal of Physical Chemistry [A-Z]?',
        r'Biophysical Journal',
        r'Biophysical Letters',
        r'Journal of the American Chemical Society',
        r'Proceedings of the National Academy of Sciences',
        r'PNAS',
        r'Nature Methods',
        r'Nature Communications',
        r'Nature Chemistry',
        r'Nature',
        r'Science',
        r'eLife',
        r'PLOS Computational Biology',
        r'PLOS ONE',
        r'Proteins[:\s]',
        r'Protein Science',
        r'Structure',
        r'Nucleic Acids Research',
        r'Physical Chemistry Chemical Physics',
        r'Physical Review Letters',
        r'Physical Review E',
        r'Journal of Molecular Biology',
        r'Journal of Biological Chemistry',
        r'Biochemistry',
        r'ChemPhysChem',
        r'Angewandte Chemie',
        r'Chemical Science',
        r'Journal of Computational Chemistry',
        r'Wiley Interdisciplinary Reviews',
        r'Annual Review of Biophysics',
        r'Current Opinion in Structural Biology',
    ]

    for journal in known_journals:
        m = re.search(journal, header, re.IGNORECASE)
        if m:
            return m.group(0).strip().rstrip(':').strip()

    m = re.search(r'(?:Published in|Journal)[:\s]+([A-Z][^\n]{5,60})', header, re.IGNORECASE)
    if m:
        return clean_field(m.group(1))

    return ""


# ---------------------------------------------------------------------------
# Text cleaning
# ---------------------------------------------------------------------------

def clean_markdown(md_text):
    """Clean pymupdf4llm output for token efficiency."""
    lines = md_text.split('\n')

    lines = [l for l in lines if not is_junk_line(l)]
    lines = [l for l in lines if '==>' not in l]

    if len(lines) > 50:
        freq = {}
        for line in lines:
            stripped = line.strip()
            if 3 < len(stripped) < 80:
                freq[stripped] = freq.get(stripped, 0) + 1
        repeated = {k for k, v in freq.items() if v >= 3}
        lines = [l for l in lines if l.strip() not in repeated]

    lines = [l for l in lines if not re.match(r'^\s*-?\s*\d{1,4}\s*-?\s*$', l)]

    cleaned = []
    blank_count = 0
    for line in lines:
        if line.strip() == '':
            blank_count += 1
            if blank_count <= 2:
                cleaned.append(line)
        else:
            blank_count = 0
            cleaned.append(line)

    return '\n'.join(cleaned)


def trim_references(md_text):
    """Remove the references/bibliography section to save tokens."""
    patterns = [
        r'\n\s*#{0,3}\s*\**\s*References\s*\**\s*\n',
        r'\n\s*#{0,3}\s*\**\s*Bibliography\s*\**\s*\n',
        r'\n\s*#{0,3}\s*\**\s*Works Cited\s*\**\s*\n',
        r'\n\s*\**\s*REFERENCES\s*\**\s*\n',
    ]

    for pattern in patterns:
        m = re.search(pattern, md_text)
        if m:
            body = md_text[:m.start()]
            refs = md_text[m.start():]
            return body.rstrip(), refs.strip()

    return md_text, ""


# ---------------------------------------------------------------------------
# BibTeX matching
# ---------------------------------------------------------------------------

def load_bibtex(bib_path):
    """Parse a .bib file and return a list of entries with citekey, title, authors, year."""
    entries = []
    with open(bib_path, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    for m in re.finditer(r'@\w+\{(\w+)\s*,\s*(.*?)\n\}', content, re.DOTALL):
        citekey = m.group(1)
        body = m.group(2)

        title_m = re.search(r'title\s*=\s*[\{"](.+?)[\}"]', body, re.IGNORECASE)
        author_m = re.search(r'author\s*=\s*[\{"](.+?)[\}"]', body, re.IGNORECASE)
        year_m = re.search(r'year\s*=\s*[\{"]?(\d{4})[\}"]?', body, re.IGNORECASE)

        entries.append({
            'citekey': citekey,
            'title': title_m.group(1) if title_m else '',
            'author': author_m.group(1) if author_m else '',
            'year': year_m.group(1) if year_m else '',
        })

    return entries


def match_bibtex(title, authors, year, bib_entries):
    """Fuzzy-match a paper to a BibTeX entry. Returns citekey or empty string."""
    if not bib_entries:
        return ""

    title_lower = re.sub(r'[^\w\s]', '', title.lower())
    title_words = set(title_lower.split())

    best_score = 0
    best_key = ""

    for entry in bib_entries:
        score = 0
        entry_title = re.sub(r'[^\w\s]', '', entry['title'].lower())
        entry_words = set(entry_title.split())

        if title_words and entry_words:
            overlap = len(title_words & entry_words) / max(len(title_words), len(entry_words))
            score += overlap * 10

        if year != "Unknown" and entry['year'] == year:
            score += 2

        if authors != "Unknown" and entry['author']:
            author_words = set(re.sub(r'[^\w\s]', '', authors.lower()).split())
            bib_author_words = set(re.sub(r'[^\w\s]', '', entry['author'].lower()).split())
            if author_words & bib_author_words:
                score += 3

        if score > best_score:
            best_score = score
            best_key = entry['citekey']

    return best_key if best_score >= 6 else ""


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def convert_one_pdf(pdf_path, output_path, bib_entries=None, keep_references=False):
    """Convert a single PDF to cleaned markdown with YAML frontmatter.

    Returns a metadata dict for the index, or None on failure.
    """
    import pymupdf4llm

    md_text = pymupdf4llm.to_markdown(pdf_path)
    original_filename = os.path.basename(pdf_path)

    title = extract_title(md_text)
    authors = extract_authors(md_text, title)
    year = extract_year(md_text, original_filename)
    abstract = extract_abstract(md_text)
    keywords = extract_keywords(md_text)
    doi = extract_doi(md_text)
    journal = extract_journal(md_text)

    citekey = ""
    if bib_entries:
        citekey = match_bibtex(title, authors, year, bib_entries)

    cleaned = clean_markdown(md_text)

    body, refs = trim_references(cleaned)
    if keep_references:
        body = cleaned

    fm_lines = ['---']
    fm_lines.append(f'title: "{title}"')
    fm_lines.append(f'authors: "{authors}"')
    fm_lines.append(f'year: "{year}"')
    if journal:
        fm_lines.append(f'journal: "{journal}"')
    if doi:
        fm_lines.append(f'doi: "{doi}"')
    if citekey:
        fm_lines.append(f'citekey: "{citekey}"')
    if keywords:
        fm_lines.append(f'keywords: "{keywords}"')
    fm_lines.append(f'source_pdf: "{original_filename}"')
    fm_lines.append('---')
    fm_lines.append('')

    if abstract:
        fm_lines.append('## Abstract')
        fm_lines.append('')
        fm_lines.append(abstract)
        fm_lines.append('')
        fm_lines.append('---')
        fm_lines.append('')

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(fm_lines) + '\n' + body)

    return {
        'title': title,
        'authors': authors,
        'year': year,
        'journal': journal,
        'doi': doi,
        'abstract': abstract[:800],
        'citekey': citekey,
        'keywords': keywords,
        'source_pdf': original_filename,
        'md_filename': os.path.basename(output_path),
        'body_tokens_approx': len(body.split()),
    }


# ---------------------------------------------------------------------------
# Index generation
# ---------------------------------------------------------------------------

def build_index(metadata_list, output_path):
    """Build index.md from extracted metadata."""
    n = len(metadata_list)
    lines = [
        '# Literature Wiki — Index',
        '',
        f'This index covers **{n} papers**. Each entry has structured metadata and the abstract.',
        'To read a full paper, open the linked `.md` file in `papers/`.',
        '',
        '**How to use this wiki:**',
        '- Read this index to understand what each paper covers',
        '- Use Grep to search across all papers for specific terms, methods, or systems',
        '- Read individual `papers/*.md` files only when you need full detail',
        '',
        '---',
        '',
    ]

    for meta in sorted(metadata_list, key=lambda m: (m['year'], m['title'])):
        lines.append(f'## {meta["title"]}')
        lines.append('')
        lines.append(f'**Authors:** {meta["authors"]}  ')
        lines.append(f'**Year:** {meta["year"]}  ')
        if meta.get('journal'):
            lines.append(f'**Journal:** {meta["journal"]}  ')
        if meta.get('doi'):
            lines.append(f'**DOI:** {meta["doi"]}  ')
        if meta.get('citekey'):
            lines.append(f'**Citekey:** `{meta["citekey"]}`  ')
        if meta.get('keywords'):
            lines.append(f'**Keywords:** {meta["keywords"]}  ')
        lines.append(f'**File:** `papers/{meta["md_filename"]}`  ')
        lines.append('')
        if meta['abstract']:
            lines.append(f'> {meta["abstract"]}')
            lines.append('')
        lines.append('---')
        lines.append('')

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))

    return output_path


# ---------------------------------------------------------------------------
# State tracking (incremental mode)
# ---------------------------------------------------------------------------

def load_state(state_path):
    if os.path.exists(state_path):
        with open(state_path, 'r') as f:
            return json.load(f)
    return {}


def save_state(state_path, state):
    with open(state_path, 'w') as f:
        json.dump(state, f, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert academic PDFs to a structured markdown wiki.'
    )
    parser.add_argument('pdf_dir', help='Directory containing PDF files')
    parser.add_argument('wiki_dir', help='Output wiki directory')
    parser.add_argument('--recursive', '-r', action='store_true',
                        help='Scan pdf_dir recursively for PDFs')
    parser.add_argument('--bib', type=str, default=None,
                        help='Path to a .bib file for citation key matching')
    parser.add_argument('--keep-references', action='store_true',
                        help='Keep the references section (default: trim it)')
    parser.add_argument('--force', action='store_true',
                        help='Re-convert all files, ignoring incremental state')

    args = parser.parse_args()

    pdf_dir = args.pdf_dir
    wiki_dir = args.wiki_dir

    if not os.path.isdir(pdf_dir):
        print(f"ERROR: PDF directory not found: {pdf_dir}")
        sys.exit(1)

    papers_dir = os.path.join(wiki_dir, 'papers')
    os.makedirs(papers_dir, exist_ok=True)

    if args.recursive:
        pdf_files = []
        for root, dirs, files in os.walk(pdf_dir):
            for f in files:
                if f.lower().endswith('.pdf'):
                    pdf_files.append(os.path.join(root, f))
    else:
        pdf_files = [os.path.join(pdf_dir, f) for f in os.listdir(pdf_dir)
                     if f.lower().endswith('.pdf')]

    if not pdf_files:
        print(f"No PDF files found in: {pdf_dir}")
        sys.exit(0)

    pdf_files.sort()
    print(f"Found {len(pdf_files)} PDF files.\n")

    state_path = os.path.join(wiki_dir, '.wiki_state.json')
    state = load_state(state_path) if not args.force else {}

    bib_entries = None
    if args.bib:
        if os.path.isfile(args.bib):
            bib_entries = load_bibtex(args.bib)
            print(f"Loaded {len(bib_entries)} BibTeX entries from {args.bib}\n")
        else:
            print(f"WARNING: .bib file not found: {args.bib}\n")

    succeeded = []
    skipped = []
    failed = []
    metadata_list = []

    for i, pdf_path in enumerate(pdf_files, 1):
        pdf_filename = os.path.basename(pdf_path)
        md_filename = sanitise_filename(pdf_filename)
        output_path = os.path.join(papers_dir, md_filename)

        current_hash = file_hash(pdf_path)
        if pdf_path in state and state[pdf_path].get('hash') == current_hash and not args.force:
            cached_meta = state[pdf_path].get('metadata')
            if cached_meta:
                metadata_list.append(cached_meta)
                skipped.append(pdf_filename)
                print(f"  [{i}/{len(pdf_files)}] SKIP (unchanged): {pdf_filename}")
                continue

        try:
            meta = convert_one_pdf(pdf_path, output_path, bib_entries, args.keep_references)

            file_size = os.path.getsize(output_path)
            if file_size < 500:
                print(f"  [{i}/{len(pdf_files)}] WARNING (possible scan, little text): {pdf_filename}")
            else:
                print(f"  [{i}/{len(pdf_files)}] OK: {pdf_filename} -> {md_filename}")

            succeeded.append(pdf_filename)
            metadata_list.append(meta)
            state[pdf_path] = {'hash': current_hash, 'metadata': meta}

        except Exception as e:
            print(f"  [{i}/{len(pdf_files)}] FAILED: {pdf_filename} | Error: {e}")
            failed.append((pdf_filename, str(e)))

    save_state(state_path, state)

    index_path = os.path.join(wiki_dir, 'index.md')
    build_index(metadata_list, index_path)

    print(f"\n{'='*50}")
    print(f"  Conversion complete")
    print(f"{'='*50}")
    print(f"  New/updated: {len(succeeded)}")
    print(f"  Skipped:     {len(skipped)}")
    print(f"  Failed:      {len(failed)}")
    print(f"  Total index: {len(metadata_list)} papers")
    print(f"  Wiki:        {wiki_dir}")
    print(f"  Index:       {index_path}")

    if failed:
        print(f"\n  Failed files:")
        for name, err in failed:
            print(f"    - {name}: {err}")

    total_papers = len(metadata_list)
    pdf_token_est = total_papers * 10000
    index_tokens = total_papers * 400
    print(f"\n  Token estimate:")
    print(f"    Reading all PDFs directly:  ~{pdf_token_est:,} tokens")
    print(f"    Reading index.md only:      ~{index_tokens:,} tokens")
    print(f"    Savings:                    ~{pdf_token_est - index_tokens:,} tokens ({((pdf_token_est - index_tokens)/max(pdf_token_est,1))*100:.0f}%)")


if __name__ == '__main__':
    main()
