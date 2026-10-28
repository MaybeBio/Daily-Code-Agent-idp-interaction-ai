#!/usr/bin/env python3
"""
Concatenates multiple InterProScan XML files (protein-matches format)
into a single merged XML file.

Usage:
    python concat_xml.py                          # merges all *.xml in cwd -> merged_output.xml
    python concat_xml.py -o my_output.xml         # custom output filename
    python concat_xml.py -p "chunk_*.xml"         # custom glob pattern
    python concat_xml.py a.xml b.xml -o out.xml   # explicit file list
"""

import glob
import sys
import argparse
from lxml import etree

NAMESPACE = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas"
NS = {"ipr": NAMESPACE}


def parse_args():
    parser = argparse.ArgumentParser(description="Merge InterProScan XML files.")
    parser.add_argument(
        "files", nargs="*",
        help="Explicit list of XML files to merge (overrides --pattern)"
    )
    parser.add_argument(
        "-p", "--pattern",
        default="*.xml",
        help="Glob pattern for input files (default: '*.xml')"
    )
    parser.add_argument(
        "-o", "--output",
        default="merged_output.xml",
        help="Output filename (default: merged_output.xml)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.files:
        input_files = sorted(args.files)
    else:
        input_files = sorted(glob.glob(args.pattern))

    if not input_files:
        print(f"No files matched. Exiting.")
        sys.exit(1)

    print(f"Found {len(input_files)} file(s) to merge:")
    for f in input_files:
        print(f"  {f}")

    first_tree = etree.parse(input_files[0])
    first_root = first_tree.getroot()

    merged_root = etree.Element(first_root.tag, attrib=first_root.attrib, nsmap=first_root.nsmap)

    total_proteins = 0

    for filepath in input_files:
        tree = etree.parse(filepath)
        root = tree.getroot()

        proteins = root.findall("ipr:protein", NS)
        print(f"  {filepath}: {len(proteins)} protein(s)")
        total_proteins += len(proteins)

        for protein in proteins:
            merged_root.append(protein)

    out_tree = etree.ElementTree(merged_root)
    out_tree.write(
        args.output,
        xml_declaration=True,
        encoding="UTF-8",
        pretty_print=True
    )

    print(f"\nDone. Merged {total_proteins} protein(s) from {len(input_files)} file(s).")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()

