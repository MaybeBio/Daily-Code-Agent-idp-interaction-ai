"""
xlsx_to_csv.py — standard-library conversion of every sheet of an .xlsx file to CSV (no openpyxl needed). Cell values
are copied as text, shared strings resolved; inline strings and numbers kept as written in the file.

  python3 -I xlsx_to_csv.py <file.xlsx> <out_dir>   ->   <out_dir>/<stem>__<sheet name>.csv
"""
import csv
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
RELS = '{http://schemas.openxmlformats.org/package/2006/relationships}'


def col_index(ref):
    letters = re.match(r'[A-Z]+', ref).group(0)
    k = 0
    for ch in letters:
        k = k * 26 + ord(ch) - 64
    return k - 1


def main(path, out):
    os.makedirs(out, exist_ok=True)
    z = zipfile.ZipFile(path)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS):
            shared.append(''.join(t.text or '' for t in si.iter(f'{{{NS["m"]}}}t')))
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels')).iter(f'{RELS}Relationship')}
    stem = os.path.splitext(os.path.basename(path))[0]
    for sh in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets', NS):
        target = rels[sh.get(f'{{{NS["r"]}}}id')].lstrip('/')
        target = target if target.startswith('xl/') else 'xl/' + target
        rows = []
        for row in ET.fromstring(z.read(target)).iter(f'{{{NS["m"]}}}row'):
            vals = {}
            for c in row.findall('m:c', NS):
                t, v = c.get('t'), c.find('m:v', NS)
                if t == 's' and v is not None:
                    val = shared[int(v.text)]
                elif t == 'inlineStr':
                    val = ''.join(x.text or '' for x in c.iter(f'{{{NS["m"]}}}t'))
                else:
                    val = v.text if v is not None else ''
                vals[col_index(c.get('r'))] = val
            rows.append([vals.get(k, '') for k in range(max(vals) + 1)] if vals else [])
        name = re.sub(r'[^A-Za-z0-9_.-]+', '_', sh.get('name'))
        dest = os.path.join(out, f'{stem}__{name}.csv')
        with open(dest, 'w', newline='') as fh:
            csv.writer(fh).writerows(rows)
        print(f'{dest}: {len(rows)} rows', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
