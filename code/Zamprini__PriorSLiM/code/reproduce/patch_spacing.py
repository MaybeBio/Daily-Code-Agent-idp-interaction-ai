"""Apply Figure 2/4 spacing changes to the specified input deck."""
import argparse, copy, hashlib, json, os, shutil, subprocess, tempfile
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from lxml import etree as E
SOURCE_SHA256 = '3cc4169cf0cec73dd34799d3c2dea03e0b7d17a69a1831d75ba699430fc31764'
KNOWN_PATCHED_SHA256 = '1c195e8cc40e7b22732de617ba585aff12d1f350dd198591cbc01803b82563ae'

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path)
    ap.add_argument('output', nargs='?', type=Path)
    ap.add_argument('--check', action='store_true', help='Verify exact input identity only; no changes.')
    ap.add_argument('--arial-font-dir', type=Path, help='Directory containing user-licensed Arial TTF/OTF files.')
    ap.add_argument('--inkscape', default='inkscape', help='Inkscape executable or absolute path.')
    a = ap.parse_args()
    digest = hashlib.sha256(a.input.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        detail = 'already patched' if digest == KNOWN_PATCHED_SHA256 else 'unknown or modified version'
        ap.error(f'Refusing {detail}. Input SHA256: {digest}. No files changed.')
    if a.check:
        print('SUPPORTED: input deck checksum matches. No files written.')
        return
    if a.output is None or a.arial_font_dir is None:
        ap.error('output and --arial-font-dir are required unless using --check')
    if a.output.exists() or a.output.resolve() == a.input.resolve():
        ap.error('Output must be a new file; overwriting input or existing output is forbidden')
    inkscape = shutil.which(a.inkscape)
    if not inkscape:
        ap.error('Inkscape executable not found')
    from fontTools.ttLib import TTFont
    fontdir = a.arial_font_dir.resolve()
    arial = False
    for f in fontdir.rglob('*'):
        if f.suffix.lower() not in ('.ttf', '.otf'):
            continue
        try:
            with TTFont(f) as font:
                families = {n.toUnicode() for n in font['name'].names if n.nameID in (1, 16)}
                arial = arial or 'Arial' in families
        except Exception:
            pass
    if not arial:
        ap.error('No Arial font family found in --arial-font-dir; font substitution is not permitted')
    with tempfile.TemporaryDirectory(prefix='figure-spacing-') as tmp:
        work = Path(tmp)
        stage = work / 'pptx'
        with ZipFile(a.input) as z:
            for item in z.infolist():
                p = (stage / item.filename).resolve()
                if not p.is_relative_to(stage.resolve()):
                    ap.error('Unsafe archive entry')
            z.extractall(stage)
        config = E.Element('fontconfig')
        E.SubElement(config, 'dir').text = str(fontdir)
        E.SubElement(config, 'cachedir').text = str(work / 'fontcache')
        fc = work / 'fonts.conf'
        fc.write_bytes(E.tostring(config, xml_declaration=True, encoding='UTF-8'))
        os.environ['FONTCONFIG_FILE'] = str(fc)
        D = stage / 'ppt/media'
        NS = 'http://www.w3.org/2000/svg'
        changes = {}

        def shift(e, dx=0, dy=0):
            e.set('transform', f'translate({dx} {dy}) ' + e.get('transform', ''))

        def find(r, id):
            return r.xpath('//*[@id=$v]', v=id)[0]

        def texts(r):
            return r.findall('.//{' + NS + '}text')

        def expand(r, extra):
            a = list(map(float, r.get('viewBox').split()))
            a[3] += extra
            r.set('viewBox', ' '.join(map(str, a)))
            r.set('height', str(a[3]) + 'pt' if 'pt' in r.get('height', '') else str(a[3]))
            return a[3] / (a[3] - extra)
        for name in ['image4.svg', 'image5.svg', 'image7.svg', 'image8.svg', 'image9.svg', 'retyped_4_0.svg', 'retyped_4_1.svg', 'retyped_4_3.svg']:
            r = E.parse(str(D / name)).getroot()
            factor = 1
            for t in texts(r):
                if name in ['image4.svg', 'image5.svg']:
                    content = ''.join(t.itertext())
                    for child in list(t):
                        t.remove(child)
                    t.text = content
            if name == 'image4.svg':
                for i, lines in enumerate([['Prior', 'SLiM'], ['SLiM', 'Mine'], ['SHARK'], ['ELM']]):
                    t = find(r, f'fonttext_{i}')
                    t.attrib.pop('transform', None)
                    t.set('x', str(108.81 + i * 56.34))
                    t.set('y', '306')
                    t.set('text-anchor', 'middle')
                    t.set('font-size', str(18.654006304 * 16 / 18))
                    t.text = lines[0]
                    if len(lines) > 1:
                        c = copy.deepcopy(t)
                        c.attrib.pop('id', None)
                        c.set('y', '329')
                        c.text = lines[1]
                        t.getparent().append(c)
                for i in range(4, 11):
                    t = find(r, f'fonttext_{i}')
                    t.set('font-size', str(18.653999781 * 16 / 18))
                    t.set('x', '72')
                    t.attrib.pop('transform', None)
            if name == 'image5.svg':
                for i in range(8):
                    t = find(r, f'fonttext_{i}')
                    t.set('font-size', str(18.083221281 * 16 / 18))
                    shift(t, dy=8 if i % 2 == 0 else 14)
                for i in range(8, 18):
                    t = find(r, f'fonttext_{i}')
                    t.set('font-size', str(18.083221281 * 16 / 18))
                    t.set('x', '193')
            if name == 'image7.svg':
                for t in texts(r):
                    s = ''.join(t.itertext()).strip()
                    if s == 'Threshold':
                        shift(t, dy=7)
                factor = expand(r, 8)
            if name == 'image8.svg':
                shift(find(r, 'fonttext_9'), dy=9)
                factor = expand(r, 12)
            if name == 'image9.svg':
                for i in [24, 26]:
                    shift(find(r, f'fonttext_{i}'), dy=7)
                shift(find(r, 'fonttext_16'), dy=7)
                leg = find(r, 'legend_1')
                for e in list(leg):
                    if e.get('id') != 'text_25':
                        shift(e, dy=14)
                factor = expand(r, 18)
            if name == 'retyped_4_0.svg':
                shift(find(r, 'legend_1'), dy=9)
            if name == 'retyped_4_1.svg':
                shift(find(r, 'legend_1'), dy=9)
                factor = expand(r, 12)
            if name == 'retyped_4_3.svg':
                shift(find(r, 'legend_1'), dy=12)
                factor = expand(r, 18)
                for t in texts(r):
                    s = ''.join(t.itertext()).strip()
                    if s in ['Human IDRs (ASHI)', 'Phospho-screen WT', 'Viral IDRs', 'Swiss-Prot MOTIF']:
                        shift(t, dy=-4)
                    if s == 'loci (%)':
                        shift(t, dy=3)
            rbytes = E.tostring(r, xml_declaration=True, encoding='UTF-8')
            (D / name).write_bytes(rbytes)
            changes[name] = factor
        ns = {'p': 'http://schemas.openxmlformats.org/presentationml/2006/main', 'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}
        R = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed'
        for n in [2, 4]:
            p = stage / f'ppt/slides/slide{n}.xml'
            r = E.parse(str(p))
            rels = E.parse(str(p.parent / '_rels' / f'slide{n}.xml.rels'))
            rd = {e.get('Id'): Path(e.get('Target')).name for e in rels.getroot()}
            for pic in r.findall('.//p:pic', ns):
                ids = [e.get(R) for e in pic.iter() if e.get(R)]
                names = [rd[i] for i in ids]
                svg = next((x for x in names if x in changes), None)
                if svg:
                    ext = pic.find('.//a:xfrm/a:ext', ns)
                    ext.set('cy', str(round(int(ext.get('cy')) * changes[svg])))
                    png = next((x for x in names if x.endswith('.png')))
                    subprocess.run([inkscape, str(D / svg), '--export-type=png', '--export-width=2600', '--export-filename=' + str(D / png)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                if n == 4 and any((x in names for x in ['retyped_4_2.svg', 'retyped_4_3.svg'])):
                    off = pic.find('.//a:xfrm/a:off', ns)
                    off.set('y', str(int(off.get('y')) + 317500))
            if n == 4:
                for sp in r.findall('.//p:sp', ns):
                    txt = ''.join(sp.xpath('.//a:t/text()', namespaces=ns))
                    if txt in ['d', 'e']:
                        off = sp.find('.//a:xfrm/a:off', ns)
                        off.set('y', str(int(off.get('y')) + 317500))
            p.write_bytes(E.tostring(r, xml_declaration=True, encoding='UTF-8', standalone=True))
        a.output.parent.mkdir(parents=True, exist_ok=True)
        with ZipFile(a.output, 'x', ZIP_DEFLATED) as z:
            for p in sorted(stage.rglob('*')):
                if p.is_file():
                    z.write(p, p.relative_to(stage))
        print(json.dumps({'input_sha256': digest, 'output': str(a.output), 'changes': changes}, indent=2))
if __name__ == '__main__':
    main()
