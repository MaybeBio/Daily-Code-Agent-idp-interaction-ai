"""Supplementary figure numbering in order of first main-text citation."""
import re

ORIGINAL_TO_CURRENT = {1:1, 2:3, 3:4, 4:5, 5:6, 6:7, 7:10, 8:15}
V5_TO_CURRENT = {n: (2 if n==16 else n+1 if 2<=n<=15 else n)
                 for n in range(1,20)}
SELECTION_REORDER = {13:15, 14:16, 15:17, 16:13, 17:14}
V5_TO_CURRENT = {n: SELECTION_REORDER.get(value, value)
                 for n, value in V5_TO_CURRENT.items()}
FIGURE_REFERENCE = re.compile(
    r'(?:Figs?\.|Figures?|图)\s*~?\s*S\d+[A-F]?'
    r'(?:(?:\s*(?:--|[-–]|and|to|至|,|和)\s*)S\d+[A-F]?)*')

def renumber_references(text, mapping, changed_wrapper=None, include_labels=True):
    """Change figure references only; table numbers and page numbers stay intact."""
    def replace(match):
        before = match[0]
        after = re.sub(r'S(\d+)', lambda m:'S'+str(mapping.get(int(m[1]),int(m[1]))), before)
        return changed_wrapper(after) if changed_wrapper and after != before else after
    text = FIGURE_REFERENCE.sub(replace, text)
    return (re.sub(r'(?<=fig:S)\d+',lambda m:str(mapping.get(int(m[0]),int(m[0]))),text)
            if include_labels else text)

def original_asset_name(name):
    return re.sub(r'^FigS(\d+)(?=[A-C]?\.pdf$)',
                  lambda m:'FigS'+str(ORIGINAL_TO_CURRENT[int(m[1])]), name)
