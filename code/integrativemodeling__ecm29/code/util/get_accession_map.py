# Create the accession map used by test/mock/sitecustomize.py

import ihm.reference

codes = [ 'P62195', 'P62333', 'P17980', 'P62191', 'P43686', 'P35998', 'P48556',
          'P55036', 'O00487', 'P60896', 'Q13200', 'Q99460', 'O43242', 'O00232',
          'O00231', 'Q15008', 'P51665', 'Q9UNM6', 'Q5VYK3']

def pp(s):
    indent = 8
    width = 66
    def get_lines(s):
        for i in range(0, len(s), width):
            yield ' ' * indent + "'" + s[i:i+width] + "'"
    return '\n'.join(l for l in get_lines(s))

for code in codes:
    u = ihm.reference.UniProtSequence.from_accession(code)
    print("    '%s': {'db_code':'%s', 'sequence':\n%s},"
          % (code, u.db_code, pp(u.sequence)))
