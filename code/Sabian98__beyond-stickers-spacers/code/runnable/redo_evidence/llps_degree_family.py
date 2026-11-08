"""
llps_degree_family.py — llps_degree_plan.md: reads the PRIMARY line of the three degree summaries, applies Holm's procedure
across the three P values, and writes llps_csat/llps_degree_family.txt next to the occupancy results of the plans.
  ~/.venvs/bender/bin/python llps_degree_family.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'llps_csat')
TESTS = [('csat of hnRNPA1 LCD variants (llps_csat_plan.md)', 'llps_csat_summary_degree.txt', 'occupancy P 0.878'),
         ('PhaSePro drivers, local-shuffle unevenness (llps_unevenness_plan.md)', 'llps_unevenness_summary_degree.txt',
          'occupancy P 0.64'),
         ('disease mutations, |local response| (llps_mutation_plan.md)', 'llps_mutations_degree.txt', 'occupancy P 0.86')]


def primary(path):
    line = next(l for l in open(path) if l.startswith('PRIMARY'))
    return line.strip(), float(re.findall(r'P = ([0-9.eE+-]+)', line)[-1])


def main():
    rows = [(name, *primary(os.path.join(D, f)), occ) for name, f, occ in TESTS]
    ps = [r[2] for r in rows]
    order = sorted(range(3), key=lambda i: ps[i])
    adj, run = [0.0] * 3, 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (3 - k) * ps[i]))
        adj[i] = run
    L = ['LLPS tests with predicted contact degree (llps_degree_plan.md), Holm across the three primaries', '']
    for (name, line, p, occ), a in zip(rows, adj):
        L += [f'{name}', f'  {line}', f'  raw P {p:.4g}, Holm P {a:.4g} -> {"POSITIVE (lead, needs replication)" if a < 0.05 else "not positive"}'
              f'   [{occ}]', '']
    txt = '\n'.join(L)
    open(os.path.join(D, 'llps_degree_family.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
