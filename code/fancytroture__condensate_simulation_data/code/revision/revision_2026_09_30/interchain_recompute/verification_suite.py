#!/usr/bin/env python3
"""Record tests and completed output audits with the actual software versions."""
from datetime import datetime, timezone
import importlib.metadata
import io
import json
import platform
import sys
import unittest
from analyze import HERE, SYSTEMS


def main():
    stream=io.StringIO()
    suite=unittest.defaultTestLoader.discover(str(HERE),pattern='test_*.py')
    result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    checks=[]
    for s,l,h in SYSTEMS:
        for temp in [l,h]:
            case=f'{s}_{temp}'
            checks.append(json.loads((HERE/'raw'/case/'output_audit.json').read_text()))
            independent=json.loads((HERE/f'validation_rawunits_{case}.json').read_text())
            if not all(r['match'] for r in independent):
                raise AssertionError(f'{case}: independent reader mismatch')
    payload=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
                 python=sys.version,platform=platform.platform(),
                 packages={p:importlib.metadata.version(p) for p in ['numpy','scipy','MDAnalysis','mdtraj','matplotlib']},
                 tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
                 successful=result.wasSuccessful(),raw_output_checks=checks,
                 scope='Intermolecular Figure 3 reanalysis only; not validation of all manuscript results.')
    (HERE/'verification.json').write_text(json.dumps(payload,indent=2)+'\n')
    (HERE/'test_results.txt').write_text(stream.getvalue())
    print(stream.getvalue())
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__=='__main__':
    main()
