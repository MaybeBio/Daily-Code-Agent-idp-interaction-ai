// The Overview's tests in the browser (docs/assets/idr-stats.js with the Methylarginine Site Explorer's ptm-core.js)
// against the same tests in Python (MicroDomains idr_overview_stats.py, written to docs/data/overview.json): for the
// whole dataset the two must agree, so a selection is tested exactly as the whole dataset was.
//   node --test tests/overview.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import vm from 'node:vm';
import path from 'node:path';

const DOCS = path.resolve(path.dirname(new URL(import.meta.url).pathname), '..', 'docs');
const read = (f) => JSON.parse(fs.readFileSync(path.join(DOCS, 'data', f), 'utf8'));
const ctx = { window: {}, document: { querySelector: () => null, createElement: () => ({ getContext: () => ({ measureText: (s) => ({ width: 6 * s.length }) }) }) }, console };
vm.createContext(ctx);
for (const f of ['ptm-core.js', 'idr-stats.js']) vm.runInContext(fs.readFileSync(path.join(DOCS, 'assets', f), 'utf8'), ctx);
const C = ctx.window.PTMCore;
const D = { overview: read('overview.json'), sites: read('sites.json'), residues: read('residues.json'), proteins: read('proteins.json'), idrs: read('idrs.json') };
const ST = ctx.window.IDRStats.make(C, D);
const rows = (t) => t.rows.map((r) => Object.fromEntries(t.cols.map((c, i) => [c, r[i]])));
// overview.json keeps 4 significant digits; p values below 1e-300 are 0 on both sides
const close = (a, b, what) => {
  if (a == null || b == null || !Number.isFinite(a) || !Number.isFinite(b)) { assert.equal(a == null || !Number.isFinite(a), b == null || !Number.isFinite(b), what + ': ' + a + ' vs ' + b); return; }
  if (Math.abs(a) < 1e-290 && Math.abs(b) < 1e-290) return;
  assert.ok(Math.abs(a - b) <= 6e-4 * Math.max(Math.abs(a), Math.abs(b)) + 1e-12, what + ': browser ' + a + ', Python ' + b);
};

test('residue tests: every modification in every context', () => {
  let n = 0;
  for (const w of rows(D.overview.residue_tests)) {
    if (!w.sites) continue;
    const r = ST.residueTest(w.kind, w.context, null);
    const at = w.kind + ' in ' + w.context;
    for (const k of ['sites_in', 'sites', 'other_in', 'other', 'proteins', 'mh_strata']) assert.equal(r[k], w[k], at + ' ' + k);
    assert.equal(r.tested, !!w.tested, at + ' tested');
    for (const k of ['pooled_or', 'pooled_p', 'mh_or', 'mh_lo', 'mh_hi', 'mh_p']) close(r[k], w[k], at + ' ' + k);
    n++;
  }
  assert.ok(n > 100, n + ' tests compared');
});

test('kinase-motif tests', () => {
  for (const w of rows(D.overview.kinase_tests)) {
    const b = D.overview.kinase_rows.findIndex((x) => x[0] === w.row);
    const r = ST.kinaseTest(w.context, b, null);
    const at = w.row + ' in ' + w.context;
    for (const k of ['in_context', 'in_context_match', 'outside', 'outside_match', 'mh_strata']) assert.equal(r[k], w[k], at + ' ' + k);
    assert.equal(r.tested, !!w.tested, at + ' tested');
    for (const k of ['mh_or', 'mh_lo', 'mh_hi', 'mh_p']) close(r[k], w[k], at + ' ' + k);
  }
});

test('protein classes of the proteins with acidic or basic flanks', () => {
  const W = rows(D.overview.protein_sets);
  for (const [fl, bit] of [['acidic', 0], ['basic', 1]]) {
    const res = ST.proteinSets(ST.hasProteins(D.overview.has.indexOf(fl + '_flank')));
    for (const r of res) {
      const w = W.find((x) => x.flank === fl && x.class === r.class);
      for (const k of ['k', 'n', 'M', 'N']) assert.equal(r[k], w[k], fl + ' ' + r.class + ' ' + k);
      close(r.fold, w.fold, fl + ' ' + r.class + ' fold');
      close(r.p, w.p, fl + ' ' + r.class + ' p');
      close(r.q, w.q, fl + ' ' + r.class + ' q');
    }
    assert.ok(bit >= 0);
  }
});

test('charge blocks against random placement', () => {
  for (const w of rows(D.overview.charge_shift)) {
    if (w.kind === 'glu') continue;
    const r = ST.chargeShift(w.kind, null);
    for (const k of ['idrs', 'above', 'below']) assert.equal(r[k], w[k], w.kind + ' ' + k);
    close(r.sign_p, w.sign_p, w.kind + ' sign p');
    close(r.share_delta_up, w.share_delta_up, w.kind + ' share up');
  }
});

test('occupancy inside and outside each context', () => {
  for (const w of rows(D.overview.stoichiometry)) {
    const which = w.measure.startsWith('phospho') ? 'occ' : 'ac';
    const r = ST.stoich(which, w.context, null);
    assert.equal(r.n_in, w.n_in, w.measure + ' ' + w.context + ' n_in');
    assert.equal(r.n_out, w.n_out, w.measure + ' ' + w.context + ' n_out');
    close(r.median_in, w.median_in, w.measure + ' ' + w.context + ' median');
    if (w.n_in > 20) close(r.mwu_p, w.mwu_p, w.measure + ' ' + w.context + ' Mann-Whitney p');
  }
});

test('a selection: the same tests on a subset, sites outside it left out of both sides', () => {
  const prots = new Set(ST.hasProteins(D.overview.has.indexOf('acidic_flank')));
  const r = ST.residueTest('ck2', 'acidic_flank', { prots, keep: null });
  assert.ok(r.proteins > 0 && r.proteins <= prots.size);
  const all = ST.residueTest('ck2', 'acidic_flank', null);
  assert.ok(r.sites <= all.sites);
});
