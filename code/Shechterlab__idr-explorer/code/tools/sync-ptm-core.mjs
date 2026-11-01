// Generates docs/assets/ptm-core.js and docs/assets/ptm-core.css from the Methylarginine Site Explorer, so the IDR
// Explorer draws its overview, runs its tests and exports its figures with that site's own code and styles.
//
//   node tools/sync-ptm-core.mjs --from ../PTM_Pipeline     read a checkout (its files must be those of the commit below)
//   node tools/sync-ptm-core.mjs                           fetch the files at the commit in tools/ptm-core.json
//   node tools/sync-ptm-core.mjs --check [--from DIR]      fail if the files in docs/assets differ from what it makes
//
// What is taken (tools/ptm-core.json): the declarations named in `roots` and every declaration of explorer.js they use,
// in the order explorer.js has them and with their comments, copied verbatim into one function; the names in `host`
// stop the search and are supplied by the IDR Explorer with PTMCore.use({...}) (its data, its selection, its links and
// storage). The CSS is every rule of explorer.css whose classes all occur in the copied code or in `classes` (markup the
// IDR Explorer writes itself in the explorer's styles), the custom properties those rules use, and their keyframes. To
// take more, add a root or a class and run this again; never edit the generated files.
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import * as acorn from 'acorn';
import * as walk from 'acorn-walk';
import postcss from 'postcss';

const HERE = path.dirname(new URL(import.meta.url).pathname);
const REPO = path.resolve(HERE, '..');
const cfg = JSON.parse(fs.readFileSync(path.join(HERE, 'ptm-core.json'), 'utf8'));
const args = process.argv.slice(2);
const CHECK = args.includes('--check');
const FROM = args.includes('--from') ? args[args.indexOf('--from') + 1] : null;
const OUT_JS = path.join(REPO, 'docs/assets/ptm-core.js');
const OUT_CSS = path.join(REPO, 'docs/assets/ptm-core.css');

async function read(rel) {
  if (FROM) return fs.readFileSync(path.join(FROM, rel), 'utf8');
  const url = `https://raw.githubusercontent.com/${cfg.source.repo}/${cfg.source.commit}/${rel}`;
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return r.text();
}
const md5 = (s) => crypto.createHash('md5').update(s, 'utf8').digest('hex');

// ---- JavaScript
function extractJs(src) {
  const comments = [];
  const ast = acorn.parse(src, { ecmaVersion: 'latest', sourceType: 'script', onComment: comments, locations: true });
  const iife = ast.body.find((s) => s.type === 'ExpressionStatement' && s.expression.type === 'CallExpression' && s.expression.callee.type === 'FunctionExpression');
  if (!iife) throw new Error('explorer.js: no top-level function expression found');
  const body = iife.expression.callee.body.body;
  const names = (pat, out = []) => {
    if (!pat) return out;
    if (pat.type === 'Identifier') out.push(pat.name);
    else if (pat.type === 'ObjectPattern') pat.properties.forEach((p) => names(p.type === 'RestElement' ? p.argument : p.value, out));
    else if (pat.type === 'ArrayPattern') pat.elements.forEach((e) => names(e, out));
    else if (pat.type === 'AssignmentPattern') names(pat.left, out);
    else if (pat.type === 'RestElement') names(pat.argument, out);
    return out;
  };
  const decl = new Map();
  for (const st of body) {
    if (st.type === 'FunctionDeclaration' || st.type === 'ClassDeclaration') decl.set(st.id.name, st);
    else if (st.type === 'VariableDeclaration') for (const d of st.declarations) for (const n of names(d.id)) decl.set(n, st);
  }
  // identifiers a statement refers to (not property names or object keys); inits: only those evaluated when the
  // statement runs (not inside a function)
  const refs = (node, inits) => {
    const skip = new Set(), out = new Set();
    walk.full(node, (n) => {
      if (n.type === 'MemberExpression' && !n.computed) skip.add(n.property);
      if ((n.type === 'Property' || n.type === 'MethodDefinition' || n.type === 'PropertyDefinition') && !n.computed && !n.shorthand) skip.add(n.key);
      if (n.type === 'LabeledStatement' || n.type === 'BreakStatement' || n.type === 'ContinueStatement') skip.add(n.label);
    });
    const visit = (root) => walk.full(root, (n) => { if (n.type === 'Identifier' && !skip.has(n)) out.add(n.name); });
    if (!inits) { visit(node); return out; }
    walk.recursive(node, null, {
      Function() { /* not evaluated now */ },
      Identifier(n) { if (!skip.has(n)) out.add(n.name); },
    });
    return out;
  };
  const host = new Set(cfg.host), want = new Set(), queue = [...cfg.roots];
  while (queue.length) {
    const k = queue.pop();
    if (want.has(k) || host.has(k)) continue;
    const st = decl.get(k);
    if (!st) throw new Error(`explorer.js: no top-level declaration named ${k}`);
    want.add(k);
    for (const r of refs(st)) if (decl.has(r) && !want.has(r)) queue.push(r);
  }
  const stmts = [...new Set([...want].map((k) => decl.get(k)))].sort((a, b) => a.start - b.start);
  // a host name read while the file loads would be undefined then
  for (const st of stmts) {
    if (st.type !== 'VariableDeclaration') continue;
    const early = [...refs(st, true)].filter((r) => host.has(r));
    if (early.length) throw new Error(`explorer.js line ${st.loc.start.line}: reads ${early.join(', ')} when it loads`);
  }
  // every name a copied statement declares is exported (and none may be a host name)
  const exported = [];
  for (const st of stmts) {
    const ns = st.type === 'VariableDeclaration' ? st.declarations.flatMap((d) => names(d.id)) : [st.id.name];
    for (const n of ns) { if (host.has(n)) throw new Error(`${n} is both copied and a host name`); exported.push(n); }
  }
  // each statement with the comment lines directly above it
  const chunks = stmts.map((st) => {
    const i = body.indexOf(st), lo = i > 0 ? body[i - 1].end : 0;
    const lead = comments.filter((c) => c.start >= lo && c.end <= st.start);
    // only the block of comments that touches the statement (no blank line between)
    let start = st.start;
    for (let k = lead.length - 1; k >= 0; k--) {
      const gap = src.slice(lead[k].end, start);
      if (/\n\s*\n/.test(gap)) break;
      start = lead[k].start;
    }
    return '  ' + src.slice(start, st.end);
  });
  // every string the copied code holds, for the CSS: its class names
  const strings = new Set();
  for (const st of stmts) walk.full(st, (n) => {
    const s = n.type === 'Literal' && typeof n.value === 'string' ? n.value : n.type === 'TemplateElement' ? n.value.cooked : null;
    if (s) for (const t of s.split(/[^A-Za-z0-9_-]+/)) if (t) strings.add(t);
  });
  return { chunks, exported, strings, count: stmts.length };
}

// ---- CSS
async function extractCss(css, strings) {
  const root = postcss.parse(css);
  const isTokenRule = (rule) => rule.selectors.every((s) => /^:root\b/.test(s.trim()) && !/[.#]/.test(s.replace(/\[[^\]]*\]/g, '')));
  const used = new Set(), anims = new Set();
  const selectorOk = (s) => {
    const bare = s.replace(/\[[^\]]*\]/g, '').replace(/:(not|is|where|has)\([^)]*\)/g, '');
    const toks = [...bare.matchAll(/[.#]([A-Za-z_][A-Za-z0-9_-]*)/g)].map((m) => m[1]);
    return toks.length > 0 && toks.every((t) => strings.has(t));
  };
  root.walkRules((rule) => {
    if (rule.parent && rule.parent.type === 'atrule' && /keyframes$/.test(rule.parent.name)) return;
    if (isTokenRule(rule)) return;
    const keep = rule.selectors.filter(selectorOk);
    if (!keep.length) { rule.remove(); return; }
    rule.selectors = keep;
    rule.walkDecls((d) => {
      for (const m of d.value.matchAll(/var\(\s*(--[A-Za-z0-9_-]+)/g)) used.add(m[1]);
      if (/^animation(-name)?$/.test(d.prop)) for (const t of d.value.split(/[\s,]+/)) anims.add(t);
    });
  });
  // custom properties: those the kept rules use, and those they are made from
  const defs = new Map();
  root.walkRules((rule) => { if (isTokenRule(rule)) rule.walkDecls((d) => { if (d.prop.startsWith('--')) { if (!defs.has(d.prop)) defs.set(d.prop, []); defs.get(d.prop).push(d.value); } }); });
  let grew = true;
  while (grew) {
    grew = false;
    for (const p of [...used]) for (const v of defs.get(p) || []) for (const m of v.matchAll(/var\(\s*(--[A-Za-z0-9_-]+)/g)) if (!used.has(m[1])) { used.add(m[1]); grew = true; }
  }
  root.walkRules((rule) => {
    if (!isTokenRule(rule)) return;
    rule.walkDecls((d) => { if (d.prop.startsWith('--') && !used.has(d.prop)) d.remove(); });
    if (!rule.nodes.some((n) => n.type === 'decl')) rule.remove();
  });
  root.walkAtRules((at) => { if (/keyframes$/.test(at.name) && !anims.has(at.params.trim())) at.remove(); });
  let changed = true;
  while (changed) {
    changed = false;
    root.walkAtRules((at) => { if (at.nodes && !at.nodes.some((n) => n.type !== 'comment')) { at.remove(); changed = true; } });
  }
  root.walkComments((c) => { if (!c.next() || c.next().type === 'comment') c.remove(); });
  return root.toString().replace(/\n{3,}/g, '\n\n').trim() + '\n';
}

const js = await read(cfg.source.js), css = await read(cfg.source.css);
const { chunks, exported, strings, count } = extractJs(js);
for (const c of cfg.classes || []) strings.add(c);
const stamp = `${cfg.source.repo} ${cfg.source.js} at ${cfg.source.commit} (MD5 ${md5(js)})`;
const outJs = `/* GENERATED by tools/sync-ptm-core.mjs from ${stamp}. Do not edit: add a root to tools/ptm-core.json, or change
   the Methylarginine Site Explorer, and run it again.

   The charts, effect strips, scorecard, tooltips, exports and statistics of the Methylarginine Site Explorer (${count}
   declarations, copied verbatim), so the IDR Explorer draws, tests and exports its overview exactly as that site does.
   Supplied by the IDR Explorer with PTMCore.use({...}): ${cfg.host.join(', ')}. */
window.PTMCore = (function () {
  'use strict';
  let ${cfg.host.join(', ')};

${chunks.join('\n\n')}

  function use(o) {
${cfg.host.map((h) => `    if ('${h}' in o) ${h} = o.${h};`).join('\n')}
  }
  return { use, ${exported.join(', ')} };
})();
`;
const outCss = `/* GENERATED by tools/sync-ptm-core.mjs from ${cfg.source.repo} ${cfg.source.css} at ${cfg.source.commit} (MD5 ${md5(css)}):
   the rules of the Methylarginine Site Explorer for what docs/assets/ptm-core.js draws. Do not edit; idr.css, loaded after
   this file, holds the IDR Explorer's own styles. */
${await extractCss(css, strings)}`;

if (CHECK) {
  const bad = [[OUT_JS, outJs], [OUT_CSS, outCss]].filter(([f, s]) => !fs.existsSync(f) || fs.readFileSync(f, 'utf8') !== s);
  if (bad.length) { console.error('Out of date (run node tools/sync-ptm-core.mjs): ' + bad.map(([f]) => path.relative(REPO, f)).join(', ')); process.exit(1); }
  console.log(`ptm-core.js and ptm-core.css match ${stamp}`);
} else {
  fs.writeFileSync(OUT_JS, outJs);
  fs.writeFileSync(OUT_CSS, outCss);
  console.log(`wrote ${path.relative(REPO, OUT_JS)} (${count} declarations, ${outJs.length.toLocaleString()} bytes) and ${path.relative(REPO, OUT_CSS)} (${outCss.length.toLocaleString()} bytes) from ${stamp}`);
}
