"""Build a hand-adjudication worksheet from the v1 classifications.

The two models disagree on 15% of file-chunks, almost all in one direction, and
revising the prompt against my own reading of those disagreements would just
substitute my judgment for the models'. This produces a stratified sample for
the author to adjudicate instead.

Design notes
------------
The unit is a file-chunk, not a paper, because that is where the evidence lines
live and it makes each decision a single snippet rather than a whole archive.

Evidence lines are only about 78% verbatim --- the models sometimes truncate or
merge lines --- so each item shows the REAL file context around the closest
matching line, not just the quoted text. Adjudicating a slightly-wrong quote
would be worse than useless.

Strata are deliberately unequal. Both-agree-no is the largest population but the
least informative; llama-only is where the suspected over-flagging lives.

Writes:
  validation-worksheet.html   read this, click through the items
  validation-verdicts.csv     fallback if the in-page recording fails
"""
from __future__ import annotations

import argparse
import difflib
import html
import json
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_HTML = ROOT / 'validation-worksheet.html'
OUT_CSV = ROOT / 'validation-verdicts.csv'
SEED = 20260820

# how many items from each stratum
QUOTA = {
    'a_only': 35,       # what run A finds and run B does not
    'both_yes': 30,     # precision on agreed positives
    'b_only': 20,
    'both_no': 15,      # recall check
}
CONTEXT = 4  # lines either side of the matched line


def load(tag, prompt):
    rows = {}
    path = ROOT / 'llm-output' / f'{tag}.{prompt}.jsonl'
    with path.open(encoding='utf-8') as fh:
        for line in fh:
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not d.get('result'):
                continue
            rows[(d['file_id'], d['chunk'])] = d  # later attempt wins
    return rows


def norm(s):
    return ' '.join(s.split()).lower()


def find_context(path: Path, evidence: list[str]):
    """Locate the closest real line to the quoted evidence, with context."""
    try:
        text = path.read_bytes().decode('utf-8', errors='replace')
    except OSError:
        return None, []
    lines = text.splitlines()
    if not evidence:
        return None, []
    target = norm(evidence[0])
    if not target:
        return None, []

    best_i, best_score = None, 0.0
    for i, ln in enumerate(lines):
        n = norm(ln)
        if not n:
            continue
        if target in n or n in target:
            score = 1.0
        else:
            score = difflib.SequenceMatcher(None, target[:160], n[:160]).ratio()
        if score > best_score:
            best_i, best_score = i, score
    if best_i is None:
        return None, []
    lo = max(0, best_i - CONTEXT)
    hi = min(len(lines), best_i + CONTEXT + 1)
    return best_score, [(j + 1, lines[j], j == best_i) for j in range(lo, hi)]


def model_call_digest(path: Path, limit=14):
    """For agreed-negative files: show the lines that look like model fits, so
    a missed model is visible without reading the whole file. Display only."""
    try:
        text = path.read_bytes().decode('utf-8', errors='replace')
    except OSError:
        return []
    keys = ('reg ', 'regress', 'logit', 'probit', 'areg', 'xtreg', 'reghdfe',
            'mixed', 'lm(', 'glm(', 'lmer(', 'glmer(', 'brm(', 'feols(',
            'felm(', 'margins', 'ivreg', 'plm(', 'estimates')
    out = []
    for i, ln in enumerate(text.splitlines()):
        low = ln.lower().strip()
        if low.startswith(('*', '#', '//')):
            continue
        if any(k in low for k in keys):
            out.append((i + 1, ln))
        if len(out) >= limit:
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--a', default='llama:classify-v1',
                    help='first run as model:prompt')
    ap.add_argument('--b', default='gpt-oss:classify-v1',
                    help='second run as model:prompt')
    ap.add_argument('--field', default='het',
                    choices=['het', 'subgroup_split', 'interaction_term'],
                    help='what to compare: the overall judgment, or one method')
    ap.add_argument('--out', default='validation-worksheet',
                    help='basename for the html and csv')
    ap.add_argument('--quota', default=None,
                    help='override sample sizes, e.g. '
                         '"a_only=30,both_yes=12,b_only=8,both_no=5"')
    args = ap.parse_args()

    quota = dict(QUOTA)
    if args.quota:
        for part in args.quota.split(','):
            k, v = part.split('=')
            if k.strip() not in quota:
                raise SystemExit(f'unknown stratum {k!r}')
            quota[k.strip()] = int(v)

    global OUT_HTML, OUT_CSV
    OUT_HTML = ROOT / f'{args.out}.html'
    OUT_CSV = ROOT / f'{args.out}.csv'

    L = load(*args.a.split(':'))
    G = load(*args.b.split(':'))
    keys = sorted(set(L) & set(G))
    print(f'file-chunks in both runs: {len(keys):,}')
    print(f'comparing {args.a} (A) against {args.b} (B) on {args.field}')

    def yes(d):
        r = d['result']
        if args.field == 'het':
            return bool(r.get('estimates_heterogeneous_effects'))
        return args.field in (r.get('methods') or [])

    strata = {k: [] for k in quota}
    for k in keys:
        l, g = yes(L[k]), yes(G[k])
        s = ('both_yes' if l and g else 'both_no' if not l and not g
             else 'a_only' if l else 'b_only')
        strata[s].append(k)
    print('population:', {k: len(v) for k, v in strata.items()})

    rng = random.Random(SEED)
    picked = []
    for s, n in quota.items():
        pool = strata[s]
        picked += [(s, k) for k in rng.sample(pool, min(n, len(pool)))]
    rng.shuffle(picked)  # so the adjudicator cannot infer the stratum
    print(f'sampled {len(picked)} items')

    items = []
    for idx, (stratum, key) in enumerate(picked, start=1):
        ld, gd = L[key], G[key]
        safe = ld['dataset_doi'].replace(':', '_').replace('/', '_')
        path = ROOT / 'replication-files' / safe / ld['filename']

        ev = (ld['result'].get('evidence') or []) + \
             (gd['result'].get('evidence') or [])
        score, ctx = find_context(path, ev)
        digest = model_call_digest(path) if not ctx else []

        items.append({
            'id': idx, 'stratum': stratum,
            'doi': ld['dataset_doi'], 'filename': ld['filename'],
            'language': ld['language'], 'path': str(path),
            'chunk': f"{ld['chunk'] + 1}/{ld['n_chunks']}",
            'a_yes': yes(ld), 'b_yes': yes(gd),
            'a_methods': ld['result'].get('methods') or [],
            'b_methods': gd['result'].get('methods') or [],
            'a_ev': ld['result'].get('evidence') or [],
            'b_ev': gd['result'].get('evidence') or [],
            'match_score': score, 'context': ctx, 'digest': digest,
        })

    write_html(items)
    with OUT_CSV.open('w', encoding='utf-8', newline='') as fh:
        fh.write('id,stratum,filename,verdict\n')
        for it in items:
            fh.write(f"{it['id']},{it['stratum']},"
                     f"\"{it['filename']}\",\n")
    print(f'wrote {OUT_HTML}')
    print(f'wrote {OUT_CSV}')
    print('\nstrata in sample:', Counter(i['stratum'] for i in items))


CHEAT = """
<h2>What counts</h2>
<p>The question for every item is the same: <b>does this file fit at least one
model in which the effect of some variable is allowed to differ across units,
subgroups, or contexts?</b> Judge the code that runs, not comments.</p>
<div class="two">
<div><h3>Counts</h3><ul>
<li><code>y ~ x*z</code>, <code>x:z</code> &mdash; product term in a formula</li>
<li>Stata <code>c.x#c.z</code>, <code>i.x#i.z</code>, <code>x##z</code></li>
<li><code>gen xz = x*z</code> then <code>reg y xz x z</code> &mdash; an
    interaction built by hand</li>
<li>The same model fit on two or more contrasting subsets to compare the
    effect: <code>reg y x if female==1</code> <i>and</i>
    <code>reg y x if female==0</code></li>
<li>Varying <b>slope</b>: <code>(1 + x | group)</code>,
    <code>(0 + x | group)</code>, Stata <code>|| g: x</code></li>
<li>Causal forest, <code>bartc</code>, <code>bcf</code>, meta-learners</li>
</ul></div>
<div><h3>Does not count</h3><ul>
<li><code>feols(y ~ x | country + year)</code> &mdash; the pipe separates fixed
    effects; there is no interaction</li>
<li><code>areg</code>, <code>xtreg</code>, <code>i.country</code> as controls</li>
<li>Varying <b>intercept</b> only: <code>(1 | group)</code>, Stata
    <code>|| cntry:</code> with nothing after the colon</li>
<li>One <code>subset()</code> or <code>if</code> defining the analysis sample,
    with no comparison across subsets</li>
<li><code>margins</code> / <code>slopes()</code> on a model with no
    interaction</li>
<li>A variable whose <i>name</i> suggests a product but which is never actually
    constructed as one (e.g. <code>gen foo_x_bar = .</code>)</li>
<li><code>sum</code>, <code>ttest</code>, <code>merge</code>, data cleaning</li>
</ul></div></div>
"""


def write_html(items):
    def esc(s):
        return html.escape(str(s))

    parts = ["""<meta charset="utf-8"><title>Validation worksheet</title>
<style>
 body{font:15px/1.5 -apple-system,Segoe UI,sans-serif;max-width:1000px;
      margin:2rem auto;padding:0 1rem;color:#111}
 h1{margin-bottom:.2rem} h2{margin-top:2rem}
 .two{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem}
 .two ul{padding-left:1.1rem} .two li{margin:.3rem 0}
 code{background:#f2f2f2;padding:1px 4px;border-radius:3px;font-size:13px}
 .item{border:1px solid #ddd;border-radius:6px;padding:1rem;margin:1.5rem 0}
 .hdr{display:flex;justify-content:space-between;align-items:baseline;
      gap:1rem;flex-wrap:wrap}
 .meta{color:#666;font-size:13px}
 pre{background:#fafafa;border:1px solid #eee;border-radius:4px;padding:.6rem;
     overflow-x:auto;font-size:13px;line-height:1.45;margin:.5rem 0}
 .hit{background:#fff3bf;display:block}
 .ln{color:#999;user-select:none}
 .said{font-size:13px;color:#444;margin:.4rem 0}
 .q{margin-top:.7rem;font-size:14px}
 label{margin-right:1rem;cursor:pointer}
 .done{border-color:#7ab97a;background:#f7fdf7}
 #bar{position:sticky;top:0;background:#fff;border-bottom:1px solid #ddd;
      padding:.6rem 0;margin-bottom:1rem;z-index:9}
 button{font:inherit;padding:.35rem .8rem;cursor:pointer}
 textarea{width:100%;height:9rem;font-family:monospace;font-size:12px}
 .warn{color:#a33;font-size:13px;margin-top:.4rem}
</style>
<h1>Heterogeneous effects &mdash; hand adjudication</h1>
<p class="meta">Model labels are hidden until you answer, so the sample reads
blind. Answers save in this browser as you go.</p>
"""]
    parts.append(CHEAT)
    parts.append('<div id="bar"><b><span id="n">0</span> / %d answered</b> '
                 '&nbsp; <button onclick="dump()">Show results</button>'
                 '<div id="warn" class="warn"></div>'
                 '<div id="out"></div></div>' % len(items))

    for it in items:
        parts.append(f'<div class="item" id="i{it["id"]}">')
        parts.append(
            f'<div class="hdr"><b>#{it["id"]}</b>'
            f'<span class="meta">{esc(it["language"])} &middot; '
            f'{esc(it["filename"])} &middot; chunk {it["chunk"]}</span></div>')

        if it['context']:
            if it['match_score'] is not None and it['match_score'] < 0.95:
                parts.append('<div class="meta">The quoted line was not exact; '
                             'the closest real line is highlighted.</div>')
            parts.append('<pre>')
            for lineno, txt, is_hit in it['context']:
                row = (f'<span class="ln">{lineno:>5}</span>  '
                       f'{esc(txt[:400])}')
                parts.append(f'<span class="hit">{row}</span>' if is_hit
                             else row + '\n')
            parts.append('</pre>')
        elif it['digest']:
            parts.append('<div class="meta">Neither model flagged this file. '
                         'Lines that look like model fits:</div><pre>')
            for lineno, txt in it['digest']:
                parts.append(f'<span class="ln">{lineno:>5}</span>  '
                             f'{esc(txt[:400])}\n')
            parts.append('</pre>')
        else:
            parts.append('<div class="meta">No model-fitting lines found in '
                         'this file.</div>')

        parts.append(f'<div class="meta">Full file: <code>'
                     f'{esc(it["path"])}</code></div>')

        said = []
        for tag, y, m in (('A', it['a_yes'], it['a_methods']),
                          ('B', it['b_yes'], it['b_methods'])):
            said.append(f'model {tag}: {"yes" if y else "no"}'
                        + (f' ({", ".join(m)})' if m else ''))
        parts.append(f'<div class="said reveal" style="display:none">'
                     f'{esc(" &middot; ".join(said))}</div>'.replace(
                         '&amp;middot;', '&middot;'))

        parts.append(
            '<div class="q">Does this file fit a model letting an effect vary? '
            f'<label><input type="radio" name="v{it["id"]}" value="yes"> Yes</label>'
            f'<label><input type="radio" name="v{it["id"]}" value="no"> No</label>'
            f'<label><input type="radio" name="v{it["id"]}" value="unclear"> '
            'Can\'t tell</label></div></div>')

    parts.append("""
<script>
// localStorage is blocked on file:// origins in some browsers, so every read
// and write is guarded and the export below reads the DOM rather than the
// store. Losing an hour of adjudication to a silent storage failure would be
// considerably worse than losing the ability to resume.
const K='hetval-v1';
let store=true, saved={};
try{ saved=JSON.parse(localStorage.getItem(K)||'{}'); }
catch(e){ store=false; }

function persist(){
  if(!store) return;
  try{ localStorage.setItem(K,JSON.stringify(saved)); }
  catch(e){ store=false; warn(); }
}
function warn(){
  document.getElementById('warn').textContent =
    'This browser will not save progress for a local file. Your answers are '
    + 'kept for this tab only \\u2014 press "Show results" and copy them out '
    + 'before closing.';
}
function paint(){
  let n=0;
  document.querySelectorAll('.item').forEach(d=>{
    const id=d.id.slice(1);
    const c=d.querySelector('input[name="v'+id+'"]:checked');
    if(saved[id] && !c){
      const r=d.querySelector('input[value="'+saved[id]+'"]');
      if(r) r.checked=true;
    }
    if(d.querySelector('input[name="v'+id+'"]:checked')){
      n++; d.classList.add('done');
      const rv=d.querySelector('.reveal'); if(rv) rv.style.display='block';
    }
  });
  document.getElementById('n').textContent=n;
}
document.addEventListener('change',e=>{
  if(e.target.type!=='radio')return;
  saved[e.target.name.slice(1)]=e.target.value;
  persist();
  paint();
});
// read straight from the radios, so this works even with no storage at all
function dump(){
  let out='id,verdict\\n', n=0;
  document.querySelectorAll('.item').forEach(d=>{
    const id=d.id.slice(1);
    const c=d.querySelector('input[name="v'+id+'"]:checked');
    if(c){ out+=id+','+c.value+'\\n'; n++; }
  });
  document.getElementById('out').innerHTML=
    '<p>'+n+' answered. Select all and copy into validation-verdicts.csv.</p>'
    +'<textarea onclick="this.select()">'+out+'</textarea>';
}
if(!store) warn();
paint();
window.addEventListener('beforeunload',function(e){
  const any=document.querySelector('.item input:checked');
  if(any && !store){ e.preventDefault(); e.returnValue=''; }
});
</script>""")

    OUT_HTML.write_text('\n'.join(parts), encoding='utf-8')


if __name__ == '__main__':
    main()
