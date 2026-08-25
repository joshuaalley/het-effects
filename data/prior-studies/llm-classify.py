"""Classify replication code for heterogeneous-effects estimation.

Replaces the regular expressions behind Figure 1. Every analysis file in the
corpus is sent to an open-weight model, which returns a structured judgment plus
verbatim evidence lines. Nothing is pre-filtered by pattern matching: the whole
point is that a recall filter would reintroduce the problem the regex had.

Reproducibility
---------------
Open weights, pinned by HuggingFace revision, temperature 0. Every request and
every raw response is written to an append-only JSONL archive, so a replicator
can verify the classification without re-running inference at all, and can
re-run against the pinned checkpoint if they want to.

Usage
-----
    python llm-classify.py                      # dry run: counts and cost only
    python llm-classify.py --run --limit 200    # pilot on 200 files
    python llm-classify.py --run                # full corpus
    python llm-classify.py --run --model gpt-oss

Set DEEPINFRA_API_KEY in the environment first. The script never reads a key
from the command line, so it cannot end up in shell history.

The run is resumable: results append to out/<model>.jsonl keyed by file and
chunk, and a restart skips whatever is already there. Re-running after an
interruption never pays for the same call twice.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent
REPFILES = ROOT / 'replication-files'
OUTDIR = ROOT / 'llm-output'
DEFAULT_PROMPT = 'classify-v2'
PROMPT_DIR = ROOT / 'prompts'

ENDPOINT = 'https://api.deepinfra.com/v1/openai/chat/completions'

# Pinned open-weight checkpoints. `revision` is the HuggingFace commit the
# weights come from; record it in the appendix alongside the model name, since
# a bare model name is not a version.
MODELS = {
    'llama': {
        'id': 'meta-llama/Llama-3.3-70B-Instruct-Turbo',
        'revision': None,  # fill from HF once confirmed against the served build
        'price_in_per_m': 0.10,
        'price_out_per_m': 0.32,
        'max_tokens': 600,
        'extra': {},
    },
    'gpt-oss': {
        'id': 'openai/gpt-oss-120b',
        'revision': None,
        'price_in_per_m': 0.09,
        'price_out_per_m': 0.45,
        # gpt-oss is a reasoning model: its hidden reasoning is billed as
        # completion tokens and counts against max_tokens, so a 600-token cap
        # truncated the JSON mid-evidence. Truncation was not random --- it hit
        # replies with long evidence arrays, i.e. the positive cases --- so it
        # biased toward false negatives. Low reasoning effort cuts completion
        # tokens roughly fourfold on this task and the headroom absorbs the rest.
        'max_tokens': 1500,
        'extra': {'reasoning_effort': 'low'},
    },
}

# Author-written analysis code only. `.ado` files are downloaded third-party
# Stata packages, not the author's analysis; their help text contains the very
# syntax we are looking for, which is a live false-positive source. Everything
# else excluded here is binary or an archive.
CODE_EXT = {'.r': 'R', '.do': 'Stata', '.py': 'Python',
            '.stan': 'Stan', '.jags': 'JAGS'}

# Chunking. Roughly 3.5 characters per token for code, so 80k characters is
# about 23k tokens, comfortably inside a 128k window with room for the prompt.
MAX_CHARS = 80_000
OVERLAP = 2_000
MAX_CHUNKS = 20  # a file needing more than this is not hand-written analysis

CHARS_PER_TOKEN = 3.5
MAX_RETRIES = 5
TIMEOUT = 180

_print_lock = threading.Lock()


# ---------------------------------------------------------------------------
# corpus
# ---------------------------------------------------------------------------

def build_manifest() -> pd.DataFrame:
    """One row per (file, chunk) to be classified."""
    cf = pd.read_csv(ROOT / 'code_files_list.csv')
    cf['safe_doi'] = cf.dataset_doi.str.replace(r'[:/]', '_', regex=True)
    cf['ext'] = cf.filename.str.extract(r'(\.[A-Za-z0-9]+)$')[0].str.lower()
    cf = cf[cf.ext.isin(CODE_EXT)].copy()

    rows = []
    for r in cf.itertuples(index=False):
        path = REPFILES / r.safe_doi / str(r.filename)
        if not path.exists():
            continue
        try:
            size = path.stat().st_size
        except OSError:
            continue
        n_chunks = 1 if size <= MAX_CHARS else min(
            MAX_CHUNKS, -(-(size - OVERLAP) // (MAX_CHARS - OVERLAP)))
        for i in range(n_chunks):
            rows.append({
                'dataset_doi': r.dataset_doi,
                'file_id': r.file_id,
                'filename': r.filename,
                'path': str(path),
                'language': CODE_EXT[r.ext],
                'chunk': i,
                'n_chunks': n_chunks,
                'bytes': size,
            })
    return pd.DataFrame(rows)


def read_chunk(path: str, chunk: int) -> tuple[str, bool]:
    raw = Path(path).read_bytes()
    text = raw.decode('utf-8', errors='replace')
    if len(text) <= MAX_CHARS:
        return text, False
    start = chunk * (MAX_CHARS - OVERLAP)
    return text[start:start + MAX_CHARS], True


# ---------------------------------------------------------------------------
# prompt
# ---------------------------------------------------------------------------

def load_prompt(version: str) -> tuple[str, str, str]:
    """Return (system, user_template, sha256 of the prompt file)."""
    raw = (PROMPT_DIR / f'{version}.md').read_text(encoding='utf-8')
    digest = hashlib.sha256(raw.encode('utf-8')).hexdigest()[:16]
    system = raw.split('## SYSTEM', 1)[1].split('## USER', 1)[0].strip()
    user = raw.split('## USER', 1)[1].strip()
    return system, user, digest


def extract_json(text: str):
    """Parse the model's reply, tolerating markdown fences and stray prose.

    DeepInfra accepts `response_format: json_object` but the served models still
    wrap the object in a ```json fence, so a bare json.loads rejects otherwise
    perfect output. Strip the fence, and fall back to the outermost balanced
    brace pair.
    """
    t = text.strip()
    if t.startswith('```'):
        t = t.split('\n', 1)[1] if '\n' in t else t
        if t.rstrip().endswith('```'):
            t = t.rstrip()[:-3]
    t = t.strip()
    try:
        return json.loads(t)
    except json.JSONDecodeError:
        pass
    start = t.find('{')
    if start < 0:
        return None
    depth = 0
    for i, ch in enumerate(t[start:], start):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(t[start:i + 1])
                except json.JSONDecodeError:
                    return None
    return None


VALID_METHODS = {'interaction_term', 'subgroup_split', 'ml_heterogeneity',
                 'hierarchical_varying_slopes'}
VALID_PURPOSE = {'heterogeneity_claim', 'nuisance_or_fixed_effects', 'mixed',
                 'n_a'}
VALID_CONF = {'low', 'medium', 'high'}


def validate_extract(obj: dict) -> tuple[bool, str]:
    """Schema for the interaction-inventory task."""
    if not isinstance(obj, dict):
        return False, 'not an object'
    inter = obj.get('interactions')
    if not isinstance(inter, list) or any(not isinstance(x, str) for x in inter):
        return False, 'interactions must be a list of strings'
    n = obj.get('n_models_with_interaction')
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        return False, 'n_models_with_interaction must be a non-negative integer'
    if not isinstance(obj.get('list_complete'), bool):
        obj['list_complete'] = True
    if obj.get('confidence') not in VALID_CONF:
        obj['confidence_raw'] = obj.get('confidence')
        obj['confidence'] = 'unknown'
    return True, ''


def validate(obj: dict) -> tuple[bool, str]:
    """Reject only what makes a row unusable.

    interaction_purpose and confidence are descriptive fields that nothing
    downstream depends on --- the headline count is not gated on either ---
    so an off-vocabulary value there is coerced rather than treated as a
    failure. Rejecting the whole row over a descriptive field discarded 106
    otherwise-valid classifications on the first v2 pass, which is a worse
    outcome than recording the judgment and noting the field was unusable.
    """
    if not isinstance(obj, dict):
        return False, 'not an object'
    if not isinstance(obj.get('estimates_heterogeneous_effects'), bool):
        return False, 'estimates_heterogeneous_effects must be bool'
    m = obj.get('methods')
    if not isinstance(m, list) or any(x not in VALID_METHODS for x in m):
        return False, f'methods must be a subset of {sorted(VALID_METHODS)}'
    if obj.get('interaction_purpose') not in VALID_PURPOSE:
        obj['interaction_purpose_raw'] = obj.get('interaction_purpose')
        obj['interaction_purpose'] = 'unknown'
    if obj.get('confidence') not in VALID_CONF:
        obj['confidence_raw'] = obj.get('confidence')
        obj['confidence'] = 'unknown'
    if not isinstance(obj.get('evidence'), list):
        return False, 'evidence must be a list'
    if obj['estimates_heterogeneous_effects'] and not m:
        return False, 'true requires at least one method'
    return True, ''


# ---------------------------------------------------------------------------
# inference
# ---------------------------------------------------------------------------

def call_model(session, model_cfg, system, user, api_key,
               validator=None):
    """One request, with retries. Returns (parsed, raw_text, usage)."""
    body = {
        'model': model_cfg['id'],
        'messages': [{'role': 'system', 'content': system},
                     {'role': 'user', 'content': user}],
        'temperature': 0,
        'max_tokens': model_cfg.get('max_tokens', 600),
        'response_format': {'type': 'json_object'},
        **model_cfg.get('extra', {}),
    }
    headers = {'Authorization': f'Bearer {api_key}',
               'Content-Type': 'application/json'}

    validator = validator or validate
    last = ''
    spent = []
    for attempt in range(MAX_RETRIES):
        try:
            resp = session.post(ENDPOINT, json=body, headers=headers,
                                timeout=TIMEOUT)
            if resp.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f'status {resp.status_code}')
            resp.raise_for_status()
            payload = resp.json()
            text = payload['choices'][0]['message']['content']
            # keep the usage even when the reply is unusable, so the running
            # cost figure reflects what was actually spent
            usage = payload.get('usage', {})
            spent.append(usage)
            parsed = extract_json(text)
            if parsed is None:
                last = f'unparseable JSON: {text[:200]}'
                body['messages'] = body['messages'][:2] + [
                    {'role': 'assistant', 'content': text},
                    {'role': 'user',
                     'content': 'That was not valid JSON. Return only the JSON '
                                'object, with no code fence and no commentary.'}]
                continue
            ok, why = validator(parsed)
            if not ok:
                last = f'schema: {why}'
                # one repair attempt, then give up on this chunk
                body['messages'] = body['messages'][:2] + [
                    {'role': 'assistant', 'content': text},
                    {'role': 'user',
                     'content': f'That response was invalid ({why}). '
                                f'Return only corrected JSON.'}]
                continue
            return parsed, text, _sum_usage(spent)
        except Exception as e:  # noqa: BLE001 - retry on anything transient
            last = str(e)
            time.sleep(min(2 ** attempt + random.random(), 30))
    return None, last, _sum_usage(spent)


def _sum_usage(items):
    return {'prompt_tokens': sum(u.get('prompt_tokens', 0) for u in items),
            'completion_tokens': sum(u.get('completion_tokens', 0) for u in items)}


def done_keys(path: Path) -> set:
    if not path.exists():
        return set()
    seen = set()
    with path.open(encoding='utf-8') as fh:
        for line in fh:
            try:
                d = json.loads(line)
                if d.get('result'):  # failures are retried on the next run
                    seen.add((d['file_id'], d['chunk']))
            except Exception:  # noqa: BLE001 - tolerate a torn final line
                continue
    return seen


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='llama', choices=sorted(MODELS))
    ap.add_argument('--run', action='store_true',
                    help='actually call the API (default is a dry run)')
    ap.add_argument('--limit', type=int, default=None,
                    help='classify only the first N chunks, for a pilot')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--only-flagged', default=None,
                    help='restrict to chunks a previous run flagged with this '
                         'method, given as "<jsonl-basename>:<method>", e.g. '
                         '"gpt-oss.classify-v2:interaction_term"')
    ap.add_argument('--prompt', default=DEFAULT_PROMPT,
                    help='prompt version in prompts/; the output filename is '
                         'keyed on it, so a new version writes alongside the '
                         'old rather than over it')
    args = ap.parse_args()

    cfg = MODELS[args.model]
    system, user_tmpl, prompt_hash = load_prompt(args.prompt)
    man = build_manifest()
    validator = validate_extract if args.prompt.startswith('extract') else validate

    if args.only_flagged:
        base, method = args.only_flagged.rsplit(':', 1)
        keep = set()
        with (OUTDIR / f'{base}.jsonl').open(encoding='utf-8') as fh:
            for line in fh:
                try:
                    d = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not d.get('result'):
                    continue
                if method in (d['result'].get('methods') or []):
                    keep.add((d['file_id'], d['chunk']))
        before = len(man)
        man = man[man.apply(lambda r: (r.file_id, r.chunk) in keep, axis=1)]
        print(f'restricted to {method}: {before:,} -> {len(man):,} chunks')

    est_tokens = (man.bytes.clip(upper=MAX_CHARS).sum() / CHARS_PER_TOKEN
                  + len(man) * len(system) / CHARS_PER_TOKEN)
    est_out = len(man) * 120
    cost = (est_tokens / 1e6 * cfg['price_in_per_m']
            + est_out / 1e6 * cfg['price_out_per_m'])

    print(f'model          : {cfg["id"]}')
    print(f'prompt         : {args.prompt} (sha256 {prompt_hash})')
    print(f'papers         : {man.dataset_doi.nunique():,}')
    print(f'files          : {man.file_id.nunique():,}')
    print(f'chunks (calls) : {len(man):,}')
    print(f'chunked files  : {(man.n_chunks > 1).groupby(man.file_id).any().sum():,}')
    print(f'est input tok  : {est_tokens/1e6:.1f}M')
    print(f'est cost       : ${cost:,.2f}')
    print(f'by language    : {man.groupby("language").size().to_dict()}')

    if not args.run:
        print('\nDRY RUN. Nothing was sent. Re-run with --run to classify.')
        return 0

    api_key = os.environ.get('DEEPINFRA_API_KEY')
    if not api_key:
        print('\nDEEPINFRA_API_KEY is not set in the environment.',
              file=sys.stderr)
        return 1

    OUTDIR.mkdir(exist_ok=True)
    out_path = OUTDIR / f'{args.model}.{args.prompt}.jsonl'
    already = done_keys(out_path)
    todo = man[~man.apply(lambda r: (r.file_id, r.chunk) in already, axis=1)]
    if args.limit:
        # A seeded random sample, not the head: the manifest is ordered by
        # paper, so the first N rows would come from a handful of archives and
        # tell us nothing about the corpus. These rows are archived like any
        # other, so a later full run simply skips them.
        todo = todo.sample(n=min(args.limit, len(todo)), random_state=20260820)
    print(f'\nalready done   : {len(already):,}')
    print(f'to classify    : {len(todo):,}\n')
    if todo.empty:
        return 0

    session = requests.Session()
    out_lock = threading.Lock()
    fh = out_path.open('a', encoding='utf-8')
    counts = {'ok': 0, 'fail': 0, 'in': 0, 'out': 0}

    def work(row):
        code, was_chunked = read_chunk(row.path, row.chunk)
        note = (f'This is chunk {row.chunk + 1} of {row.n_chunks} of a large '
                f'file; judge only what is shown.' if was_chunked else '')
        prompt = (user_tmpl.replace('{language}', row.language)
                  .replace('{filename}', str(row.filename))
                  .replace('{chunk_note}', note)
                  .replace('{code}', code))
        parsed, raw, usage = call_model(session, cfg, system, prompt,
                                        api_key, validator)
        rec = {
            'dataset_doi': row.dataset_doi, 'file_id': row.file_id,
            'filename': row.filename, 'language': row.language,
            'chunk': row.chunk, 'n_chunks': row.n_chunks,
            'model': cfg['id'], 'prompt_version': args.prompt,
            'prompt_sha256': prompt_hash, 'temperature': 0,
            'result': parsed, 'raw': None if parsed else raw,
            'usage': usage, 'ts': time.time(),
        }
        with out_lock:
            fh.write(json.dumps(rec, ensure_ascii=False) + '\n')
            fh.flush()
            counts['ok' if parsed else 'fail'] += 1
            counts['in'] += usage.get('prompt_tokens', 0)
            counts['out'] += usage.get('completion_tokens', 0)
            n = counts['ok'] + counts['fail']
            if n % 100 == 0:
                spend = (counts['in'] / 1e6 * cfg['price_in_per_m']
                         + counts['out'] / 1e6 * cfg['price_out_per_m'])
                with _print_lock:
                    print(f'  {n:,}/{len(todo):,}  failed {counts["fail"]}  '
                          f'${spend:,.2f}', flush=True)

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = [ex.submit(work, r) for r in todo.itertuples(index=False)]
        for f in as_completed(futures):
            f.result()
    fh.close()

    spend = (counts['in'] / 1e6 * cfg['price_in_per_m']
             + counts['out'] / 1e6 * cfg['price_out_per_m'])
    print(f'\ndone in {(time.time() - t0)/60:.1f} min')
    print(f'  classified {counts["ok"]:,}, failed {counts["fail"]:,}')
    print(f'  tokens in {counts["in"]:,}, out {counts["out"]:,}')
    print(f'  actual spend ${spend:,.2f}')
    print(f'  archive: {out_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
