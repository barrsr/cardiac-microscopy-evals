"""Run each eval against coding agents from three model makers and grade the answers.

Each run gets a fresh working directory holding only the task's data/ snapshot,
the task prompt and a Python environment. The agent answers in the eval's
<EVAL_ANSWER> block, which is graded against the eval's own grader config.
Full transcripts are written to results/trajectories/<eval id>/<model>/r<run>/.

Agents:
  claude:<model>  Claude Code (`claude -p`), isolated from the operator's own
                  CLAUDE.md, memory, skills, MCP servers and user settings
  codex:[model]   OpenAI Codex CLI (`codex exec`), workspace-write sandbox
  gemini:<model>  Gemini CLI (`gemini -p`), macOS seatbelt sandbox; needs GEMINI_API_KEY

usage:
  python harness/run_models.py --evals A2,B2 --models claude:haiku,claude:sonnet,claude:opus,codex: --runs 3
Fetch the data first: python scripts/fetch_all.py
"""
import argparse, concurrent.futures as cf, json, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = os.environ.get('CME_PYTHON', sys.executable)  # interpreter the agents are told to use
EVALS = {'A': ROOT / 'evals' / 'task_A_beat_rate', 'B': ROOT / 'evals' / 'task_B_sarcomere_length',
         'A2': ROOT / 'evals' / 'task_A2_pacing_capture', 'B2': ROOT / 'evals' / 'task_B2_sarcomere_maturity'}
OUT = ROOT / 'results'

PREAMBLE = (
    'You are a computational biologist working in the current directory, which '
    'contains the data for this task under ./data. A Python environment with numpy, '
    'scipy, scikit-image, opencv, tifffile, pandas and matplotlib is available as '
    f'{PY}; ffmpeg and ffprobe are on the PATH. Work only inside this directory. '
    'Inspect the data, do the analysis, and finish with the answer block exactly as specified.\n\n')


def grade(ev, text):
    """Grade the last answer block with the eval's numeric_tolerance / label_set_jaccard / all_of config."""
    m = re.findall(r'<EVAL_ANSWER>\s*(\{.*?\})\s*</EVAL_ANSWER>', text or '', re.S)
    if not m:
        return None, False, 'no answer block'
    try:
        ans = json.loads(m[-1])
    except Exception as e:
        return None, False, f'unparseable: {e}'
    g = ev['grader']
    ok = True
    for ch in (g['config']['children'] if g['type'] == 'all_of' else [g]):
        cfg = ch['config']
        if ch['type'] == 'label_set_jaccard':
            f = cfg.get('answer_field', 'cell_types_predicted')
            if f not in ans or not isinstance(ans[f], list):
                return ans, False, f'missing {f}'
            a, b = set(map(str, ans[f])), set(map(str, cfg['ground_truth_labels']))
            ok &= len(a & b) / len(a | b) >= cfg.get('scoring', {}).get('pass_threshold', 0.9)
            continue
        for k, gt in cfg['ground_truth'].items():
            tol = cfg['tolerances'][k]
            if k not in ans or not isinstance(ans[k], (int, float)):
                return ans, False, f'missing {k}'
            lim = tol['value'] * (abs(gt) if tol['type'] == 'relative' else 1)
            ok &= abs(ans[k] - gt) <= lim
    return ans, ok, ''


def command(vendor, name, prompt, work):
    if vendor == 'codex':
        return ['codex', 'exec', '--skip-git-repo-check', '--sandbox', 'workspace-write', '--json',
                '-o', str(work / '.last_message.txt'), '-C', str(work)] + (['-m', name] if name else []) + [prompt]
    if vendor == 'gemini':
        return ['gemini', '-p', prompt, '--approval-mode', 'yolo', '--sandbox', '--skip-trust',
                '--output-format', 'json'] + (['-m', name] if name else [])
    return ['claude', '-p', prompt, '--model', name, '--output-format', 'stream-json', '--verbose',
            '--setting-sources', 'project', '--strict-mcp-config', '--disable-slash-commands',
            '--settings', json.dumps({'claudeMdExcludes': [str(Path.home() / '.claude' / '**')],
                                      'autoMemoryEnabled': False}),
            '--allowedTools', f'Read Glob Grep Write Bash({PY}:*) Bash(ls:*) Bash(ffprobe:*) Bash(ffmpeg:*) Bash(file:*)']


def final_message(vendor, stdout, work):
    if vendor == 'codex':
        p = work / '.last_message.txt'
        return p.read_text() if p.exists() else stdout, None
    if vendor == 'gemini':
        try:
            return json.loads(stdout[stdout.find('{'):]).get('response', ''), None
        except Exception:
            return stdout, None
    final, cost = '', None
    for line in stdout.splitlines():
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get('type') == 'result':
            final, cost = d.get('result', ''), d.get('total_cost_usd')
    return final, cost


def one(key, model, run):
    ev = json.loads((EVALS[key] / 'eval.json').read_text())
    vendor, _, name = model.partition(':')
    slug = model.replace(':', '_')
    work = Path(tempfile.mkdtemp(prefix=f'cme_{key}_{slug}_{run}_'))
    shutil.copytree(EVALS[key] / 'data', work / 'data')
    tdir = OUT / 'trajectories' / ev['id'] / slug / f'r{run}'
    tdir.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, 'MPLBACKEND': 'Agg'}  # no GUI windows from headless agents
    t0 = time.time()
    try:
        p = subprocess.run(command(vendor, name, PREAMBLE + ev['task'], work), cwd=work,
                           capture_output=True, text=True, timeout=1500, env=env)
        stdout, stderr = p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        stdout = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or '')
        stderr = 'timeout'
    (tdir / 'trajectory.jsonl').write_text(stdout + ('' if vendor == 'claude' else '\n--- stderr ---\n' + stderr))
    final, cost = final_message(vendor, stdout, work)
    ans, ok, why = grade(ev, final)
    rec = {'eval': ev['id'], 'model': model, 'run': run, 'answer': ans, 'pass': bool(ok),
           'note': why, 'seconds': round(time.time() - t0), 'cost_usd': cost}
    (tdir / 'result.json').write_text(json.dumps(rec, indent=1))
    shutil.rmtree(work, ignore_errors=True)
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--evals', default='A2,B2')
    ap.add_argument('--models', default='claude:haiku,claude:sonnet,claude:opus,codex:')
    ap.add_argument('--runs', type=int, default=3)
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args()
    jobs = [(k, m, r) for k in a.evals.split(',') for m in a.models.split(',') for r in range(1, a.runs + 1)]
    with cf.ThreadPoolExecutor(a.jobs) as pool:
        for rec in pool.map(lambda j: one(*j), jobs):
            print(json.dumps(rec), flush=True)
