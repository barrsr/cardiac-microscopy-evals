"""Assemble the results page from the graded runs on disk.

Reads every results/trajectories/<eval>/<model>/r*/result.json, condenses each
transcript to its commands and final message, embeds the task figures, and writes
docs/index.html, which GitHub Pages serves.
"""
import base64, html, json, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TRAJ = ROOT / 'results' / 'trajectories'

TASKS = {  # eval id -> folder, short key
    'brightfield_contractility_microbundle_pacing_capture_v1': ('task_A2_pacing_capture', 'A2'),
    'spinningdisk_hipsc_cm_actn2_striation_panel_sarcomere_length_v2': ('task_B2_sarcomere_maturity', 'B2'),
    'brightfield_contractility_microbundle_beat_rate_v1': ('task_A_beat_rate', 'A1'),
    'spinningdisk_hipsc_cm_actn2_sarcomere_length_v1': ('task_B_sarcomere_length', 'B1'),
}
MODELS = [  # folder name, label, maker
    ('claude_haiku', 'Claude Haiku', 'Anthropic'), ('claude_sonnet', 'Claude Sonnet', 'Anthropic'),
    ('claude_opus', 'Claude Opus', 'Anthropic'), ('codex_gpt-6-luna', 'GPT-6 Luna', 'OpenAI'),
]


SCRUB = [(re.compile(r'/Users/[^/\s]+/\.claude/jobs/[^/\s]+/tmp/lbenv/bin/python3?'), 'python'),
         (re.compile(r'(~|/Users/[^/\s"\']+)?/\.claude/jobs/[^/\s]+/tmp/lbenv'), '<env>'),
         (re.compile(r'/private/tmp/claude-\d+/[^\s"\']+'), '<tmp>'),
         (re.compile(r'~/\.claude/jobs/[^/\s"\']+/tmp'), '<scratch>'),
         (re.compile(r'(/private)?/var/folders/[^\s"\']*?/T/(lb|cme)_[A-Za-z0-9_-]+'), '<workdir>'),
         (re.compile(r'(/private)?/var/folders/[^\s"\']+'), '<tmp>'),
         (re.compile(r'/Users/[^/\s"\']+'), '~')]


def scrub(s):
    for rx, rep in SCRUB:
        s = rx.sub(rep, s)
    return s


def clip(s, n):
    s = scrub(s or '')
    return s if len(s) <= n else s[:n] + ' …'


def condense(path, folder):
    """Return a list of (kind, text) steps from one transcript."""
    raw = path.read_text(errors='ignore') if path.exists() else ''
    steps = []
    if folder.startswith('codex'):
        for line in raw.splitlines():
            try:
                d = json.loads(line)
            except Exception:
                continue
            it = d.get('item') or {}
            if d.get('type') != 'item.completed':
                continue
            if it.get('type') == 'command_execution':
                cmd = re.sub(r"^/bin/zsh -lc ", '', it.get('command', ''))
                steps.append(('cmd', clip(cmd, 1400)))
                if it.get('aggregated_output'):
                    steps.append(('out', clip(it['aggregated_output'], 600)))
            elif it.get('type') == 'agent_message':
                steps.append(('say', clip(it.get('text', ''), 1600)))
    elif folder.startswith('gemini'):
        try:
            body = raw.split('\n--- stderr ---')[0]
            d = json.loads(body[body.find('{'):])
            steps.append(('say', clip(d.get('response', ''), 2400)))
        except Exception:
            steps.append(('say', clip(raw, 1200)))
    else:
        for line in raw.splitlines():
            try:
                d = json.loads(line)
            except Exception:
                continue
            if d.get('type') == 'assistant':
                for c in d['message']['content']:
                    if c.get('type') == 'tool_use':
                        inp = c.get('input', {})
                        steps.append(('cmd', clip(inp.get('command') or json.dumps(inp), 1400)))
                    elif c.get('type') == 'text' and c['text'].strip():
                        steps.append(('say', clip(c['text'], 1600)))
            elif d.get('type') == 'user':
                for c in d['message'].get('content', []) if isinstance(d['message'].get('content'), list) else []:
                    if c.get('type') == 'tool_result':
                        t = c.get('content')
                        t = t if isinstance(t, str) else ' '.join(x.get('text', '') for x in t or [] if isinstance(x, dict))
                        steps.append(('out', clip(t, 600)))
    return steps


def runs():
    out = []
    for ev, (folder, key) in TASKS.items():
        for mdir, label, maker in MODELS:
            for rj in sorted((TRAJ / ev / mdir).glob('r*/result.json')):
                r = json.loads(rj.read_text())
                tj = (rj.parent / 'trajectory.jsonl').read_text(errors='ignore')
                if 'exhausted your daily quota' in tj or ('status: 503' in tj and not r.get('answer')):
                    continue  # provider quota or outage, not a model answer
                out.append({'task': key, 'model': label, 'maker': maker, 'run': r['run'],
                            'answer': r.get('answer'), 'pass': r.get('pass'), 'seconds': r.get('seconds'),
                            'steps': condense(rj.parent / 'trajectory.jsonl', mdir)})
    return out


def img(rel):
    p = ROOT / rel
    return 'data:image/png;base64,' + base64.b64encode(p.read_bytes()).decode()


if __name__ == '__main__':
    data = runs()
    figs = {
        'A2': img('evals/task_A2_pacing_capture/figures/fig1_phase_vs_pacing_grid.png'),
        'B2': img('evals/task_B2_sarcomere_maturity/figures/fig1_striated_vs_immature.png'),
        'A1': img('evals/task_A_beat_rate/figures/fig1_motion_trace.png'),
        'B1': img('evals/task_B_sarcomere_length/figures/fig1_profile_axis_vs_myofibril.png'),
    }
    tpl = (HERE / 'template.html').read_text()
    page = (tpl.replace('/*__RUNS__*/[]', json.dumps(data))
               .replace('/*__FIGS__*/{}', json.dumps(figs)))
    head = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            '<meta name="description" content="Microscopy benchmark tasks with biology traps, built from public cardiac data, '
            'run as agents on Claude, GPT and Gemini, with every transcript.">\n'
            '<style>body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style>\n')
    i = page.index('<div class="wrap">')  # title, fonts and styles belong in <head>
    (ROOT / 'docs' / 'index.html').write_text(head + page[:i] + '</head>\n<body>\n' + page[i:] + '\n</body>\n</html>\n')
    print(len(data), 'runs;', round(len(page) / 1e6, 2), 'MB')
