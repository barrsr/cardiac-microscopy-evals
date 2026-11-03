# Offbeat

**Cardiac microscopy evals.** Benchmark tasks for AI agents that analyze microscopy data, built from public cardiac datasets. Each task has one right answer that is fixed by the experiment itself or by expert annotation, and a trap built from how paced cardiac tissue and immature cardiomyocytes behave. Every task was run as an autonomous agent on Anthropic and OpenAI models, and every run's transcript is in this repository.

**Results page with run-by-run transcripts: https://barrsr.github.io/cardiac-microscopy-evals/**

The tasks follow the [SpatialBench](https://github.com/latchbio/spatialbench) eval format: a task prompt that ends in a fixed answer block, a grader config, and notes that give the ground truth, the traps, a separation table of approaches that were actually run, and the tolerance rationale.

## Tasks

The two tasks are A2 and A3. A1, B1 and B2 are earlier versions, kept as the development record.

| Task | Data | Question | Answer |
|---|---|---|---|
| **A2** pacing capture | 11 brightfield time-lapses of hiPSC-cardiomyocyte microbundles under 1 Hz field stimulation | Which tissues show 1:1 capture of the stimulus? | 8 captured; Ex00636, Ex00637 and Ex00638 are not |
| **A3** capture vs own rhythm | 19 brightfield time-lapses of hiPSC-cardiomyocyte microbundles under 1 Hz field stimulation, three imaging sessions | Which tissues show 1:1 capture of the stimulus? | 6 captured (Ex00560, Ex00561, Ex00567, Ex00569, Ex00570, Ex00572); 13 are not |
| B2 striation panel (set aside) | 8 ACTN2-mEGFP confocal projections of day-32 hiPSC-cardiomyocytes, with masks | Which cells have striated myofibrils, and what is their median sarcomere length? | cells B, D and E; 1.73 µm |
| A1 beat rate (first version) | one tissue from A2's set | Beat rate in bpm | 60.0 bpm |
| B1 sarcomere length (first version) | one ACTN2-mEGFP confocal z-stack | Sarcomere length in µm | 1.72 µm |

Both tasks ask which paced tissues are captured, and both traps come from what capture means. In A2, three tissues are not captured, and two of them still average close to 60 bpm, so an analysis that accepts anything near 1 Hz counts them. A2 can be passed by accepting any tissue whose mean rate is within 2 to 3% of 60 bpm on the file timestamps. On A3 that rule returns 7 or 8 captured tissues instead of 6, because two of its tissues beat as regularly as several captured ones but a few percent faster than the stimulus, which means they are beating on their own rhythm. A1 and B1 were the first versions, with traps in the signal processing, such as two motion peaks per beat, a stated 30 Hz frame rate that disagrees with the 32.9 fps file timestamps, and measuring along the image axis instead of the myofibril. B2 moved the sarcomere trap into the biology, but every model passed it, so it was set aside.

## Results

Each cell counts passing runs out of completed runs. Every run was an independent agent session with only the task data and a Python environment.

| Model | A2 capture | A3 capture vs own rhythm | A1 (v1) | B1 (v1) | B2 (set aside) |
|---|---|---|---|---|---|
| Claude Haiku | 3/3 | 0/3 | 3/3 | 3/3 | 3/3 |
| Claude Sonnet | 2/3 | 1/3 | 3/3 | 3/3 | 3/3 |
| Claude Opus | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 |
| GPT-6 Luna (Codex CLI) | 1/3 | 0/3 | 0/3 | 3/3 | 3/3 |

On A3, every failing run counted Ex00575 and Ex00576, which beat steadily at 0.96 to 0.98 s per beat, as captured by the 1 Hz stimulus; only Claude Opus held every interval to the 1.000 s period in all three runs. On A2, Claude Sonnet's miss and one of GPT-6 Luna's two misses counted the tissues beating at 57 and 62 bpm as captured; GPT's other miss read the sheet's 30 Hz and called every tissue captured. On A1, every GPT-6 Luna run answered 55 to 56 bpm for a tissue paced at 60, close to what the depositors' 30 Hz frame rate gives against the file's 32.9 fps timestamps. B1 and B2 were passed by every model in every run. Claude ran through Claude Code, isolated from any personal configuration, and GPT through the Codex CLI in its workspace sandbox; a short fixed environment note plus the task prompt was the only input. Gemini 3 Flash was set up through the Gemini CLI, but the free API tier's daily quota did not allow a full set of runs, so it is not reported.

## Repository layout

```
evals/<task>/eval.json        task prompt, grader config, ground truth, traps and separation table
evals/<task>/ground_truth.py  reproduces the ground truth and every naive answer in the notes
evals/<task>/fetch_data.py    downloads the source data and rebuilds the exact snapshot, with checksums
evals/<task>/SOURCE.md        citation, license, and what was selected and why
evals/<task>/NOTES.md         summary of the task and its open questions
evals/<task>/figures/         figures used on the results page
evals/task_A3_*/make_video.py builds the Figure 1 clip from the ground-truth beat times
evals/task_B*/provenance/     source metadata and the B2 snapshot builder
harness/run_models.py         runs the evals as agents (Claude Code, Codex CLI, Gemini CLI) and grades them
scripts/fetch_all.py          fetches and verifies the data for all five tasks
results/trajectories/         one folder per eval, model and run: result.json and the full transcript
site/                         builds docs/index.html, the results page
```

## Reproduce

```bash
python -m venv .venv && . .venv/bin/activate
pip install numpy scipy scikit-image opencv-python-headless tifffile matplotlib pandas requests
python scripts/fetch_all.py                         # about 2.4 GB; checksums verified
python evals/task_A2_pacing_capture/ground_truth.py # prints the per-tissue table and the separation table
python harness/run_models.py --evals A2,A3 --models claude:sonnet,codex:gpt-6-luna --runs 3
python site/build_page.py                           # rebuilds docs/index.html from results/
```

The harness needs the `claude` and `codex` command-line tools signed in (and the `gemini` CLI with `GEMINI_API_KEY` for Gemini runs). Claude runs are isolated from the operator's own Claude Code setup (no personal CLAUDE.md, memory, skills, MCP servers or user settings).

## Data and licenses

- **A1, A2, A3:** Kobeissi, Lejeune et al., "Strain gauge platforms: Time-lapse microscopy dataset of engineered cardiac microbundles," Dryad [10.5061/dryad.sqv9s4nbg](https://doi.org/10.5061/dryad.sqv9s4nbg), CC0. The platform is described in Zhang et al., Science Advances 2021. Files are used unmodified.
- **B1, B2:** Allen Institute for Cell Science, Gerbin et al., Cell Systems 2021, [10.1016/j.cels.2021.05.001](https://doi.org/10.1016/j.cels.2021.05.001). The Allen Institute's [terms of use](https://www.allencell.org/terms-of-use.html) allow noncommercial research use with attribution, so this repository ships scripts that rebuild the snapshots from the public bucket instead of the images themselves. Use in a commercial benchmark would need the Allen Institute's written permission.

Code in this repository is released under the MIT License.

Ryan Barrs, PhD
