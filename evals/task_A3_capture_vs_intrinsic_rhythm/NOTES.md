# Task A3: 1:1 capture versus an intrinsic rhythm under 1 Hz pacing

The data is nineteen unmodified brightfield TIFF stacks (150 frames each) from condition 023 of the CC0 Kobeissi, Lejeune et al. microbundle deposit (Dryad 10.5061/dryad.sqv9s4nbg, Zenodo 11120809): L550 patient-derived hiPSC-cardiomyocytes carrying an HCM-causing MYBPC3 variant, which the depositors document as recorded at day 7 "with 1 Hz stimulation", in three imaging sessions. The agent decides which tissues follow the stimulus 1:1 and returns `{"n_captured": int, "not_captured": [sorted IDs]}`.

**Ground truth: n_captured = 6 (Ex00560, Ex00561, Ex00567, Ex00569, Ex00570, Ex00572); the other 13 are not captured.** The ground truth rests on the depositors' 1 Hz condition label plus analysis: ground_truth.py times each beat at its contraction upstroke against a relaxed reference frame, on the per-plane MetaMorph timestamps (32.91 frames/s). The six captured tissues have every interval within 0.010 s of 1.000 s and drift of 0.005 s or less. The closest non-captured tissue (Ex00575) has every interval 0.024 to 0.035 s short and drifts 0.085 s ahead of the stimulus grid in three beats; its smallest deviation per interval exceeds the largest captured one by a factor of 2.4, and its largest by 3.5.

**What the task tests.** A tissue whose own rate is slightly above the pacing rate is not captured by the stimulus, yet beats just as regularly as a captured one. Ex00575 (0.976, 0.965, 0.974 s; 61.8 bpm) and Ex00576 (0.960, 0.961, 0.959 s; 62.5 bpm, CV 0.1%) are the trap cases: within 5% of 1 Hz on every interval and very regular, but every beat arrives before the next stimulus could, so neither beat is stimulus-driven. Ex00565 averages 59.1 bpm and ends only 0.048 s off the grid, yet its intervals are 0.844, 0.883 and 1.321 s. Ex00562 and Ex00574 beat regularly at about 0.91 s. Measured on the data: "every interval within 5%" returns 8 captured, a 2% mean-rate window returns 7, a 5% window returns 11, drift alone returns 8, and the displacement spectrum returns 10. Timing beats at the displacement peak instead of the upstroke rejects Ex00570 under a strict window (5). Trusting the sheet's 30 Hz returns 2 (strict) or 12 (10% window).

**Excluded recordings.** Condition 023 has 23 files. Ex00568, Ex00573 and Ex00578 are left out because the capture call flips between the upstroke and peak beat markers, and Ex00577 because its contraction is too faint to time. eval.json, section 2, gives the measured intervals.

**Grader.** `all_of`: `numeric_tolerance` on n_captured (absolute 0) and `label_set_jaccard` on not_captured (pass_threshold 1.0, `answer_field` set to not_captured). Checked against latch-eval-tools 0.4.54 graders. The GT answer passes in any order. Answers keeping Ex00575, keeping Ex00575 and Ex00576, keeping Ex00565, rejecting Ex00570, or giving the right set with a wrong count all fail. The harness grade and the latch-eval-tools grade agree on all twelve runs below.

**Limitations.**
- The margin is 3.5x, narrower than task A2's 8x. The answer needs beat times precise to about a third of a frame, which the upstroke gives (19 captured intervals within 0.010 s of 1.000 s) and the displacement peak does not in the low-amplitude tissues.
- Two other rules also pass, both on upstroke times with the timestamp time base: task A2's criterion (mean within 2%, every interval within 5%, drift under 0.10 s) and a 2% mean-rate window combined with CV under 3%. The task separates on whether regular beating a few percent faster than the stimulus is accepted as capture.
- Ex00572 is the faintest captured tissue. Three failing runs (two GPT-6 Luna, one Haiku) also called it not captured, each alongside the intrinsic-rhythm errors; no run failed on Ex00572 alone. Its upstroke intervals (1.003, 1.001, 0.999 s) and peak drift (+0.007 s) both show it locked to the stimulus.
- Why Ex00575 and Ex00576 are not captured (intrinsic rate above 1 Hz, or a stimulus near threshold in those wells) cannot be told from the images. The answer does not depend on it: neither follows the 1 Hz stimulus.
- Each recording holds 3 to 6 beats.

**Model results** (autonomous agents, three runs each through harness/run_models.py; transcripts in results/trajectories/).

| Model | Passed | Failing runs |
|---|---|---|
| Claude Opus | 3/3 | none |
| Claude Sonnet | 1/3 | r1 kept Ex00556, Ex00575, Ex00576 (n=9); r2 kept Ex00562, Ex00574, Ex00575, Ex00576 (n=10) |
| Claude Haiku | 0/3 | r1 kept Ex00575, Ex00576 (n=8); r2 kept Ex00574, Ex00575, Ex00576 (n=9); r3 kept Ex00556, 562, 563, 564, 574, 575, 576 and rejected Ex00572 (n=12) |
| GPT-6 Luna (Codex CLI) | 0/3 | r1 kept Ex00557, 562, 563, 564, 566, 574, 575, 576 and rejected Ex00572 (n=13); r2 kept Ex00558, 562, 564, 575, 576 and rejected Ex00572 (n=10); r3 kept Ex00558, 562, 564, 575, 576 (n=11) |

Every failing run counted Ex00575 and Ex00576 as captured, which is the intrinsic-rhythm trap. The Claude transcripts state the rules that let them through: intervals of 0.94 to 1.0 s, a dominant period of 0.85 to 1.15 s, or a 1 Hz spectral power fraction; the GPT-6 Luna runs give no criterion beyond the answer block. All three Opus runs rejected Ex00575 and Ex00576 because their intervals sit steadily 3 to 4% below 1.000 s while the captured tissues from the same imaging session sit at 1.000 s. Sonnet r3 passed on the same reasoning ("an intrinsic rate of about 1.04 Hz rather than pacing at 1 Hz"), although it read the frame interval as 30 ms rather than the 30.4 ms in the timestamps.

**Designs considered and dropped.** Two other tasks on this deposit were built and set aside because their answers were not unambiguous. A twitch-force task (pillar deflection times the stated stiffness) was fair only on a hand-picked tissue: on others the two pillars deflect by up to 3x different amounts (Ex00633, 26 and 9 um; Ex00635, 21 and 10 um), and three trackers disagree by about 12%. A paced-versus-spontaneous task could not rely on the condition labels: two tissues labeled spontaneous in condition 028 (Ex00676 and Ex00678) beat at exactly 1.000 s with no drift.
