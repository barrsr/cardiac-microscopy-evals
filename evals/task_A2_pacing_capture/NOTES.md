# Task A2: which paced microbundles show 1:1 capture

The data is all eleven unmodified brightfield TIFF stacks (Ex00628 to Ex00638, 150 frames each) from condition 025 of the CC0 Chen/Lejeune microbundle deposit (Dryad 10.5061/dryad.sqv9s4nbg, Zenodo 11120809). The depositors document the condition as recorded at day 7 under 1 Hz electrical stimulation, and the prompt says so. The agent decides which tissues follow the stimulus 1:1 and returns `{"n_captured": int, "not_captured": [sorted IDs]}`.

**Ground truth: n_captured = 8; not captured = Ex00636, Ex00637, Ex00638.** ground_truth.py builds one event per beat from displacement against a relaxed reference frame and times each beat with the per-plane MetaMorph timestamps. The eight captured tissues have all 27 of their intervals within 0.010 s of 1.000 s and drift of 0.012 s or less from a 1.000 s grid. The three others each have at least one interval 0.079 s or more off and drift of 0.137 s or more:
- Ex00636 averages 57.3 bpm; its intervals lengthen beat by beat (1.007, 1.056, 1.079 s) and it slides 0.14 s off the grid.
- Ex00637 averages 62.1 bpm and is irregular (0.914, 1.122, 0.892, 0.935 s). Intervals shorter than 1 s cannot be stimulus-driven.
- Ex00638 beats regularly at 50.1 bpm.

**What the task tests.** Capture means beats at the stimulus period, interval by interval, not a rate near 1 Hz. The 30 Hz frame rate on the metadata sheet is wrong (32.91 frames/s from the timestamps), and frame differencing gives two peaks per beat. Measured on the data: a 5% or 10% rate window returns 10 captured; a "rate at or above 1 Hz" rule returns 9; rhythm regularity alone returns 9; trusting 30 Hz returns 0 (strict) or 9 (10% window); the frame-difference spectrum peaks at a harmonic in four captured tissues, and the displacement spectrum puts Ex00637 at exactly 1.00 Hz. 2:1 block and alternans are not present in this condition and are not claimed.

**Grader.** `all_of`: `numeric_tolerance` on n_captured (absolute 0) and `label_set_jaccard` on not_captured (pass_threshold 1.0, `answer_field` set to not_captured). Checked against latch-eval-tools 0.4.54 graders: the GT answer passes, and n=10, n=9 and a two-ID list fail.

**Limitations.**
- A mean-rate window of 2% or 3% on the correct time base passes. With three or four intervals per recording, the drifting and irregular tissues already land 3.5 to 4.5% from 60 bpm. The task therefore separates on the time base and the width of the capture window more than on rhythm analysis.
- Three of the four batch_001 tissues (recorded 2023-11-29) do not capture; Ex00635 from the same session does. Whether Ex00636 to Ex00638 were under-stimulated, sat in a device whose electrodes were not driving, or failed to capture at threshold cannot be told from the images. The answer does not depend on it: none of the three follows the 1 Hz stimulus.
- Each recording holds 4 to 5 beats, which is enough to separate the classes here with an eightfold margin but too short to characterize rhythm further.

**Model results** (autonomous agents, three runs each; transcripts in results/trajectories/). Claude Opus and Haiku passed 3 of 3 each and Claude Sonnet 2 of 3; GPT-6 Luna passed 1 of 3. Sonnet's miss and one GPT miss answered n_captured = 10, counting Ex00636 (57.3 bpm) and Ex00637 (62.1 bpm) as captured, which is the near-rate trap in the separation table; the other GPT miss called all eleven tissues captured.
