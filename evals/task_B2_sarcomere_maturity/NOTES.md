# Task B2 notes

The data are eight day-32 hiPSC-cardiomyocytes from the Allen Institute ACTN2-mEGFP dataset (Gerbin et al. 2021), all from one imaging session at 0.1195257 um per pixel. Each cell is the authors' maximum projection cropped with an 8 um margin, plus its mask. Three cells are striated: both experts scored them 4 or 5, and the classifier puts 0.87 to 0.93 of their area in organized z-disks. The other five have puncta and fibers with no z-disks: experts scored them 1 or 2, and the z-disk fraction is 0.02 or less. No score-3 cell was used. The agent has to name the striated cells and report their median sarcomere length.

Ground truth is cell_B, cell_D and cell_E, n = 3, median 1.73 um. Each cell's value is the mean of the two v1 methods: B 1.731, D 1.695, E 1.773. The image agrees with the author call for every cell. Spectral peak prominence is 8.3 to 13.1 in the striated cells and 0.8 to 1.4 in the rest, whose radial spectra have no peak at all.

The traps that fire are about running the measurement without asking whether the cell has sarcomeres:
- Aligned profiles return 1.96 to 2.35 um for the non-striated cells, and the all-cell median is 2.04.
- The authors' GLCM peak distance returns 1.32 to 1.68 um for them, and the all-cell median is 1.56.
- On unmasked crops, the 2D FFT picks up striated neighbors in cell_A and cell_H. A range test then calls five cells striated.

Grading needs the exact set, n_striated = 3, and a median within 1.73 +/- 0.12 um.

**Model results** (autonomous agents, three runs each; transcripts in results/trajectories/). Claude Haiku, Sonnet and Opus and GPT-6 Luna passed every run (12 of 12); every run masked the cells, fit a 2D spectrum and dropped the five cells with no sarcomere peak. As built, the panel does not separate these models. A harder variant would add expert-called score-3 cells, at the cost of a less clean ground truth.

The Allen terms permit noncommercial redistribution only.
