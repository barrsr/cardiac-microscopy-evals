# Source: task A data

**Dataset.** Kobeissi H, Gao X, DePalma SJ, Ewoldt JK, Wang MC, Das SL, Jilberto J, Nordsletten D, Baker BM, Chen CS, Lejeune E. *Strain gauge platforms: Time-lapse microscopy dataset of engineered cardiac microbundles.* Dryad, DOI [10.5061/dryad.sqv9s4nbg](https://doi.org/10.5061/dryad.sqv9s4nbg); mirrored as Zenodo record [11120809](https://zenodo.org/records/11120809), published 2024-05-06.

**License.** CC0 1.0 (public domain dedication), as listed on the Zenodo record.

**Platform reference.** Zhang K, Cloonan PE, Sundaram S, et al., Chen CS. Plakophilin-2 truncating variants impair cardiac contractility by disrupting sarcomere stability and organization. *Science Advances* 7(42):eabh3995 (2021). Analysis software: MicroBundlePillarTrack, https://github.com/HibaKob/MicroBundlePillarTrack.

**File used.** `Type1_condition_025.zip` (654.9 MB, MD5 012aca6b7861e1685fc3f38bc0585d20) -> `condition_025/batch_000/Ex00628.tif`.
Copied unmodified to `data/Ex00628.tif`: 79,444,634 bytes, MD5 a1d8dac26038f1c5c6798d3c5e32b709, 150 frames, 512 x 512, uint16, MetaMorph 7.10.5 stream. Nothing was trimmed, cropped or re-encoded, so every page keeps its MetaMorph acquisition timestamp.

**Documented conditions** (from the deposit's `Type1_experimental_conditions.pdf` and `Type1_common_experimental_metadata.pdf`):
- Condition 025: L333 patient-derived stem cells carrying an HCM-causing MYBPC3 mutation, with ventricular cardiac fibroblasts (Lonza), DMEM high glucose + 10% FBS, pillar stiffness 2.667 N/m, "taken at Day 7 with 1 Hz stimulation".
- Platform: brightfield, Nikon Eclipse Ti, Evolve EMCCD, 4 um/pixel, "Imaging frequency 30 Hz".
- `metadata_type1.csv`: Ex00628 is condition 25, batch 0, comment 0 ("No comment"). The `included_in_plot` field is 0 for all of condition 025.

**Time base found in the file.** The per-plane `acquisition-time-local` stamps run from 15:01:40.036 on 2023-11-28 over 4.527 s for 150 frames (intervals 30 to 31 ms, mean 32.91 frames/s). The header records a 30 ms exposure in frame-transfer mode. The 30 Hz on the metadata sheet does not match the acquisition rate; the 1 Hz pacing confirms the timestamps (see NOTES.md).
