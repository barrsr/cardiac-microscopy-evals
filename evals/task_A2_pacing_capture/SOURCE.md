# Source: task A2 data

**Dataset.** Kobeissi H, Gao X, DePalma SJ, Ewoldt JK, Wang MC, Das SL, Jilberto J, Nordsletten D, Baker BM, Chen CS, Lejeune E. *Strain gauge platforms: Time-lapse microscopy dataset of engineered cardiac microbundles.* Dryad, DOI [10.5061/dryad.sqv9s4nbg](https://doi.org/10.5061/dryad.sqv9s4nbg); mirrored as Zenodo record [11120809](https://zenodo.org/records/11120809), published 2024-05-06.

**License.** CC0 1.0 (public domain dedication), as listed on the Zenodo record.

**Platform reference.** Zhang K, Cloonan PE, Sundaram S, et al., Chen CS. Plakophilin-2 truncating variants impair cardiac contractility by disrupting sarcomere stability and organization. *Science Advances* 7(42):eabh3995 (2021). Analysis software: MicroBundlePillarTrack, https://github.com/HibaKob/MicroBundlePillarTrack.

**Archive used.** `Type1_condition_025.zip` (654,914,176 bytes, MD5 012aca6b7861e1685fc3f38bc0585d20, matching the archive used for task A). All eleven TIFFs in it were copied byte for byte into `data/`, with the two batch subfolders flattened. Nothing was trimmed, cropped or re-encoded, so every page keeps its MetaMorph acquisition timestamp.

| File | Archive path | Bytes | MD5 | First frame acquired |
|------|--------------|-------|-----|----------------------|
| Ex00628.tif | condition_025/batch_000/ | 79,444,634 | a1d8dac26038f1c5c6798d3c5e32b709 | 2023-11-28 15:01:40.036 |
| Ex00629.tif | condition_025/batch_000/ | 79,444,634 | 6297d51e5222bfe3836fe26a2758b697 | 2023-11-28 15:01:01.251 |
| Ex00630.tif | condition_025/batch_000/ | 79,444,934 | 0f7f1525f1ce74bc21bed80c5eba4695 | 2023-11-28 15:00:03.014 |
| Ex00631.tif | condition_025/batch_000/ | 79,444,934 | a9a15a39e2b52ab76d200e0b0bf66fd5 | 2023-11-28 14:59:19.911 |
| Ex00632.tif | condition_025/batch_000/ | 79,444,634 | 08da3e8f86812fcddf5a6122b3be7a5b | 2023-11-28 15:02:17.328 |
| Ex00633.tif | condition_025/batch_000/ | 79,444,934 | f82adafd815314c96b43f6f706597751 | 2023-11-28 15:06:10.770 |
| Ex00634.tif | condition_025/batch_000/ | 79,444,934 | ef6d71faf5b59e3c5d5874d7ed0686e2 | 2023-11-28 15:06:47.665 |
| Ex00635.tif | condition_025/batch_001/ | 79,446,434 | 9be9371ba75065ddc888285048950393 | 2023-11-29 16:24:46.278 |
| Ex00636.tif | condition_025/batch_001/ | 79,446,434 | a4ab40c931291f06577abb4c5f13074f | 2023-11-29 16:28:00.109 |
| Ex00637.tif | condition_025/batch_001/ | 79,446,134 | 51f847c6244c3f8b66f592a2a7ed57a9 | 2023-11-29 16:29:07.596 |
| Ex00638.tif | condition_025/batch_001/ | 79,446,134 | 442af8713d7a78ff4a073f6cf2bfa233 | 2023-11-29 16:29:32.211 |

Every file: 150 frames, 512 x 512, uint16, MetaMorph 7.10.5 stream, "Exposure: 30 msec", sensor mode FT.

**Documented conditions** (deposit files `Type1_experimental_conditions.pdf`, `Type1_common_experimental_metadata.pdf`, `metadata_type1.csv`, `Comments_key.pdf`, all read directly):
- Condition 025: L333 patient-derived stem cells carrying an HCM-causing MYBPC3 mutation, with ventricular cardiac fibroblasts (Lonza), DMEM high glucose + 10% FBS, pillar stiffness 2.667 N/m, "taken at Day 7 with 1 Hz stimulation". Condition 026 is the same line recorded while beating spontaneously.
- Platform: brightfield, Nikon Eclipse Ti, Evolve EMCCD, 4 um/pixel, "Imaging frequency 30 Hz".
- `metadata_type1.csv`: all eleven files carry comment code 0 ("No comment") and `included_in_plot` = 0. The depositors flag nothing about pacing for any of them.

**Time base found in the files.** In every file the per-plane `acquisition-time-local` stamps span 4.527 s for 150 frames (steps of 30 or 31 ms, mean 32.91 frames/s). The 30 Hz on the metadata sheet does not match the acquisition rate; the eight tissues that follow the 1 Hz stimulus sit at 1.000 s per beat on the timestamp time base (see eval.json, section 2).
