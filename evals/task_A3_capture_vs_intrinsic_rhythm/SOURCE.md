# Source: task A3 data

**Dataset.** Kobeissi H, Gao X, DePalma SJ, Ewoldt JK, Wang MC, Das SL, Jilberto J, Nordsletten D, Baker BM, Chen CS, Lejeune E. *Strain gauge platforms: Time-lapse microscopy dataset of engineered cardiac microbundles.* Dryad, DOI [10.5061/dryad.sqv9s4nbg](https://doi.org/10.5061/dryad.sqv9s4nbg); mirrored as Zenodo record [11120809](https://zenodo.org/records/11120809), published 2024-05-06.

**License.** CC0 1.0 (public domain dedication), as listed on the Zenodo record.

**Platform reference.** Zhang K, Cloonan PE, Sundaram S, et al., Chen CS. Plakophilin-2 truncating variants impair cardiac contractility by disrupting sarcomere stability and organization. *Science Advances* 7(42):eabh3995 (2021). Analysis software: MicroBundlePillarTrack, https://github.com/HibaKob/MicroBundlePillarTrack.

**Archive used.** `Type1_condition_023.zip` (1,258,653,935 bytes, MD5 1dc07af0f95876443ef555d8968a6e12, matching the MD5 published on the Zenodo record). Nineteen of its twenty-three TIFFs were copied byte for byte into `data/`, with the three batch subfolders flattened. Nothing was trimmed, cropped or re-encoded, so every page keeps its MetaMorph acquisition timestamp.

| File | Archive path | Bytes | MD5 | First frame acquired |
|------|--------------|-------|-----|----------------------|
| Ex00556.tif | condition_023/batch_000/ | 79,444,634 | f852aba118a6e3635ed722c78e267aba | 2023-11-28 15:21:31.730 |
| Ex00557.tif | condition_023/batch_000/ | 79,444,334 | eba14f60b2e1199fe329cf90cb00e9bd | 2023-11-28 15:22:28.156 |
| Ex00558.tif | condition_023/batch_000/ | 79,444,634 | 1534e95e0591974b8ffa9db9aaf70bb5 | 2023-11-28 15:23:09.904 |
| Ex00559.tif | condition_023/batch_000/ | 79,444,634 | 39c5dccfcfd0024e502acf234a23a088 | 2023-11-28 15:23:48.017 |
| Ex00560.tif | condition_023/batch_001/ | 79,446,134 | 60d56436f70cfbc3e25652cdaa535b45 | 2023-11-29 16:36:58.282 |
| Ex00561.tif | condition_023/batch_001/ | 79,446,434 | 7a1b7af2f2c3ab9a2d716d63dbc5e08b | 2023-11-29 16:39:45.915 |
| Ex00562.tif | condition_023/batch_001/ | 79,446,134 | 2c63c39958af34fe9dd233d656609c3c | 2023-11-29 16:45:06.051 |
| Ex00563.tif | condition_023/batch_001/ | 79,446,134 | ecc30082b78797a7195438feec775d5c | 2023-11-29 16:51:00.325 |
| Ex00564.tif | condition_023/batch_001/ | 79,446,134 | f9114e1401a774da5ca889a5a422dc8d | 2023-11-29 16:53:56.721 |
| Ex00565.tif | condition_023/batch_001/ | 79,446,134 | b21ee8c8f7c7b2e6d9c0e130ee08084e | 2023-11-29 16:54:32.838 |
| Ex00566.tif | condition_023/batch_002/ | 79,444,634 | 54f191740d22346676f9b4faf6051b52 | 2023-12-12 14:20:26.449 |
| Ex00567.tif | condition_023/batch_002/ | 79,445,234 | 7a6231ed2e58f531856e737be24af90c | 2023-12-12 14:52:53.982 |
| Ex00569.tif | condition_023/batch_002/ | 79,445,234 | f33becfea170485932f670972963374d | 2023-12-12 14:57:49.085 |
| Ex00570.tif | condition_023/batch_002/ | 79,444,634 | a8322af5513e517a39b240d7b99e602e | 2023-12-12 14:59:39.095 |
| Ex00571.tif | condition_023/batch_002/ | 79,444,634 | 58978aa85ada6ff5f478a6dff1dbf74b | 2023-12-12 14:22:00.199 |
| Ex00572.tif | condition_023/batch_002/ | 79,444,934 | 9bd82c9ebaff21e54f7241e00254e505 | 2023-12-12 14:22:41.367 |
| Ex00574.tif | condition_023/batch_002/ | 79,445,234 | 0a680ffc9ef8a5e8dead40efccbd4c35 | 2023-12-12 14:32:25.718 |
| Ex00575.tif | condition_023/batch_002/ | 79,444,634 | ed71d33e8c0e1b95e3f270d0a11611fe | 2023-12-12 14:39:12.964 |
| Ex00576.tif | condition_023/batch_002/ | 79,444,634 | 46c5bc14dcae12305459d5721e174591 | 2023-12-12 14:42:10.927 |

Every file: 150 frames, 512 x 512, uint16, MetaMorph stream, "Exposure: 30 msec", sensor mode FT.

**Left out.** Ex00568, Ex00573, Ex00577 and Ex00578 (all condition_023/batch_002/) are not in the snapshot. In Ex00568, Ex00573 and Ex00578 the capture call depends on whether beats are timed at the contraction upstroke or at the displacement peak; in Ex00577 the contraction is too faint to time. eval.json, section 2, gives the measured values.

**Documented condition** (deposit files `Type1_experimental_conditions.pdf`, `Type1_common_experimental_metadata.pdf`, `metadata_type1.csv`, `Comments_key.pdf`, all read directly). Condition 23 in `Type1_experimental_conditions.pdf` reads:

> L550 patient-derived stem cells harboring an HCM causing mutation in MYBPC3 and vCFs (Lonza) cultured in DMEM high glucose + 10% FBS, 1% P/S, 1% glutamax, 1% NEAA supplemented with .033mg/ml aprotinin and 150ug/mL L-ascorbic acid, pillar stiffness of 2.667 N/m, taken at Day 7 with 1 Hz stimulation

Condition 24 is the same line "taken at Day 7 while beating spontaneously".
- Platform (`Type1_common_experimental_metadata.pdf`): brightfield, Nikon Eclipse Ti, Evolve EMCCD, 4 μm/pixel, imaging frequency 30 Hz.
- `metadata_type1.csv`: all nineteen files carry `included_in_plot` = 0. Ex00562 carries comment code 4 ("Tracked outputs contain drift"); the other eighteen carry code 0 ("No comment"). The depositors flag nothing about pacing for any of them.

**Time base found in the files.** In every file the per-plane `acquisition-time-local` stamps span 4.526 to 4.528 s for 150 frames (steps of 30 or 31 ms, 32.91 to 32.92 frames/s). The 30 Hz on the metadata sheet does not match the acquisition rate; the six tissues that follow the 1 Hz stimulus sit at 1.000 s per beat on the timestamp time base (see eval.json, section 2).
