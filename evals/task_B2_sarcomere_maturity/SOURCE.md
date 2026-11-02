# Source

**Paper.** Gerbin KA, Grancharova T, Donovan-Maiye RM, et al. Cell states beyond transcriptomics: Integrating structural organization and gene expression in hiPSC-derived cardiomyocytes. *Cell Systems* 12(6):670-687.e10 (2021). doi:10.1016/j.cels.2021.05.001. Preprint: bioRxiv 10.1101/2020.05.26.081083.

**Data package.** Allen Institute for Cell Science, AWS Open Data, `s3://allencell/aics/integrated_transcriptomics_structural_organization_hipsc_cm/` (anonymous access).

**License.** Allen Institute for Cell Science terms of use (https://www.allencell.org/terms-of-use.html). Use, copying, redistribution and derivative works are permitted for research and noncommercial purposes with attribution. Commercial redistribution requires written permission from the Allen Institute. No Creative Commons license applies.

## Files used

| Role | Object in the bucket |
|---|---|
| 2D fields (max projections + segmentation) | `2d_segmented_fields/2D_fov_tiff_path/<hash>_5500000014_63X_20190816_S1_P<n>_<well>_annotations_corrected.tiff`, channel 1 (`488`, ACTN2-mEGFP) and channel 9 (`cell`); channel map in `2d_segmented_fields/supporting_files/channel_defs.json` |
| FOV index | `2d_segmented_fields/metadata.csv` |
| Expert scores, classifier fractions, peak distance | `revised_manuscript_plots/data.csv` (columns `mh_structure_org_score`, `kg_structure_org_score`, `Fraction cell area ...`, `Peak distance (μm)`) |
| Pixel calibration check | `raw_images/FISH/5500000014/20190816_M02_001_Capture 3 - Capture 3 - B2 [144] Montage - 11_C5.tif` and `... B4 [144] Montage - 23_C5.tif` (TIFF tags read by HTTP range request) |

The exact FOV, label and crop window of each panel cell are in `provenance/panel_author_annotations.csv`. `provenance/build_snapshot.py` downloads the inputs and rebuilds `data/`.

## What was selected and why

Eight day-32 cells, all from one imaging session (20190816_M02_001, plate 5500000014, wells B2, B4 and B5), so the pixel calibration is identical and age cannot hint at the answer.

- Three striated cells: both experts scored 4 or 5 and the classifier puts 0.87 to 0.93 of the area in organized z-disks.
- Five cells without z-disks: both experts scored 1 or 2 and the z-disk fraction is 0.00 to 0.02. Their ACTN2 is organized puncta, fibers and disorganized puncta (z-bodies on premyofibril-like fibers and scattered puncta). They were chosen to cover scores 1/1, 1/2 and 2/2, and to avoid cells cut by the field edge. Several sit next to striated neighbors, so their crops contain real sarcomeres outside the mask (enough to move the unmasked FFT onto a sarcomere peak for cell_A and cell_H).
- No score-3 cell was used, because score 3 is the mixed state where the call is a judgment.

The score scale as summarized by work built on this dataset (Le et al., SarcNet, arXiv 2405.17926): 1 sparse, disorganized puncta; 2 denser, more organized puncta; 3 puncta mixed with fibers and z-discs; 4 regular but not aligned z-discs; 5 almost organized and aligned z-discs.

## Pixel size

The session 20190816_M02_001 raw files carry 0.1195257 um/px in the SlideBook calibration tag (65007) and the TIFF XResolution, the value used in Task B v1. The source OME-XML and the paper round it to 0.120 um. The snapshot carries 0.1195257 um.

A second session from the same day, 20190816_M02_002, carries 0.1243068 um/px in the same tags. It was excluded from the panel so that every file has the same calibration.

## Projection

The image for each cell is the authors' own maximum projection (channel `488` of the 2D field). In Task B v1 the maximum projection of the raw C5 stack for FOV B4 Montage-23 of this session was pixel-for-pixel identical to that channel.
