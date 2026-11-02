# Source

**Paper.** Gerbin KA, Grancharova T, Donovan-Maiye RM, et al. Cell states beyond transcriptomics: Integrating structural organization and gene expression in hiPSC-derived cardiomyocytes. *Cell Systems* 12(6):670-687.e10 (2021). doi:10.1016/j.cels.2021.05.001. Preprint: bioRxiv 10.1101/2020.05.26.081083.

**Data package.** Allen Institute for Cell Science, AWS Open Data, `s3://allencell/aics/integrated_transcriptomics_structural_organization_hipsc_cm/` (anonymous access).

**License.** Allen Institute for Cell Science terms of use (https://www.allencell.org/terms-of-use.html). Use, copying, redistribution and derivative works are permitted for research and noncommercial purposes with attribution. Commercial redistribution requires written permission from the Allen Institute. No Creative Commons license applies.

## Files used

| Role | Object in the bucket |
|---|---|
| Raw 3D ACTN2-mEGFP stack | `raw_images/FISH/5500000014/20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 23_C5.tif` (and its `.tif.xml`) |
| Cell segmentation | `2d_segmented_fields/2D_fov_tiff_path/3dc13480_5500000014_63X_20190816_S1_P23_B4_annotations_corrected.tiff`, channel 9 (`cell`) |
| Author per-cell features | `revised_manuscript_plots/data.csv`; `automated_local_and_global_structure/metadata.csv` (CellId `fov-323-cell-3`) |

Channel C5 is the 488 nm / 525 nm channel (`488/TL 50um Dual`). Its maximum projection is pixel-for-pixel identical to channel 1 (`488`) of the authors' 2D field, which confirms it is ACTN2-mEGFP.

## What was selected and why

One cell: plate 5500000014, well B4, field Montage-23, cell 3, day 32 of differentiation, fixed (RNA-FISH plate). Two expert annotators scored it 5 and 4 (mean 4.5 on the 1-5 organization scale), and the authors' classifier labels 77% of its area as organized z-disks. It is a large cell (7,227 um2) with long, regular myofibrils that run about 44 degrees from the image x axis, which is the condition under which axis-aligned measurements fail.

The snapshot crop is rows 250:1250 and columns 600:1630 of the 1736 x 1776 field, all 50 z planes, pixel values unchanged. The mask is the authors' label 3 set to 1; it holds 501,844 pixels, the same as the authors' `Total_Area` for this cell.

## Pixel size

- Source OME-XML (`.tif.xml` and embedded): PhysicalSizeX/Y = 0.120 um, PhysicalSizeZ = 0.330 um.
- Same file, SlideBook private tag 65007 and TIFF XResolution/YResolution: 0.1195257 um per pixel.
- Paper methods: 3i spinning disk, 63x/1.2 NA W C-Apochromat, 0.83x tube lens (52.29x), CSU-W1, Prime BSI sCMOS, "pixel size 0.120 um in X-Y and 0.33 um in Z".
- The authors' per-cell peak distances are multiples of 0.12 um, consistent with that calibration.

The snapshot carries the calibrated 0.1195257 um; the rounded 0.120 um changes the answer by 0.4%. The original SlideBook XML is kept in `provenance/original_slidebook_ome_metadata.xml`.

## Published reference for this cell

`revised_manuscript_plots/data.csv` gives Peak distance = 1.80 um (15 px at 0.12 um, GLCM correlation method of Sutcliffe et al. 2018, highest-intensity slice) and the paper describes "1.5 - 2 um mEGFP-alpha-actinin-2 band spacing ... for most cells".
