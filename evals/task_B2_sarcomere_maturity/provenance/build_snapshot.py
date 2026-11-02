"""Build the Task B2 data snapshot from the Allen Institute public bucket.

Selects eight day-32 cells from one imaging session (20190816_M02_001, plate 5500000014),
crops each from the authors' 2D field (channel '488' = max projection of ACTN2-mEGFP)
with an 8 um margin around the cell's bounding box, and writes an OME-TIFF image and a
target-cell mask per cell. Pixel values are unchanged.

usage: python build_snapshot.py [cache_dir]
"""
import json, sys, urllib.parse, urllib.request
from pathlib import Path
import numpy as np, pandas as pd, tifffile

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'data'
CACHE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / '_cache'
BUCKET = 'https://allencell.s3.amazonaws.com/aics/integrated_transcriptomics_structural_organization_hipsc_cm/'
PX = 0.1195257          # um/px, SlideBook calibration tag 65007 and TIFF XResolution of the M02_001 raw files
MARGIN_UM = 8.0

# (panel id, FOV name, label in the authors' cell segmentation). Ids from numpy default_rng(20261102).permutation(8).
PANEL = [
    ('cell_A', '20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 13', 1),
    ('cell_B', '20190816_M02_001_Capture 3 - Capture 3 - B2 [144] Montage - 17', 4),
    ('cell_C', '20190816_M02_001_Capture 3 - Capture 3 - B5 [144] Montage - 2', 5),
    ('cell_D', '20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 14', 5),
    ('cell_E', '20190816_M02_001_Capture 3 - Capture 3 - B2 [144] Montage - 11', 5),
    ('cell_F', '20190816_M02_001_Capture 3 - Capture 3 - B5 [144] Montage - 9', 5),
    ('cell_G', '20190816_M02_001_Capture 3 - Capture 3 - B5 [144] Montage - 15', 4),
    ('cell_H', '20190816_M02_001_Capture 3 - Capture 3 - B4 [144] Montage - 27', 4),
]


def fetch(key):
    dst = CACHE / Path(key).name
    if not dst.exists() or dst.stat().st_size == 0:
        CACHE.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(BUCKET + urllib.parse.quote(key), dst)
    return dst


def main():
    seg = pd.read_csv(fetch('2d_segmented_fields/metadata.csv'))
    seg['fov'] = seg.original_fov_location.str.split('/').str[-1]
    tab = pd.read_csv(fetch('revised_manuscript_plots/data.csv'), low_memory=False)
    tab = tab[tab.Dataset == 'OldFish'].copy()
    tab['fov'] = tab.fov_path.str.split('/').str[-1]
    tab['napariCell_ObjectNumber'] = pd.to_numeric(tab.napariCell_ObjectNumber)
    OUT.mkdir(exist_ok=True)
    rows = []
    for cid, fov, lab in PANEL:
        s = seg[seg.fov == fov].iloc[0]
        a = tifffile.imread(fetch('2d_segmented_fields/' + s['2D_fov_tiff_path']))
        assert a.shape[0] == 10
        img, cell = a[1], a[9]          # channel_defs.json: '488' = 1, 'cell' = 9
        m = cell == lab
        ys, xs = np.where(m)
        M = int(round(MARGIN_UM / PX))
        H, W = m.shape
        y0, y1, x0, x1 = max(ys.min() - M, 0), min(ys.max() + M + 1, H), max(xs.min() - M, 0), min(xs.max() + M + 1, W)
        crop, mk = img[y0:y1, x0:x1], m[y0:y1, x0:x1].astype(np.uint8)
        md = dict(axes='YX', PhysicalSizeX=PX, PhysicalSizeXUnit='µm', PhysicalSizeY=PX, PhysicalSizeYUnit='µm')
        tifffile.imwrite(OUT / f'{cid}_actn2.ome.tif', crop, ome=True, compression='zlib',
                         metadata=dict(md, Name=f'{cid} ACTN2-mEGFP maximum projection'))
        tifffile.imwrite(OUT / f'{cid}_mask.ome.tif', mk, ome=True, compression='zlib',
                         metadata=dict(md, Name=f'{cid} mask (1 = target cell, 0 = other)'))
        t = tab[(tab.fov == fov) & (tab.napariCell_ObjectNumber == lab)].iloc[0]
        rows.append(dict(cell_id=cid, fov_id=int(s.fov_id), fov=fov, label=lab, crop_rows=f'{y0}:{y1}', crop_cols=f'{x0}:{x1}',
                         mask_px=int(mk.sum()), field_2d=s['2D_fov_tiff_path'], age_days=int(t['Cell age']),
                         expert_mh=t.mh_structure_org_score, expert_kg=t.kg_structure_org_score,
                         frac_organized_zdisks=t['Fraction cell area organized z-disks'],
                         frac_organized_puncta=t['Fraction cell area organized puncta'],
                         frac_disorganized_puncta=t['Fraction cell area disorganized puncta'],
                         frac_fibers=t['Fraction cell area fibers'],
                         frac_diffuse=t['Fraction cell area diffuse/other'],
                         frac_background=t['Fraction cell area background'],
                         author_peak_distance_um=t['Peak distance (μm)'], cell_area_um2=t['Cell area (μm^2)']))
        print(cid, fov, lab, crop.shape, int(mk.sum()))
    pd.DataFrame(rows).to_csv(HERE / 'panel_author_annotations.csv', index=False)


if __name__ == '__main__':
    main()
