

import sys,os,glob,shutil,json
import numpy as np
import pandas as pd 
import anndata as ad
import scanpy as sc
from matplotlib import pyplot as plt


data_dir = '/data/zhongz2/data/pointllm/'

with open(os.path.join(data_dir, 'PointLLM_brief_description_660K_filtered.json'), 'r') as fp:
    alldata = json.load(fp)

merfish1 = pd.read_csv('/data/zhongz2/download/data/abc_atlas/metadata/MERFISH-C57BL6J-638850-CCF/20230630/views/cell_metadata_with_parcellation_annotation.csv', index_col=0, low_memory=False)

merfish2 = pd.read_csv('/data/zhongz2/download/data/abc_atlas/metadata/MERFISH-C57BL6J-638850-CCF/20230830/views/cell_metadata_with_parcellation_annotation.csv', index_col=0, low_memory=False)

merfish = pd.read_csv('/data/zhongz2/download/data/abc_atlas/metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv', index_col=0, low_memory=False)
adata = ad.read_h5ad('/data/zhongz2/download/data/abc_atlas/expression_matrices/MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-log2.h5ad', backed='r')

patch_size = 100
section_list = merfish['brain_section_label'].unique()
num_patches = 1000



for si, section in enumerate(section_list):
    break

    merfish_sub = merfish[merfish['brain_section_label']==section]
    adata_sub = adata[adata.obs['brain_section_label']==section].to_memory()
    genes = adata[adata.obs['brain_section_label']==section].X.toarray()[:, :500]
    valid_cell_ids = np.where((merfish_sub['x'].notna())&(merfish_sub['y'].notna()))[0]
    merfish_sub = merfish_sub[(merfish_sub['x'].notna())&(merfish_sub['y'].notna())]
    merfish_sub['x'] = merfish_sub['x']*1e3
    merfish_sub['y'] = merfish_sub['y']*1e3

    adata_sub = adata_sub[merfish_sub.index]
    adata_sub.obs = merfish_sub
    genes = genes[valid_cell_ids]
    xy = merfish_sub[['x','y']].values
    adata_sub.obsm['spatial'] = xy
    axes = sc.pl.spatial(adata_sub, color='class_color', spot_size=15, show=False)

    cell_types = merfish_sub['class'].values
    xmin,ymin = xy.min(axis=0)
    xmax,ymax = xy.max(axis=0)

    x1s = (xmax-patch_size-xmin) * np.random.sample(size=num_patches)+xmin
    y1s = (ymax-patch_size-ymin) * np.random.sample(size=num_patches)+ymin

    patches = []
    for x1, y1 in zip(x1s, y1s):
        x2,y2 = x1+patch_size, y1+patch_size

        if np.random.rand()<0.2:
            axes[0].plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], color='k')

        inds = np.where((xy[:, 0] >= x1) & (xy[:, 0] <= x2) & (xy[:, 1] >= y1) & (xy[:, 1] <= y2))[0]
        if len(inds)>0:
            patches.append((xy[inds], genes[inds], cell_types[inds]))
            print(len(inds), len(np.unique(cell_types[inds])))

    plt.savefig(f'{section}.png')
    plt.close('all')

    break
