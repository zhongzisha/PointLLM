



import sys,os,glob,shutil,json
import numpy as np
import pandas as pd
# import openslide
import tifffile
from spatialdata_io import xenium, visium
import spatialdata as sd
import squidpy as sq
import scanpy as sc
from PIL import Image
# import pyvips
import cv2
from matplotlib import pyplot as plt
# Import Python libraries
# This script was tested with zarr v2.13.6
import zarr

# Function to open a Zarr file
def open_zarr(path: str) -> zarr.Group:
    store = (zarr.ZipStore(path, mode="r")
    if path.endswith(".zip")
    else zarr.DirectoryStore(path)
    )
    return zarr.group(store=store)

data_dir = '/Users/zhongz2/down/Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE_outs'
data_dir = '/Users/zhongz2/down/Xenium_V1_FFPE_Human_Breast_ILC_outs'
data_dir = '/Users/zhongz2/down/Xenium_V1_humanLung_Cancer_FFPE_outs/Xenium_V1_humanLung_Cancer_FFPE_xe_outs'
data_dir = '/Users/zhongz2/down/Xenium_V1_humanLung_Cancer_FFPE_outs'

data_dir = '/Users/zhongz2/down/Visium_FFPE_Human_Breast_Cancer_1.3.0'

data_dir = '.'

# For example, use the above function to open the cells Zarr file, which contains segmentation mask Zarr arrays
root = open_zarr(os.path.join(data_dir, "cell_feature_matrix.zarr.zip"))
root1 = open_zarr(os.path.join(data_dir, "cells.zarr.zip"))


# Look at group array info and structure
root.info
root.tree() # shows structure, array dimensions, data types


# Read in secondary analysis Zarr arrays
root = open_zarr(os.path.join(data_dir, "analysis.zarr.zip"))
# Examples for exploring file contents
# How to show a slice of the clustering_index arrays
root["cell_groups"][0]["indices"][0:9]
# How to show attributes
root["cell_groups"].attrs["group_names"]


# Read in transcripts Zarr arrays
root = open_zarr(os.path.join(data_dir, "transcripts.zarr.zip"))
# Examples for exploring file contents
# How to show array info
root['grids'][0]['0,0']['gene_identity'].shape
root['grids'][0]['0,0']['quality_score'][0:9]
root['grids'][0]['0,0']['location'][0:9,]
# How to show array attributes
root.attrs['major_version']
root['density']['gene'].attrs['gene_names'][0:9]


with tifffile.TiffFile(os.path.join(data_dir, 'morphology.ome.tif')) as tif:
    for i in range(len(tif.pages)):

        image_data = tif.pages[i].asarray()
        image_data[image_data>2048] = 2048
        image_data = (255*(image_data.astype(np.float32)/2048)).astype(np.uint8)

        cv2.imwrite(f'dapi_{i}.jpg', image_data)


    image_data = tif.pages[7].asarray() 
    image_data[image_data>2048] = 2048
    image_data = (255*(image_data.astype(np.float32)/2048)).astype(np.uint8)
    # metadata = page.tags 

    # Process image data and metadata
    # print(image_data.shape)  # Print image dimensions
    # print(metadata)  # Access metadata like pixel size, channel names, et

    # for k,v in metadata.items():
    #     print(v.name, v.value)


sdata = xenium(data_dir,
    n_jobs=8,
    cell_boundaries=True,
    nucleus_boundaries=True,
    morphology_focus=True,
    cells_as_circles=True)
adata = sdata['table']

slide = openslide.open_slide(os.path.join(data_dir, 'Xenium_V1_FFPE_Human_Breast_ILC_he_image.ome.tif'))
tif = tifffile.TiffFile(os.path.join(data_dir, 'Xenium_V1_FFPE_Human_Breast_ILC_he_image.ome.tif'))
r = tif.pages[0].asarray() 
g = tif.pages[1].asarray() 
b = tif.pages[2].asarray() 
tif.close()
he = np.stack([b,g,r], axis=2)



with open(os.path.join(data_dir, 'experiment.xenium')) as f:
    specs = json.load(f)

pixel_size = specs['pixel_size']
scale = 1/pixel_size

patch_size = 224
num_patches = 1000

xy = np.copy(adata.obsm['spatial']) * scale

xmin,ymin = xy.min(axis=0)
xmax,ymax = xy.max(axis=0)

cxs = (xmax-xmin) * np.random.sample(size=num_patches)+xmin+patch_size//2
cys = (ymax-ymin) * np.random.sample(size=num_patches)+ymin+patch_size//2

patche_inds = []
all_coords = []
save_dir = './patches2'
os.makedirs(save_dir, exist_ok=True)
for cx, cy in zip(cxs, cys):
    x1, y1 = cx-patch_size/2, cy-patch_size/2
    x2, y2 = cx+patch_size/2, cy+patch_size/2

    inds = np.where((xy[:, 0] >= x1) & (xy[:, 0] <= x2) & (xy[:, 1] >= y1) & (xy[:, 1] <= y2))[0]
    if len(inds)>0:
        all_coords.append([cx, cy])
        patche_inds.append(inds)

        # img = slide.read_region((int(x1), int(y1)), level=0, size=(patch_size, patch_size)).convert('RGB')
        img = he[int(y1):int(y1+patch_size), int(x1):int(x1+patch_size), :]
        dapi = image_data[int(y1):int(y1+patch_size), int(x1):int(x1+patch_size)]
        dapi = np.stack([dapi, dapi, dapi], axis=2)
        cv2.imwrite(os.path.join(save_dir, f'x{cx}_y{cy}.jpg'), np.concatenate([img, dapi], axis=1))

if False:
    im = np.array(slide.read_region((0,0),0,slide.level_dimensions[0]))[:,:,:3]  # the image looks weird
    im = np.ascontiguousarray(im)
    for x, y in all_coords:
        cv2.rectangle(im, (int(x-patch_size/2), int(y-patch_size/2)), (int(x+patch_size/2), int(y+patch_size/2)), (0, 255, 0), 5)


    img=Image.fromarray(im)

    img_vips = pyvips.Image.new_from_array(img)

    save_filename = 'shown.tif'
    img_vips.tiffsave(save_filename, compression="jpeg",
        tile=True, tile_width=512, tile_height=512,
        pyramid=True,  bigtiff=True)



















