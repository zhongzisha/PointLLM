

import sys,os,glob,shutil,json
import numpy as np
import pandas as pd
import openslide
import tifffile
from spatialdata_io import xenium

data_dir = '/Users/zhongz2/down/Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE_outs'
data_dir = '/Users/zhongz2/down/Xenium_V1_FFPE_Human_Breast_ILC_outs'

# slide = openslide.open_slide(os.path.join(data_dir, 'morphology.ome.tif'))

cells = pd.read_csv(os.path.join(data_dir, 'cells.csv.gz'), low_memory=False)
cells_pq = pd.read_parquet(os.path.join(data_dir, 'cells.parquet'))


tif = tifffile.TiffFile(os.path.join(data_dir, 'morphology.ome.tif'))
# Access image data and metadata
pages = tif.pages 
for page in pages:
    image_data = page.asarray() 
    metadata = page.tags 

    # Process image data and metadata
    print(image_data.shape)  # Print image dimensions
    print(metadata)  # Access metadata like pixel size, channel names, et

    for k,v in metadata.items():
        print(v.name, v.value)
    break

tif.close()

# 
# ImageWidth 48511
# ImageLength 41405

sdata = xenium(data_dir, cells_as_circles=True)
# sdata.write("main.zarr")

slide = openslide.open_slide(os.path.join(data_dir, 'Xenium_V1_FFPE_Human_Breast_ILC_he_image.ome.tif'))

with open(os.path.join(data_dir, 'experiment.xenium')) as f:
    specs = json.load(f)


# Import Python libraries
# This code uses python v3.12.0, tifffile v2023.9.26, matplotlib v3.8.2
import tifffile
import matplotlib.pyplot as plt

# Option 1: Load full resolution image channels
# The following may produce a warning: 'OME series cannot read multi-file pyramids'. This is because tifffile does not support loading a pyramidal multi-file OME-TIFF file. Only the full resolution (level=0) data will load for all channels in the directory.
fullres_multich_img = tifffile.imread(
"morphology_focus/morphology_focus_0000.ome.tif", is_ome=True, level=0, aszarr=False)

# Examine shape of array (number of channels, height, width), e.g. (4, 40867, 31318)
fullres_multich_img.shape

# Extract number of channels, e.g. 4
n_ch = fullres_multich_img.shape[0]

# Plot each channel
fig, axes = plt.subplots(ncols=n_ch, nrows=1, squeeze=False)
for i in range(n_ch):
    axes[0, i].imshow(fullres_multich_img[i], cmap="gray")
    axes[0, i].set_title(f"Channel: {i}")
plt.savefig('tifffile_fullres_four_channels.png')

# Option 2: Load a single channel image at any resolution, e.g., level=0 or level=1. Note 'is_ome' is set to False.
# Load one of the multi-file OME-TIFF files as a regular TIFF file at full resolution.
fullres_img_tiff = tifffile.imread(
    "morphology_focus/morphology_focus_0000.ome.tif", is_ome=False, level=0)

# Now load the file at downsampled resolution
downsampled_img = tifffile.imread(
    "morphology_focus/morphology_focus_0000.ome.tif", is_ome=False, level=1)
# Examine shape of array (number of channels, height, width)
downsampled_img.shape

# Plot the full resolution and downsampled images side-by-side
fig, axes = plt.subplots(ncols=2, nrows=1, squeeze=False)
axes[0, 0].imshow(fullres_img_tiff, cmap="gray")
axes[0, 0].set_title(f"Full resolution: {fullres_img_tiff.shape}")
axes[0, 1].imshow(downsampled_img, cmap="gray")
axes[0, 1].set_title(f"Downsampled: {downsampled_img.shape}")
plt.savefig('example_fullres_downsample.png')









