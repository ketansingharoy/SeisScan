# SeisScan
Scan coherent seismic signals in continuous data.

# Ceate a conda environment and install dependencies
conda create -n seisscan
conda activate seisscan
conda install -c conda-forge obspy dask jupyterlab pandas cartopy pyproj utm

# Getting started
### Prepare database
p1_prepare_db.ipynb

### Map of the study area
p2_map.ipynb

### Record section of waveform
p3_plot_repord_section.ipynb

### Compute local similarity

### Prepare Travel time lookup table
p5_traveltime_lookup_table.ipynb

### Perform backprojection
p6_backprojection.ipynb

