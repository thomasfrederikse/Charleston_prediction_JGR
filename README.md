# Supplement for 'A hybrid dynamical approach for seasonal prediction of sea-level anomalies: a pilot study for Charleston, South Carolina'
(c) 2022 All Rights Reserved

Authors: Thomas Frederikse [1], Tong Lee [1], Ou Wang [1], Ben Kirtman [2], Emily Becker [2], Ben Hamlington [1], Daniel Limonadi [1], Duane Waliser [1]

[1] NASA Jet Propulsion Laboratory, California Institute of Technology, Pasadena, USA

[2] University of Miami Rosenstiel School of Marine and Atmospheric Science, Miami, FL, USA

Please cite 'A hybrid dynamical approach for seasonal prediction of sea-level anomalies: a pilot study for Charleston, South Carolina' when using these scripts.

Next to these scripts, a data supplement with the adjoint sensitivities, the predictions, and statistics is available on Zenodo: 

This supplement contains the following directories:

# `Scripts`
The `Scripts` directory contains the [Julia](https://www.julialang.org) scripts used to compute the convolution and determine the predictive skill of the various projections. These scripts rely on some external datasets, such as the ECCO forcings [ECCO](https://www.ecco-group.org/) and the CCSM4 forcings part of the North-American multi-model ensemble [(NMME)](https://www.ncei.noaa.gov/products/weather-climate-models/north-american-multi-model). 

## Pre-processing scripts
These scripts perfrom various pre-processing tasks.
- `prepare_obs.jl` Read and prepare the tide-gauge and altimetry observations.
- `regrid_CCSM_nn.jl` Re-grid the CCSM4 forcing data onto the ECCO LLC90 grid using a simple nearest-neighboor approach.
- `regrid_functions.jl` Various functions used by `regrid_CCSM_nn.jl`.
- `read_CCSM_ssh_projections.jl`Read the direct SSH predictions from CCSM4.

## Processing scripts
- `compute_convolution.jl` # Compute the convolution between various forcings and the adjoint sensitivities.
- `create_persistence_projections.jl` # Create the damped persistence predictions.

## Post-processing scripts
- `postprocess_convolution.jl` # Post-process the convolution and other predictions and compute statistics (ACC, RMSE etc.).
- `convert_sens_to_2d.jl` # Re-grid adjoint sensitivities from LLC90 to a 2-dimensional grid for plotting purposes

# `GMT`
This directory contains the figures from the paper and supplement and the GMT scripts and input data to reproduce them. You can download and build GMT from [here](https://www.generic-mapping-tools.org/). All plots have been made using GMT version 6.2. 
