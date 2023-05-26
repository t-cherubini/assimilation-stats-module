# Validation routines
This repository contains the software (python) used to analyze WRF output data for a JGR publication. 
The goal of this software is to compare WRF output from different assimilation runs: 1) Assimilation of conventional observation only (m-conv); 2) Assimilation of conventional and Transformed Retrievals (TR, m-oper); 2) Assimilation of conventional, TR and Microwave Radiances (TRMW, m-foper)'
The WRF ouput were first interpolated on common grids to either GFS or ECMWF analyses data and on fixed vertical pressure levels by means of the UPP software (https://dtcenter.org/community-code/unified-post-processor-upp-wrf), then accessed with the python scripts in this repository. The netcdf files are read using the Xarray libraries (https://docs.xarray.dev/en/stable/). Available routines have been manipulated and changed in time and not cleaned or well documented. Please, use at your own risks. Always double check your results.

https://zenodo.org/badge/latestdoi/631101105
