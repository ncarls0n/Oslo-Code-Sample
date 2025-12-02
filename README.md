# Oslo-Code-Sample
This repository contains a sample of some of the codes I have worked on recently during my PhD research at CITA. Please see the readme for a breakdown of what is in this repository.

The codes in this repository are representative of a few proprietary codes that we have been using at CITA. I explain my contributions to each code in the following sections. We plan to make public versions of some of these codes in the next 6 months, but this work has not been completed yet, so I've done my best to explain the code structure of each subdirectory.

## Peak Patch

In the directory `Oslo-Code-Sample/pp/` contains a subset of the [Peak Patch](https://mocks.cita.utoronto.ca/index.php/Peak_Patch_and_WebSky#The_Peak_Patch_simulations) dark matter halo code. This Fortran-based code was first introduced by [Bond and Meyers in 1996](https://ui.adsabs.harvard.edu/abs/1996ApJS..103....1B/abstract) and updated for massively parallel computing by [Stein et al., 2018](https://arxiv.org/abs/1810.07727). 

My primary contributions to the Peak Patch simulations have been numerous coding updates to clean up the code and improve efficiency, as well as developing the novel primordial non-Gaussianity (PNG) models. Code efficiency updates are a bit challenging to isolate, as I've made many changes all through the source code `Oslo-Code-Sample/pp/src/` over the last few years, but there is a reasonable correlation that anywhere there are comments, I have made changes because most of the comments in this code base were added by me in an effort to imporve readability. The PNG models I have added can be found in `Oslo-Code-Sample/pp/src/modules/RandomField/RandomField.f90`, lines 290-525 outline each of the PNG models. These then reference several functions and subroutines that I added such as `zetang_sli`, `T_interp`, `chi_to_zeta`. Peak Patch is highly modular, so any updates ultimately result in cascading changes. In particular, PNG model updates have required additional updates to `Oslo-Code-Sample/pp/src/modules/RandomField/gaussian_field.f90`, `Oslo-Code-Sample/pp/src/modules/RandomField/pktable.f90` and many other scripts.

## Peak Patch post-processing pipeline in python

I have also produced a large Python class for post-processing Peak Patch and WebSky data, `Oslo-Code-Sample/peakpatchtools/peakpatchtools.py` which was made wholly by me. With this code, Peak Patch simulation results can be loaded into python object of class `PeakPatch`, from here many manipulations can be done. I produced a Jupyter notebook for an undergraduate summer student `Oslo-Code-Sample/JupyterNotebooks/example.ipynb` which provides a nice overview of how `peakpatchtools.py` works. 

## Fortran code to measure the powerspectrum from a large discretized field

Another cleaner example of a code wholly written by me is a script to measure the power spectrum from a large discretized field in parallel on our HPC clusters. 
