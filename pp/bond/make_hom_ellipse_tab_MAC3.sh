#!/bin/sh
gfortran -O2 intreal_types.f90 params.f90 input_parameters.f90 Dlin_params.f90 homeosubs.f90 psubs_Dlinear.f90 run_hom_sphere_tab.f90 cosmoparams.f90 Solvers.f90 run_hom_ellipse_tab.f90 -o run_hom_ellipse_tab_MAC.x
exit
