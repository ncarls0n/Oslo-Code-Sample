#!/bin/sh
gfortran -O2 \
    ../src/modules/External/intreal_types.f90 \
    ../src/modules/GlobalVariables/cosmoparams.f90 \
    ../src/modules/globalvariables/params.f90 \
    ../src/modules/globalVariables/input_parameters.f90 \
    ../src/cosmology/Dlin_params.f90 \
    ../src/modules/Solvers/Solvers.f90 \
    ../src/modules/HomogeneousEllipsoid/homeosubs.f90 \
    ../src/cosmology/psubs_Dlinear.f90 \
    run_hom_ellipse_tab.f90 \
    -o run_hom_ellipse_tab_MAC.x
exit
