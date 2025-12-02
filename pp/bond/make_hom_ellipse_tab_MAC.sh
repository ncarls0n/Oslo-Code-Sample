#!/bin/sh
gfortran -O2 intreal_types.f90 homeosubs.f90 psubs_Dlinear.f90 hom_sphere_tab.f90 run_hom_ellipse_tab.f90 -o run_hom_ellipse_tab_MAC.x
exit
