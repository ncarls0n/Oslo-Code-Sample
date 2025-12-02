README.txt

!!! Description !!!
cosmo_pert/ is a set of code used to make background and perturbative calculations during inflation. It includes
code to evolve backgrounds, ballistic trajectories, and power spectra.

!!! Contents !!!
The files in this project are organized as follows:
----------------
params.f90           - some widely used parameter declarations
----------------
util.f90             - definitions of utility subroutines for matrix/array handling
----------------
gl_integrator.f90    - genreal purpose Guass-Legendre(GL) integrator
butcher_table.f90    - parameter declarations for the GL integrator
newton_root.f90      - root finding by Newton's method using the GL integrator
gl_macros.h          - macro definitions for the GL integrator and convergence testing
----------------
potential.f90        - definitions of potential functions
eom_bg_cosmic.f90    - background equations of motion(EoM) in cosmic time
eom_pert_cosmic.f90  - perturbative EoM in cosmic time
corr_cosmic.f90      - subroutines to set initial conditions for perturbative EoM and compute spectra from
                       integration variables
macros.h             - macro defintions used by the EoM modules
----------------
evolve_corr.f90      - program to compute and evolve spectra for a range of modes
Makefile_corr
----------------
evolve_pert.f90      - program to evolve the background along with a single pertubative mode
Makefile_pert
----------------
ballistic_ensemble.f90  - program to evolve an ensemble of background solutions
Makefile_ens
---------------
evolve.f90           - program to evolve the background system
Makefile

!!! Future Development !!!
Write and EoM module to inculde the effects of metric perturbations.
