These codes generate the power spectrum and transfer functions for a field $\chi$ transverse to the inflaton during inflation which contributes to a generic multi-field inflation model. The program that computes spectra is in evolve_pert.f90, the program that determines the parametric form of the transfer function $\Delta\Phi = T[\chi]$ is in  

The basic idea is to set initial spectra for the fields subhorizon where they follow approximately the Minkowski space correlations, then evolve them forward k-by-k along with some cosmological background equations until they have the desired k/(aH). The equations of motion are stored in eom_bg_cosmic.f90 and eom_pert_cosmic.f90. The details of the evolution are kind of involved, if you are interested I'll dig up my notes on them.

The makefile is Makefile_corr, it is set for gfortran and you will need to modify the library paths.

The form of potential is in potential.f90, if you what to match a particular \eplsion(H) you will need to do that here.

Everything is calculated in dimensionless quantities where:
\bar\phi = \phi/M_{PL},
\bar x = x/\mu,
\bar t = t/\mu,
\bar V = V/(\mu^2*M_{Pl}^2),
and likewise with derivatives. The mass \mu is defined as M_{Pl}/mpl, (mpl is a parameter in params.f90).

The integrator is in gl_integrator.f90, gl_macros.h, and butcher_table.f90.

Please let me know if you need any help getting this up and running.




Procedure:

Before and after the instability, the potential V as a function of the inflaton $\phi$ is monotonically decreasing, and in the transverse field, $\chi$, the potential is just a well (like a quadratic or similar). During the instability, trajectories separate, and $V(\chi)$ looks like a Mexican hat potential, while $V(\phi)$ continues its slow decrease.

The parameters we want to control for are the "strength" of the instability (something akin to the effective $f_NL$ in the toy model) and the k range. The strength of the instability can be roughly controlled with the following parameters in `potential.f90`:

    real(dl), parameter, dimension(nfld) :: m2 = (/1._dl,1._dl/) ! m_phi^2,m_chi^2
    real(dl), parameter :: phi_p      = 8.525_dl
    real(dl), parameter :: phi_w      = 0.1_dl
    real(dl), parameter :: lambda_chi = 1.6e5
    real(dl), parameter :: vev        = 0.1_dl

The paraameter `m2` is an array encoding the squared masses of the $\phi$ and $\chi$ fields. This is in units of $\mu$, which is set to $M_{Pl}/10^5$ in params.f90 (where $M_{Pl}$ is the reduced planck mass). Heavier fields will see stronger k dependence at low k.

The parameter `phi_p` is the value of the inflaton $\phi$ at the middle of the instability, and `phi_w` is the width (in $\phi$) of the instability, so these together determine when the instability occurs and how longlived it is in $\phi$. These (as with all other field parameters in this code) are in units of $M_{Pl}^2$.

The parameter `lambda_chi` is related to the maximum depth of the transverse minima in $V(\phi,\chi)$. That is, for the value of $\phi$ and $\chi$ at which the transverse minima are deepest, $\phi=\phi_m$ and $\chi=\pm\chi_m$, `lambda_Chi` is related to $V(\phi_m,0)-|V(\phi_m,\pm\chi_m)|$. 

The parameter `vev` is the vacuum expectation value, which is the width in $\chi$ of the transverse instability $|V(\phi_m,\pm\chi_m)-V(\phi_m,0)|$.

The Fourier-space k modes are set in `evolve_corr.f90` by parameters:

    integer, parameter :: nstep = 2**20     ! maximum number of steps
    integer, parameter :: nk = 10           ! number of k modes
    real(dl), parameter :: dk = 1._dl       ! fundamental k mode
    real(dl) :: k2                          ! square wavenumber

The parameter `nstep` is the maximum number of steps that the integrator takes to find a convergent value. The resultant $\chi$ power spectrum will be an array of length `nk` correpsonding to wavenumbers `sqrt(k2)` separated by `dk`. This code simply works in Fourier modes so these wavenumbers aren't the spatial Fourier wavenumbers $k$ with units of 1/Mpc, however, since we don't know exactly when during inflation these instabilities would take place (or even precisely how long inflation is) the relation between `sqrt(k2)` and $k$ is a free parameter, meaning e.g. we could set it so that the entire `sqrt(k2)` range is in the large scale structure band.

To see the instability, we want dk << H ~ 3.434 (in planck units) and nk > (\lambda_\chi v^2)^.5/dk 

Once these parameters are set, load the intel fortran compiler module on any CITA computer and run the module:

    module load intel
    make clean -f Makefile_corr
    make -f Makefile_corr

This compiles an executable `corr_test` which can then be run

    ./corr_test

This will output a 2D array `corr.out` which is specified in `evolve_corr.f90` and written out in the subroutine:

    subroutine make_output()
        write(n_file,'(32(ES22.15,2X))') y(SYS_BG_I), sqrt(k2), corr(:,:)
    end subroutine make_output

The first two columns are $\alpha = \ln(a)$ and $H$ in planck units. These are the first two elements of the array y, which has dimension (6).

The important parts are `sqrt(k2)` the seventh column, which as I've said is the wavenumber, and `corr(:,:)`, which is the correlation matrix with the two-point correlation functions for the inflaton $phi$ and the transverse field $chi$ and their canonical momenta. The diagonal is the autocorrelation functions (in k-space) with the third entry being the $\chi$ power spectrum, which is what I'm after.
           _                                                             _
          |  < |φ|^2 >(k)  <π_φ φ^*  >(k)  <χ φ^*  >(k)  <π_χ φ^*  >(k)   |
corr(k) = |  <φ π_φ^*>(k)  < |π_φ|^2 >(k)  <χ π_φ^*>(k)  <π_χ π_φ^*>(k)   |
          |  <φ χ^*  >(k)  <π_φ χ^*  >(k)  < |χ|^2 >(k)  <π_χ χ^*  >(k)   |
          |_ <φ π_χ^*>(k)  <π_φ π_χ^*>(k)  <χ π_χ^*>(k)  < |π_χ|^2 >(k)  _|

These are output for each k as an array of dimension 16:

corr(k) => ( <|φ|^2>(k), <π_φ φ^*>(k), <χ φ^*>(k), <π_χ φ^*>(k), <φ π_φ^*>(k), ... <|π_χ|^2>(k) )

So each column of the array in `corr.out` will be for a new wavenumber `sqrt(k2)`, with the wavenumber specified in the -17th element and the $\chi$ power spectrum specified in the -6th element.

I've written a python script `









So the ansatz once we have <|chi^2|>(k) is to select the spatial $k$ range and then scale it so it's in the same form as the Peak Patch power spectra, and get zeta using a transfer function from Tom's lattice simulations. From there, we can measure $\sigma_{8,exp}$ from it and then taking $P(k) = \sigma_{8,Planck}/\sigma_{8,exp} P_{\zeta\zeta}(k)$.





For my purposes, we're not interested in $\chi$ itself, but its effect on the $\zeta$ field. A transfer function takes $\Delta\phi_f$, the change in the inflaton due to the instability, to $\zeta$, and $\chi$ couples to $\Delta\phi_f$ in a nearly perfect quadratic. The form of this quadratic can be found by running an ensemple of ballistic trajectories through a portion of the potential. This is done numerically in `ballistic_ensemble.f90`. To run this:

    module load intel
    make clean -f Makefile_ens
    make -f Makefile_ens

The output from this is put in a data file `ball_ens.out`. I've written a script `read_ball_ens.py` which fits that curve to a quadratic to the exact form of

$\Delta\phi_f = a \chi^2 + b \chi + c$

Typically, `a` is of order 1, `b` is consistent with 0, and `c` is of order 1. We ignore c though because we only ever work with the mean-zero feilds $\chi-\langle\chi\rangle$ and so forth, so this term will be tossed out anyway, ultimately giving us a form for $\zeta$:

$   \zeta = \zeta_G + T[ a\chi^2 ]   $

