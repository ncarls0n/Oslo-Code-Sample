module input_parameters

  use cosmoparams ! src/modules/GlobalVariables/cosmoparams.f90

  character(len=512) &
       filein    ,& ! 
       fileout   ,& ! 
       densfilein,& ! name of file for cosmological density field *
       pkfile    ,& ! file (with path) containing tabulated CAMB power spectra used to generate cosmological matter density fields
       Tkfile    ,& ! file (with path) containing PING non-Gaussianity \Delta\phi \to \Delta\zeta transfer function
       filterfile,& ! file (with path) containing filter bank
       fielddir     ! path where cosmological fields are saved for this run once generated

  integer &
       ireadfield      ,& ! boolean switch, 1 to read in initial conditions fields, 0 to generate from a power spectrum *
       ibatch          ,& ! 
       ioutshear       ,& ! boolean switch, 
       num_redshifts   ,& ! boolean switch, Peak Patch has an option of running multiple single-redshift boxes at a range of redshifts, this is used to set how many redshift steps are used (only if ievol is 0). *
       iseedFFT        ,& ! 
       nlx,nly,nlz     ,& ! number of parallel cubes along each axis (currently only nlx=nly=nlz is supported) *
       nbuff           ,& ! Thickness of buffer in voxels *
       next            ,& ! Sidelength of simulation cube in voxels (including buffers) *
       ievol           ,& ! boolean switch, 1 for lightcone runs, 0 for single-redshift runs *
       ihard           ,& ! 
       debug           ,& ! 
       wsmooth         ,& ! 
       ioutfield       ,& ! boolean switch, determines whether to save IC cosmological fields *
       NonGauss        ,& ! selects the early-universe non-Gaussianity model to use *
       ng_seed         ,& ! independent random number seed for non-Gaussian source fields *
       ilpt            ,& ! order of Lagrangian perturbation theory to use (1 for first order, 2 for second order) *
       iwrap           ,& ! 
       iwant_field_part,& ! 
       largerun           ! boolean switch, to split up halo catalogues for very large runs larger runs *

  real &
       global_redshift ,& ! redshift for single redshift run (minimum if num_redshifts>1) *
       maximum_redshift,& ! maximum redshift if num_redshifts>1 *
       fhdmclus        ,& ! 
       dcore_box       ,& ! parallelization cube sidelength (in Mpc, excluding buffers) *
       dL_box          ,& ! simulation sidelength (in Mpc, excluding buffers) *
       cenx,ceny,cenz  ,& ! coordinates of observer (in Mpc) for lightcone runs *
       zinit_fac_ell   ,& ! 
       dFrho           ,& ! 
       rmax2rs         ,& ! 
       fNL             ,& ! non-Gaussianity order for f_NL type early-universe non-Gaussianity *
       A_nG,B_nG,R_nG  ,& ! non-Gaussianity parameters, see python/peak-patch.py for details *
       H_e                ! Hubble parameter at end of instability in PING model non-Gaussianity *

  ! Other parameters

  real    fbuffer_box
  real    biasold,biasnew

  integer nboxes, idocore, idoboxf, ifmt
  integer min_core_box, max_core_box
  logical verbose
  integer*8 ncxm,ncxp,ncym,ncyp,nczm,nczp

end module input_parameters

! *Starred variables are set in the parameter file for a Peak Patch run, see example/param/param.param for details.
