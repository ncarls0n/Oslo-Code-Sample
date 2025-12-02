module bbps_profile

  use cosmology

  implicit none

  integer bbps_model

  real, dimension(3), parameter :: P0_Am     = (/ 18.100,  7.490, 20.7   /)
  real, dimension(3), parameter :: P0_alpham = (/  0.154,  0.226, -0.074 /)
  real, dimension(3), parameter :: P0_alphaz = (/ -0.758, -0.957, -0.743 /)
  
  real, dimension(3), parameter :: xc_Am     = (/  0.497,    0.710,  0.438 /)
  real, dimension(3), parameter :: xc_alpham = (/ -0.00865, -0.0833, 0.011 /)
  real, dimension(3), parameter :: xc_alphaz = (/  0.731,    0.853,  1.01  /)

  real, dimension(3), parameter :: beta_Am     = (/ 4.35,   4.19,  3.82   /)
  real, dimension(3), parameter :: beta_alpham = (/ 0.0393, 0.048, 0.0375 /)
  real, dimension(3), parameter :: beta_alphaz = (/ 0.415,  0.615, 0.535  /)

  real, dimension(3), parameter :: bbps_delta  = (/  200,  500,  200 /)
  real, dimension(3), parameter :: bbps_m2c    = (/ 0.26, 0.15, 0.26 /)

  real, parameter :: alpha = 1
  real, parameter :: gamma = -0.3

  real, dimension(3), parameter :: P0_Am_rho     = (/ 4.0e3,  1.9e4, 1.5e-4 /)
  real, dimension(3), parameter :: P0_alpham_rho = (/  0.29,  0.09, 0.14 /)
  real, dimension(3), parameter :: P0_alphaz_rho = (/ -0.66, -0.95, -1.32 /)

  real, dimension(3), parameter :: alpha_Am_rho     = (/  0.88,    0.70,  0.68 /)
  real, dimension(3), parameter :: alpha_alpham_rho = (/ -0.03, -0.017, -0.02 /)
  real, dimension(3), parameter :: alpha_alphaz_rho = (/  0.19,    0.27, 0.29  /)

  real, dimension(3), parameter :: beta_Am_rho     = (/ 3.83,   4.43,  6.4   /)
  real, dimension(3), parameter :: beta_alpham_rho = (/ 0.04, 0.005, 0.028 /)
  real, dimension(3), parameter :: beta_alphaz_rho = (/ -0.025,  0.037, -0.25  /)

  real, parameter :: xc_rho = 0.5
  real, parameter :: gamma_rho = 0.2

  real, dimension(3), parameter :: P0_Am_dm     = (/ 4.1e3,  2.8e3, 1.5e-4 /)
  real, dimension(3), parameter :: P0_alpham_dm = (/  -0.23,  -0.09, 0.14 /)
  real, dimension(3), parameter :: P0_alphaz_dm = (/ -1.09, -0.81, -1.32 /)

  real, dimension(3), parameter :: alpha_Am_dm     = (/  0.77,    0.9,  0.68 /)
  real, dimension(3), parameter :: alpha_alpham_dm = (/ 0.12, 0.07, -0.02 /)
  real, dimension(3), parameter :: alpha_alphaz_dm = (/  0.37,    0.34, 0.29  /)

  real, dimension(3), parameter :: beta_Am_dm     = (/ 3.75,   3.76,  6.4   /)
  real, dimension(3), parameter :: beta_alpham_dm = (/ 0.007, 0.004, 0.028 /)
  real, dimension(3), parameter :: beta_alphaz_dm = (/ -0.002,  0.024, -0.25  /)

  real, parameter :: xc_dm = 0.5
  real, parameter :: gamma_dm = 1.0


contains

  real function bbps_ptilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! The Compton y-parameter for a given halo is given by
    !   y = y0 * Mhalo * [Omegam(1+z)^3+Omegal] * int Ptilde * dltilde
    ! where
    !   Mhalo  = halo mass in Msun
    !   Ptilde = dimensionless pressure <--> Pressure = Ptilde*P_Delta
    !   ltilde = distance along los in units of virial radius rvir
    !   y0     = 3*sigma_T*(100 km/sec)^2/(8pi*me*c^2) in 1/Msun
    ! see 
    !   Battaglia et al. 2012, ApJ, 758, 75 
    ! section 4.1 for the precise definitions of P_Delta and Ptilde
    ! ----------------------------------------------------------------

    integer i

    real mh,z,x
    real A0,alpham,alphaz
    real P0,xc,beta

    i = bbps_model

    ! P0

    A0     = P0_Am(i)
    alpham = P0_alpham(i)
    alphaz = P0_alphaz(i)

    P0     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! xc

    A0     = xc_Am(i)
    alpham = xc_alpham(i)
    alphaz = xc_alphaz(i)

    xc     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! beta

    A0     = beta_Am(i)
    alpham = beta_alpham(i)
    alphaz = beta_alphaz(i)

    beta   = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! Pressure

    bbps_ptilde  = P0*(x/xc)**gamma*(1+(x/xc)**alpha)**(-1.0*beta)

    return

  end function bbps_ptilde

  real function bbps_rhotilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! The gas profile from BBPS
    ! ---------------------------------------------------------------- 

    integer i

    real mh,z,x
    real A0,alpham,alphaz
    real P0,beta,alpha_rho

    ! debug
    real r200c

    i = bbps_model

    ! P0

    A0     = P0_Am_rho(i)
    alpham = P0_alpham_rho(i)
    alphaz = P0_alphaz_rho(i)

     P0     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! xc

    A0     = alpha_Am_rho(i)
    alpham = alpha_alpham_rho(i)
    alphaz = alpha_alphaz_rho(i)

     alpha_rho     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! beta

    A0     = beta_Am_rho(i)
    alpham = beta_alpham_rho(i)
    alphaz = beta_alphaz_rho(i)

    beta   = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! gas density
    bbps_rhotilde  = P0*(x/xc_rho)**(-1.0*gamma_rho)*(1+(x/xc_rho)**alpha_rho)**(-beta)
    bbps_rhotilde  = bbps_rhotilde * (1+omegal/omegam/(1+z)**3) !now in units of rhobar not rhocrit

    !bbps_rhotilde  = P0*(x/xc_rho)**(-1.0*gamma_rho)*(1+(x/xc_rho)**alpha_rho)**(-1.0*(beta-gamma_rho)/alpha_rho)
    ! SIS instead
!    r200c = (3.*mh/4./3.14159/200./rhocrit(z))**(1./3.) ! comoving Mpc
!    r200c = 1e3 * r200c / (1+z) ! proper kpc
!    bbps_rhotilde = 0.172 * mh / (4.*3.14159/3.*r200c**3) ! msun/kpc^3
!    bbps_rhotilde = bbps_rhotilde / 1e10 /3/x**2 ! 1e10 h^2 msun/kpc^3

   return

  end function bbps_rhotilde

  real function bbps_dmtilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! The DM profile from BBPS
    ! ---------------------------------------------------------------- 

    integer i

    real mh,z,x
    real A0,alpham,alphaz
    real P0,beta,alpha_dm

    i = bbps_model

    ! P0

    A0     = P0_Am_dm(i)
    alpham = P0_alpham_dm(i)
    alphaz = P0_alphaz_dm(i)

    P0     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! xc

    A0     = alpha_Am_dm(i)
    alpham = alpha_alpham_dm(i)
    alphaz = alpha_alphaz_dm(i)

    alpha_dm     = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! beta

    A0     = beta_Am_dm(i)
    alpham = beta_alpham_dm(i)
    alphaz = beta_alphaz_dm(i)

    beta   = A0*(mh/1e14)**alpham*(1+z)**alphaz

    ! Pressure

    bbps_dmtilde  = P0*(x/xc_dm)**(-1.0*gamma_dm)*(1+(x/xc_dm)**alpha_dm)**(-1.0*beta)
    bbps_dmtilde  = bbps_dmtilde * (1+omegal/omegam/(1+z)**3) !now in units of rhobar not rhocrit
!    bbps_dmtilde  = P0*(x/xc_dm)**(-1.0*gamma_dm)*(1+(x/xc_dm)**alpha_dm)**(-1.0*(beta-gamma_dm)/alpha_dm)

    return

  end function bbps_dmtilde

  real function nfw_dmtilde(mh,z,x)

    implicit none

    ! ----------------------------------------------------------------
    ! NFW profile for Lensing. 
    ! Uses concentration parameter from Duffy 2008
    ! units are rho/\bar{rho}
    ! ---------------------------------------------------------------- 


    real mh,z,x
    real C,dcrit

    C = C_duffy(mh,z)
    dcrit = deltacrit(z)

    nfw_dmtilde = dcrit/3 * 1./( log(1+C) - C/(1+C) ) * C**2/(x*(1+C*x)**2)  

    return

  end function nfw_dmtilde

  real function C_duffy(mh,z)

    implicit none

    ! ----------------------------------------------------------------
    ! Concentration parameter from Duffy 2008
    ! ----------------------------------------------------------------

    real mh, z

    C_duffy = 7.85 * (mh / (2e12/h))**(-0.081) * (1+z)**(-0.71)

    return

  end function C_duffy

  real function W_Kappa(z,chi,chist)

    implicit none

    ! ----------------------------------------------------------------                       
    ! Lensing Kernel                                                                         
    ! ----------------------------------------------------------------                       

    real z, chi, chist

    W_Kappa = 3./2 * omegam * (h * 100/ckms)**2 * (1+z) * chi*(1-chi/chist)

    return

  end function W_Kappa

end module bbps_profile
