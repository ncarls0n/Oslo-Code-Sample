module integrate_profiles

  use cosmology
  use bbps_profile

  implicit none

  abstract interface
     function func_ptr(aa,bb,cc)
       real :: func_ptr
       real, intent (in) :: aa, bb, cc
     end function func_ptr
  end interface

  procedure (func_ptr), pointer :: f_ptr => null ()

  real, parameter :: ds = 1e-2, dx = 1e-2, xmax = 4

contains

  real function integrate_profile(theta,z,mh,halo_profile)

    implicit none

    integer i,n,halo_profile
    real theta,z,mh
    real Delta,chi,rvir,b,s,r,x,s0
    real y,rmax

    if(halo_profile==1) f_ptr => bbps_ptilde
    if(halo_profile==2) f_ptr => bbps_rhotilde
    if(halo_profile==3) f_ptr => bbps_dmtilde

    rmax = 4

    integrate_profile = 0

    Delta = bbps_delta(bbps_model)

    chi = rofz(z)
    rvir = (3*mh/4/pi/Delta/rhocrit(z))**(1./3.)     

    b = chi * sin(theta) / rvir 

    if(b>rmax) return

    s0 = sqrt(rmax**2-b**2) ! this is for a spherical cut at rmax
!    s0 = rmax              ! this is for a cylindrical cut at rmax

    n = int(s0/ds) + 1

    y = 0
    do i=1,n
       s = (i-0.5) * ds
       x = sqrt(s**2+b**2)
       y = y + f_ptr(mh,z,x) !integrate over profile specified by pointer
    enddo

    y = y * ds

    ! factor of two is because we only integrate half of it (int_0^rmax)
    integrate_profile = 2*y 

    return 

  end function integrate_profile

  !OLD profile, kept here for now
  real function integrate_ptilde_ell(ell,z,mh)

    implicit none

    integer i,n
    real ell,z,mh
    real Delta,chi,rvir,x
    real l2yl2

    Delta = bbps_delta(bbps_model)
    chi = rofz(z)
    rvir = (3*mh/4/pi/Delta/rhocrit(z))**(1./3.)     

    l2yl2 = 0
    x = 0
    do while(x<xmax)

       x     = x + dx/2
       l2yl2 = l2yl2 + x * bbps_ptilde(mh,z,x) * sin(ell*x*rvir/chi)
       x     = x + dx/2

    enddo
    l2yl2 = l2yl2**2

    integrate_ptilde_ell = l2yl2

    return 

  end function integrate_ptilde_ell

end module integrate_profiles
