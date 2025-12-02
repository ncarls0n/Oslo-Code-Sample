! corr_mod.f90

! Module to interface to the linear calculation of correlations.

! to do: make root finding functions private?

module corr_mod
#include "pert_macros.h"
  use params
  use gl_integrator
  use newton_root
  use eom_bg_cosmic
  use eom_pert_cosmic
  use corr_cosmic
  
  implicit none

  integer, parameter :: kos = 2**2            ! oversampling of k modes
  integer, parameter :: nkos = kos*nn         ! number of oversampled modes
  real(dl), parameter :: dkos = dk/dble(kos)  ! spacing of over sampled k modes
  real(dl), parameter :: corr_norm = sqrt(nvol/dx**3/mpl**2) 
  
  integer, parameter :: step_lim = 2**12      ! maximum number of steps to take
  real(dl), private :: dt0 = 1._dl/dble(2**4) ! base time step
  real(dl), private :: dt                     ! adapted time step
  
  real(dl), dimension(SYS_DIM_PERT) :: y_pert
  real(dl), private, dimension(SYS_DIM_PERT) :: z           ! temp state of system
  real(dl), private, dimension(2,SYS_DIM_PERT,nu_gl) :: g   ! Gauss-Legendre g vector
  real(dl), dimension(nkos,2*nfld,2*nfld) :: corr           ! correlation matrix

contains

  ! Subroutine to do the full calculation of setting corr
  ! to do:
  subroutine init_corr(y_bg0, alpha_init, phi_init, alpha_fin, h_fac)
    real(dl), intent(in) :: y_bg0(:)
    real(dl), intent(in) :: alpha_init, phi_init
    real(dl), intent(in) :: alpha_fin
    real(dl), intent(in) :: h_fac        ! number of efolds before horizon crossing to initialize a mode
    
    integer :: i, j
    real(dl) :: k2
    real(dl), dimension(SYS_DIM_BG) :: y_init
    
    ! Set background ICs and evolve to start of lat calculation
    y_pert(:) = 0._dl
    y_pert(SYS_BG_I) = y_bg0(:)
    call set_hub_cosm(y_pert)
    call adapt_dt(y_pert, 0._dl)
    g = 0._dl
    call gl_solve_g(y_pert(SYS_BG_I), deriv_bg_cosm, z(SYS_BG_I), g(:,SYS_BG_I,:), dt, SYS_DIM_BG, eps_0, c_max_0)
    do i=1,step_lim
       call gl_integrate(y_pert(SYS_BG_I), deriv_bg_cosm, z(SYS_BG_I), g(:,SYS_BG_I,:), dt, &
            SYS_DIM_BG, eps_0, c_max_0)
       if (y_pert(3) .lt. phi_init) exit
    end do
    call gl_newton_root(y_pert(SYS_BG_I), deriv_bg_cosm, get_phi, get_dphi, phi_init, z(SYS_BG_I), g, dt, &
         SYS_DIM_BG, eps_0, c_max_0, eps_0, c_max_0/4)
    ! set y_init with \alpha = \alpha_{init} when \phi = \phi_{init}
    y_init = y_bg0
    y_init(ALPHA_I) = y_init(ALPHA_I) - y_pert(ALPHA_I) + alpha_init
    call set_hub_cosm(y_init)

    ! Loop over modes to set corr
    corr = 0._dl
    do i=1, nkos ! put in max and min (can start at kos if below is set to zero)
       k2 = dble(i**2)*dkos**2
       ! evolve background to fixed k/(aH)
       y_pert(SYS_BG_I) = y_init(:)
       g = 0._dl
       call adapt_dt(y_pert,0._dl)
       do j=1,step_lim
          call gl_integrate(y_pert(SYS_BG_I), deriv_bg_cosm, z(SYS_BG_I), g(:,SYS_BG_I,:), dt, &
               SYS_DIM_BG, eps_0, c_max_0)
          if (get_hor(y_pert) .gt. 0.5_dl*log(k2)-h_fac) exit
       end do
       call gl_newton_root(y_pert(SYS_BG_I), deriv_bg_cosm, get_hor, get_dhor, 0.5_dl*log(k2)-h_fac, z(SYS_BG_I), g, dt, &
            SYS_DIM_BG, eps_0, c_max_0, eps_0, c_max_0/2)
       y_init(:) = y_pert(SYS_BG_I)
     
       ! evolve pert to initialization point (phi_init or alpha = 0)
       call init_pert_cosm(y_pert, k2, sqrt(k2))
       call adapt_dt(y_pert, sqrt(k2))
       g = 0._dl
       call gl_solve_g(y_pert, deriv_pert_cosm, z, g, dt, SYS_DIM_PERT, eps_0, 2*c_max_0, k2)
       do j=1,step_lim
          call gl_integrate(y_pert, deriv_pert_cosm, z, g, dt, SYS_DIM_PERT, eps_0, c_max_0+4, k2)
          call adapt_dt(y_pert, sqrt(k2))
          if (y_pert(ALPHA_I) .gt. alpha_fin) exit
          !if (y_pert(3) .gt. phi_fin) exit
       end do
       call gl_newton_root(y_pert, deriv_pert_cosm, get_alpha, d_alpha_cosm, alpha_fin, z, g, dt, &
            SYS_DIM_PERT, eps_0, c_max_0, eps_0, c_max_0/2, k2)
       !call gl_newton_root(y_pert, deriv_pert_cosm, get_phi, get_dphi, phi_fin, z, g, dt, &
       !     SYS_DIM_PERT, eps_0, c_max_0, eps_0, c_max_0/2, k2)
       call get_2pt_cosm(y_pert, corr(i,:,:), sqrt(k2), k2)
    end do
  end subroutine init_corr

  ! to do: make this private
  subroutine adapt_dt(y,k)
    real(dl), intent(in) :: y(:)
    real(dl), intent(in) :: k

    integer :: i
    real(dl) :: time_scale
    real(dl), dimension(nfld,nfld) :: m2
    real(dl), dimension(nfld) :: eig
    real(dl), dimension((nfld+2)*nfld) :: work

    m2 = bg_ddv_mat(y(PHI_I))
    call dsyev('N','L',nfld,m2,nfld,eig,work,(nfld+2)*nfld,i)

    ! eigenvalues are sorted, so only need to check max and min
    time_scale = max(y(HUB_I), sqrt(abs(k**2/exp(2._dl*y(ALPHA_I)) + eig(1))), &
         sqrt(abs(k**2/exp(2._dl*y(ALPHA_I)) + eig(nfld))))    
    dt = dt0/time_scale

    ! check future timescale
    m2 = bg_ddv_mat(y(PHI_I) + dt*d_f_cosm(y))
    call dsyev('N','L',nfld,m2,nfld,eig,work,(nfld+2)*nfld,i)
    time_scale = max(time_scale, sqrt(abs(k**2/exp(2._dl*(y(ALPHA_I)+dt*d_alpha_cosm(y))) + eig(1))), &
         sqrt(abs(k**2/exp(2._dl*(y(ALPHA_I)+dt*d_alpha_cosm(y))) + eig(nfld))))
    dt = dt0/time_scale
  end subroutine adapt_dt

  ! Functions for root finding
  function get_alpha(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(ALPHA_I)
  end function get_alpha

  function get_phi(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(3)
  end function get_phi

  function get_dphi(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(3+nfld)
  end function get_dphi

  ! \ln(aH)
  function get_hor(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(ALPHA_I) + log(y(HUB_I))
  end function get_hor

  ! \dot{\ln(aH)}
  function get_dhor(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(HUB_I) + d_hub_cosm(y)/y(HUB_I)
  end function get_dhor

  ! Output the correlation matrix
  subroutine write_corr(n_file)
    integer, intent(in) :: n_file
    integer :: i
    
    open(unit=n_file, file='corr.out')
    do i=1,nkos
       write(n_file, '(30(ES22.15, 2x))') i*dkos, corr(i,:,:)*corr_norm**2
    end do
    close(n_file)
  end subroutine write_corr
  
end module corr_mod
