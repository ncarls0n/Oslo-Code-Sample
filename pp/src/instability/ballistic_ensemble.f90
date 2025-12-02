! ballistic_ensemble.f90

! Program to run an ensemble of ballistic trajectories through a portion of the potential.

program ballistic_ensemble
#include "macros.h"
#include "gl_macros.h"
  use params
  use gl_integrator
  use newton_root
  use eom_bg_cosmic

  implicit none

  integer :: i,j

  integer, parameter :: nstep = 2**6    ! total number of steps
  integer, parameter :: stepout = 2**0  ! number of steps between outputs
  real(dl) :: dt = 1._dl/dble(2**6)     ! time step
  real(dl) :: t                         ! time

  integer, parameter :: nstep_chi = 2**8  ! number of trajectories in ensemble
  real(dl), parameter :: chi_max = 1.e-4  ! max IC excursion of ensemble
  real(dl), parameter :: chi_step = chi_max/dble(nstep_chi)  ! step size of IC in ensemble
  integer :: nstep_ens  ! number of time steps to take during the ensemble
  
  character(len=*), parameter :: f_init = 'ic.out'   ! file name for ICs
  character(len=*), parameter :: f_out = 'ball_ens.out'  ! file name for ouput
  integer, parameter :: n_file = 99

  real(dl), dimension(2,SYS_DIM_BG,nu_gl) :: g = 0._dl ! Gauss-Legendre g vector
  real(dl), dimension(SYS_DIM_BG) :: y                 ! odds are q, evens are dq
  real(dl), dimension(SYS_DIM_BG) :: z                 ! temp state of system

  real(dl), parameter, dimension(SYS_DIM_BG) :: y_0 = (/0._dl,0._dl,10._dl,0.0_dl,-0.815_dl,0._dl/)  ! ICs of y
  real(dl), dimension(SYS_DIM_BG) :: y_attract
  real(dl), dimension(SYS_DIM_BG) :: y_fid
  
  real(dl) :: start_time, end_time

  ! Initialize i/o
  call init_output(n_file)
  
  ! Initialize system to find attractor
  y(:) = y_0(:)
  call set_hub_cosm(y)
  call gl_solve_g(y, deriv_bg_cosm, z, g, dt, SYS_DIM_BG, eps_0, c_max_0)  
  do while (y(3) .gt. phi_p+phi_w)
     call gl_integrate(y, deriv_bg_cosm, z, g, dt, SYS_DIM_BG, eps_0, c_max_0)
  end do
  call gl_newton_root(y, deriv_bg_cosm, get_phi, get_dphi, phi_p+phi_w, z, g, dt, SYS_DIM_BG, eps_0, c_max_0, eps_0, c_max_0)
  y_attract(:) = y(:)
  y_attract(ALPHA_I) = 0._dl
  
  ! Run the case of ICs with $\chi=0$, $\dot{\chi}=0$ to the nearest time step where $\phi < \phi_p - \phi_w$.
  ! Count the steps and use same number of steps for other trajectories, ie output of uniform t.
  nstep_ens = 0
  y(ALPHA_I) = 0._dl
  do while (y(3) .gt. phi_p-phi_w)
     nstep_ens = nstep_ens + 1
     call gl_integrate(y, deriv_bg_cosm, z, g, dt, SYS_DIM_BG, eps_0, c_max_0)
     if (nstep_ens .gt. nstep) exit
  end do
  y_fid(:) = y(:)
  
  ! Start the clock
  call cpu_time(start_time)

  ! Loop setting ICs for $\chi$ and $\dot{\chi}$, and updating $H$ and restting $t$ and  $\alpha$.
  do i=0,nstep_chi
     t = 0
     y(:) = y_attract(:)
     y(2+nfld) = y(2+nfld) + i*chi_step
     call set_hub_cosm(y)
     g = 0._dl
     call gl_solve_g(y, deriv_bg_cosm, z, g, dt, SYS_DIM_BG, eps_0, c_max_0)
     ! Loop evolving system for determined amount of time.
     do j=1,nstep_ens
        call gl_integrate(y, deriv_bg_cosm, z, g, dt, SYS_DIM_BG, eps_0, c_max_0)
        t = t + dt
     end do
     ! Here choose what hypersurface to output on
     call gl_newton_root(y, deriv_bg_cosm, get_hub, d_hub_cosm, y_fid(HUB_I), z, g, dt, SYS_DIM_BG, eps_0, c_max_0, eps_0, c_max_0)
     ! Make output including t, y, and ICs for $\chi$ and $\dot{\chi}$.
     call make_output()
  end do
  
  ! Stop the clock
  call cpu_time(end_time)
  print*, 'elapsed time: ', end_time-start_time
  
contains

  subroutine init_output(n_file)
    integer :: n_file
    open(unit=n_file,file=f_out)
  end subroutine init_output

  subroutine make_output()
    write(n_file,'(30(ES22.15,2X))') t, y(:)
  end subroutine make_output

  function get_hub(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(HUB_I)
  end function get_hub
  
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

  function get_alpha(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(ALPHA_I)
  end function get_alpha
  
end program ballistic_ensemble
