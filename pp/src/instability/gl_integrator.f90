! gl_integrator.f90

module gl_integrator
#include "gl_macros.h"
  use params
  !use eom
  use butcher_table

  implicit none

  interface
     function d_template(y,dim,param_r) result(f)
       import
       real(dl), intent(in) :: y(:)
       integer, intent(in) :: dim
       real(dl), optional, intent(in) :: param_r
       real(dl) :: f(dim)
     end function d_template
  end interface
  
  integer, parameter :: nu_gl = GL_ORDER/2
  real(dl), parameter, dimension(nu_gl,nu_gl) :: a_gl = GL_SET_A
  real(dl), parameter, dimension(nu_gl) :: b_gl = GL_SET_B

  real(dl) :: eps_0 = 1._dl/dble(2**32)
  integer :: c_max_0 = 12
  
contains

  ! Subroutine to perfom one integration step
!  subroutine gl_integrate(y, f, z, g, h, dim, eps, c_max)
  subroutine gl_integrate(y, f, z, g, h, dim, eps, c_max, param_r)
    procedure(d_template) :: f
    real(dl) :: y(:), z(:)
    real(dl) :: g(:,:,:)
    real(dl) :: h
    integer :: dim
    real(dl), intent(in) :: eps
    integer, intent(in) :: c_max
    real(dl), optional, intent(in) :: param_r

    integer :: i

    call gl_solve_g(y,f,z,g,h,dim,eps,c_max,param_r)
    z(:) = 0._dl
    do i=1, nu_gl
       z(:) = z(:) + b_gl(i)*g(2,:,i)
    end do
    y = y + h*z
  end subroutine gl_integrate

  ! Subroutine to solve for the g vector
  ! Two ways to do this:
  ! i) recurse and compare to an epsilon
  ! ii) recurse a fixed number of times
  !  subroutine gl_solve_g(y, f, z, g, h, dim, eps, c_max)
  subroutine gl_solve_g(y, f, z, g, h, dim, eps, c_max, param_r)
    procedure(d_template) :: f
    real(dl), intent(in) :: y(:)
    real(dl) :: z(:)
    real(dl) :: g(:,:,:)
    real(dl), intent(in) :: h
    integer, intent(in) :: dim
    real(dl), intent(in) :: eps
    integer, intent(in) :: c_max
    real(dl), optional, intent(in) :: param_r
    
    integer :: i

    ! g(:,:,:) = 0._dl  ! removed reset of g to provide better initial guess
#if RECURSEFIX
    ! perform a fixed number of recursions
    do i=1,c_max
       call gl_recurse_g(y, f, z, g, h, dim, param_r)
#if CONVTEST
       print*, 'CONVERGENCE: ', maxval(abs(g(2,:,:) - g(1,:,:)))
#endif
    end do
#else
    ! recurse until an error threashold is met, or maximum recursions exceeded
#endif
#if CONVTEST
    if (maxval(abs(g(2,:,:) - g(1,:,:))) .GT. eps) then
       print*, 'CONVERGENCE FAILURE', maxval(abs(g(2,:,:) - g(1,:,:)))
    end if
#endif
  end subroutine gl_solve_g

  ! Subroutine to perform one recursion on the g vector
  subroutine gl_recurse_g(y, f, z, g, h, dim, param_r)
    procedure(d_template) :: f
    real(dl), intent(in) :: y(:)
    real(dl) :: z(:)
    real(dl) :: g(:,:,:)
    real(dl), intent(in) :: h
    integer, intent(in) :: dim
    real(dl), optional, intent(in) :: param_r
    
    integer :: i,j,k

    g(1,:,:) = g(2,:,:)
    do i=1,nu_gl
       z(:) = y(:)
       do j=1,nu_gl
          z(:) = z(:) + h*a_gl(i,j)*g(1,:,j)
       end do
       g(2,:,i) = f(z(:),dim,param_r)
    end do
  end subroutine gl_recurse_g
  
end module gl_integrator
