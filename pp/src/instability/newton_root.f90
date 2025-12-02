! newton_root.f90

module newton_root
#include "gl_macros.h"
  use params
  use gl_integrator

  implicit none

  interface
     function root_template(y) result(f)
       import
       real(dl), intent(in) :: y(:)
       real(dl) :: f
     end function root_template
  end interface

contains

  ! Subroutine to evolve the system y until intercept==root(y).
  subroutine gl_newton_root(y, f, root, droot, intercept, z, g, h, dim, eps, c_max, eps_newt, c_newt, param_r)
    procedure(d_template) :: f
    procedure(root_template) :: root, droot
    real(dl) :: intercept
    real(dl) :: y(:), z(:)
    real(dl) :: g(:,:,:)
    real(dl), intent(in) :: h
    integer, intent(in) :: dim
    real(dl), intent(in) :: eps, eps_newt
    integer, intent(in) :: c_max, c_newt
    real(dl), optional, intent(in) :: param_r

    real(dl) :: h_temp
    integer :: i
    ! iterate Newton's method c_newt times
    do i=1, c_newt
       h_temp = get_newton_h(y, root, droot, intercept)
       if (abs(h) .lt. abs(h_temp)) then
          h_temp = sign(h, h_temp)
       end if
       call gl_integrate(y, f, z, g, h_temp, dim, eps, c_max, param_r)
    end do
#if CONVTEST_NEWT
    if (abs(intercept-root(y)) .gt. eps_newt) then
       print*, 'NEWTON ROOT CONVERENCE', (root(y)-intercept)
    end if
#endif
  end subroutine gl_newton_root

  ! Subroutine to update the time step used in Newton's method.
  function get_newton_h(y, root, droot, intercept) result(h)
    procedure(root_template) :: root, droot
    real(dl), intent(in) :: intercept
    real(dl), intent(in) :: y(:)
    real :: h

    h = (intercept-root(y))/droot(y)
  end function get_newton_h
  
end module newton_root
