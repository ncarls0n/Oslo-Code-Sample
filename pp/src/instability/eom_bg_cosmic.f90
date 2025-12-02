! eom_bg_cosmic.f90

! Module for the equations of motion of $(\alpha,H,\phi^A,\dot{\phi}^A)$

module eom_bg_cosmic
#include "macros.h"
  use params
  use potential
  
  implicit none

contains

  function deriv_bg_cosm(y, dim, param_r) result(f)
    real(dl), intent(in) :: y(:)
    integer, intent(in) :: dim
    real(dl), optional, intent(in) :: param_r
    real(dl) :: f(dim)

    f(ALPHA_I) = d_alpha_cosm(y)
    f(HUB_I) = d_hub_cosm(y)
    f(PHI_I) = d_f_cosm(y)
    f(DPHI_I) = d_df_cosm(y)
  end function deriv_bg_cosm
  
  function d_alpha_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = y(HUB_I)
  end function d_alpha_cosm

  function d_hub_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f

    f = -0.5_dl*sum(y(DPHI_I)**2)
  end function d_hub_cosm

  function d_f_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f(nfld)

    f = y(DPHI_I)
  end function d_f_cosm

  function d_df_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f(nfld)

    integer i

    do i=1, nfld
       f(i) = -3._dl*y(HUB_I)*y(2+nfld+i) - bg_dv(y(PHI_I),i)
    end do
  end function d_df_cosm

  subroutine set_hub_cosm(y)
    real(dl) :: y(:)

    y(HUB_I) = sqrt((0.5_dl*sum(y(DPHI_I)**2) + bg_v(y(PHI_I)))/3._dl)
  end subroutine set_hub_cosm
  
end module eom_bg_cosmic
