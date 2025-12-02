! eom_pert_cosmic.f90

! Module for equations of motion of L_{A,B}, \dot{L}_{A,B},
! \langle\hat{\chi}_A,\dot{\hat{\chi}}_B\rangle, and \langle\dot{\hat{\chi}}_A,\dot{\hat{\chi}}_B\rangle
! or (\langle\dot{\hat{\chi}}_A,\dot{\hat{\chi}}_B\rangle - k^2/a^2).

! to do: imporove macro readability for variable conditioning

! to do: Write a subtourine to set l_mat, dl_mat, cxdx_mat, cdxdx_mat
! to do: Write a subroutine to set the effective mass matrix

#define VARCOND 
! When defined use: cdxdx = (\langle\dot{\hat{\chi}}_A,\dot{\hat{\chi}}_B\rangle - k^2/a^2)
! When not defined use: cdxdx = \langle\dot{\hat{\chi}}_A,\dot{\hat{\chi}}_B\rangle

module eom_pert_cosmic
#include "macros.h"
  use params
  use potential
  use eom_bg_cosmic
  use util_mod

  implicit none

  real(dl), dimension(nfld,nfld) :: l_mat, dl_mat, cxdx_mat, cdxdx_mat
  real(dl), dimension(nfld,nfld) :: a_mat, b_mat, c_mat, m2_mat
  
contains

  function deriv_pert_cosm(y, dim, param_r) result(f)
    real(dl), intent(in) :: y(:)
    integer, intent(in) :: dim
    real(dl), optional, intent(in) :: param_r
    real(dl) :: f(dim)

    call set_matricies(y,param_r)

    f(SYS_BG_I) = deriv_bg_cosm(y, dim)
    f(L_I) = d_l_cosm(y)
    f(DL_I) = d_dl_cosm(y)
    f(CORR_XDX_I) = d_corr_xdx_cosm(y)
    f(CORR_DXDX_I) = d_corr_dxdx_cosm(y, param_r)
  end function deriv_pert_cosm

  function d_l_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f(nfld*(nfld+1)/2)

    f = y(DL_I)
  end function d_l_cosm

  function d_dl_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f(nfld*(nfld+1)/2)

    real(dl), dimension(nfld,nfld) :: temp_mat

    ! calc L*cdxdx + (C*cxdx + (L*A - B))
    temp_mat = b_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_mat,nfld,a_mat,nfld,-1._dl,temp_mat,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,c_mat,nfld,cxdx_mat,nfld,1._dl,temp_mat,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_mat,nfld,cdxdx_mat,nfld,1._dl,temp_mat,nfld)

    ! set array from matrix
    call set_lower_ar(f,temp_mat,nfld,0)
  end function d_dl_cosm

  function d_corr_xdx_cosm(y) result(f)
    real(dl), intent(in) :: y(:)
    real(dl) :: f(nfld*(nfld-1)/2)

    call set_lower_ar(f,a_mat,nfld,1)
  end function d_corr_xdx_cosm

#ifndef VARCOND
  function d_corr_dxdx_cosm(y, k2) result(f)
    real(dl), intent(in) :: y(:)
    real(dl), optional, intent(in) :: k2
    real(dl) :: f(nfld*(nfld+1)/2)
    
    integer :: i,j
    real(dl), dimension(nfld,nfld) :: temp_mat, temp2_mat, temp3_mat

    ! calc cxdx*cxdx + cdxdx
    temp_mat = cdxdx_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cxdx_mat,nfld,cxdx_mat,nfld,1._dl,temp_mat,nfld)
    ! calc C*(cxdx*cxdx + cdxdx)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,c_mat,nfld,temp_mat,nfld,0._dl,temp2_mat,nfld)
    ! calc  L^{-1}*(C*(cxdx*cxdx + cdxdx))
    temp_mat = l_mat
    call dtrtri('L','N',nfld,temp_mat,nfld,i)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp_mat,nfld,temp2_mat,nfld,0._dl,temp3_mat,nfld)
    ! calc cxdx*cdxdx - (L^{-1}*(C*(cxdx*cxdx + cdxdx)))
    call dgemm('N','N',nfld,nfld,nfld,1._dl,cxdx_mat,nfld,cdxdx_mat,nfld,-1._dl,temp3_mat,nfld)
    ! calc -cxdx*A + (cxdx*cdxdx - (L^{-1}*(C*(cxdx*cxdx + cdxdx))))
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cxdx_mat,nfld,a_mat,nfld,1._dl,temp3_mat,nfld)
    ! symmetrize the matrix
    call symmetrize_mat(temp3_mat, temp_mat, nfld)
    ! set array from matrix
    call set_lower_ar(f,temp_mat,nfld,0)
  end function d_corr_dxdx_cosm
#endif

#ifdef VARCOND
  function d_corr_dxdx_cosm(y, k2) result(f)
    real(dl), intent(in) :: y(:)    
    real(dl), optional, intent(in) :: k2
    real(dl) :: f(nfld*(nfld+1)/2)

    integer :: i,j
    real(dl), dimension(nfld,nfld) :: temp_mat, temp2_mat, temp3_mat

    ! get l_inv
    temp_mat = l_mat
    call dtrtri('L','N',nfld,temp_mat,nfld,i)
    ! zero elements not in lower triangle
    !call fill_utri(temp_mat,nfld,1,0._dl)  ! this is to clear the work from the upper triangle
    ! calc L^{-1}*C
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp_mat,nfld,c_mat,nfld,0._dl,temp2_mat,nfld)
    ! calc L^{-1}*C*cxdx
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp2_mat,nfld,cxdx_mat,nfld,0._dl,temp_mat,nfld)
    ! calc L^{-1}*C*cxdx*cxdx
    call dgemm('N','N',nfld,nfld,nfld,1._dl,temp_mat,nfld,cxdx_mat,nfld,0._dl,temp3_mat,nfld)
    ! calc -cxdx*A - L^{-1}*C*cxdx*cxdx)
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,cxdx_mat,nfld,a_mat,nfld,-1._dl,temp3_mat,nfld)
    ! calc (cxdx - L^{-1})*cdxdx + (-cxdx*A - L^{-1}*C*cxdx*cxdx)
    ! (using the reduced cdxdx -> cdxdx - k^2/a^2)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,(cxdx_mat-temp2_mat),nfld,cdxdx_mat,nfld,1._dl,temp3_mat,nfld)
    ! calc (cxdx - L^{-1}*C + H)*k^2/a^2
    temp_mat = cxdx_mat - temp2_mat
    do i=1,nfld
       temp_mat(i,i) = temp_mat(i,i) + y(HUB_I)
    end do
    temp_mat = temp_mat*k2*exp(-2._dl*y(ALPHA_I))
    ! calc (cxdx - L^{-1}*C + H)*k^2/a^2 + (cxdx - L^{-1})*cdxdx + (-cxdx*A - L^{-1}*C*cxdx*cxdx)
    temp3_mat = temp_mat + temp3_mat
    ! symmetrize the matrix
    call symmetrize_mat(temp3_mat, temp_mat, nfld)
    ! set array from matrix
    call set_lower_ar(f,temp_mat,nfld,0)
  end function d_corr_dxdx_cosm
#endif

#ifndef VARCOND
  ! Make this a subroutine to set a_mat, b_mat, and c_mat
  subroutine set_matricies(y, k2)
    real(dl), intent(in) :: y(:)
    real(dl), intent(in) :: k2

    integer :: i,j
    real(dl), dimension(nfld,nfld) :: l_inv, temp_mat

    call set_ltri_mat(y(L_I), l_mat, nfld)
    call set_ltri_mat(y(DL_I), dl_mat, nfld)
    call set_asym_mat(y(CORR_XDX_I), cxdx_mat, nfld)
    call set_sym_mat(y(CORR_DXDX_I), cdxdx_mat, nfld)

    ! set m2_mat including the k^2/a^2 part
    m2_mat = bg_ddv_mat(y(PHI_I))
    do i=1,nfld
       m2_mat(i,i) = m2_mat(i,i) + k2*exp(-2._dl*y(ALPHA_I))
    end do

    ! calc B
    b_mat(:,:) = -(d_hub_cosm(y) + 2._dl*y(HUB_I)**2)*l_mat + y(HUB_I)*dl_mat
    call dsymm('L','L',nfld,nfld,1._dl,m2_mat,nfld,l_mat,nfld,1._dl,b_mat,nfld)
    ! calc C
    c_mat(:,:) = 2._dl*dl_mat + y(HUB_I)*l_mat

    ! calc A
    ! calc L^{-1}
    l_inv = l_mat
    call dtrtri('L','N',nfld,l_inv,nfld,i)
    ! zero elements not in lower triangle
    !call fill_utri(temp_mat,nfld,1,0._dl)  ! this is to clear the work from the upper triangle
    ! calc C*corr_xdx
    temp_mat = cxdx_mat
    call dtrmm('L','L','N','N',nfld,nfld,1._dl,c_mat,nfld,temp_mat,nfld)
    ! calc M*L - C*corr_xdx
    call dsymm('L','L',nfld,nfld,1._dl,m2_mat,nfld,l_mat,nfld,-1._dl,temp_mat,nfld)
    ! calc L^{-1}*(M*L - C*corr_xdx) - corr_dxdx
    a_mat = cdxdx_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_inv,nfld,temp_mat,nfld,-1._dl,a_mat,nfld)
    ! anti-sym a_mat
    do j=1,nfld
       a_mat(j,j) = 0._dl
       do i=j,nfld
          a_mat(i,j) = -a_mat(j,i)
       end do
    end do
  end subroutine set_matricies
#endif

#ifdef VARCOND
  ! Make this a subroutine to set a_mat, b_mat, and c_mat
  subroutine set_matricies(y, k2)
    real(dl), intent(in) :: y(:)
    real(dl), intent(in) :: k2

    integer :: i,j
    real(dl), dimension(nfld,nfld) :: l_inv, temp_mat

    call set_ltri_mat(y(L_I), l_mat, nfld)
    call set_ltri_mat(y(DL_I), dl_mat, nfld)
    call set_asym_mat(y(CORR_XDX_I), cxdx_mat, nfld)
    call set_sym_mat(y(CORR_DXDX_I), cdxdx_mat, nfld)

    ! calc B
    m2_mat = bg_ddv_mat(y(PHI_I))
    b_mat(:,:) = -(d_hub_cosm(y) + 2._dl*y(HUB_I)**2)*l_mat + y(HUB_I)*dl_mat
    call dsymm('L','L',nfld,nfld,1._dl,m2_mat,nfld,l_mat,nfld,1._dl,b_mat,nfld)

    ! calc C
    c_mat(:,:) = 2._dl*dl_mat + y(HUB_I)*l_mat

    ! calc A
    ! calc L^{-1}
    l_inv = l_mat
    call dtrtri('L','N',nfld,l_inv,nfld,i)
    ! zero elements not in lower triangle
    !call fill_utri(temp_mat,nfld,1,0._dl)  ! this is to clear the work from the upper triangle
    ! calc C*corr_xdx
    temp_mat = cxdx_mat
    call dtrmm('L','L','N','N',nfld,nfld,1._dl,c_mat,nfld,temp_mat,nfld)
    ! calc M*L - C*corr_xdx
    call dsymm('L','L',nfld,nfld,1._dl,m2_mat,nfld,l_mat,nfld,-1._dl,temp_mat,nfld)
    ! calc L^{-1}*(M*L - C*corr_xdx) - corr_dxdx
    a_mat = cdxdx_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_inv,nfld,temp_mat,nfld,-1._dl,a_mat,nfld)
    ! anti-sym a_mat
    do j=1,nfld
       a_mat(j,j) = 0._dl
       do i=j,nfld
          a_mat(i,j) = -a_mat(j,i)
       end do
    end do
  end subroutine set_matricies
#endif
  
end module eom_pert_cosmic
