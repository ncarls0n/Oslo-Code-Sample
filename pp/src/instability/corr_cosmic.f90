! corr_cosmic.f90

! Module to compute the 2pt correlator from the eom_pert_cosmic system, and set its initial conditions.

! to do: check that I have things properly in column major order for utility functions
! to do: set ICs with cdxdx -> cdxdx - k^2/a^2. Needs to be tested.
! to do: comment the BLAS routines
! to do: clean up the VARCOND macros for readability, reduce amount of repeated code
! to do: allow for setting \dot{L}_{AB} by analytically differentiating the functional form of L_{AB}
!        and setting cdxdx to give the correct determinant for the correlations of each field/momentum
!        pair.

#define VARCOND 
! When defined use: cdxdx = (\langle\dot{\hat{\chi}}_A,\dot{\hat{\chi}}_B\rangle - k^2/a^2)
! When not defined use: cdxdx = \langle\dot{\hat{\chi}}_A,\dot{\hat{\chi}}_B\rangle

module corr_cosmic
#include "macros.h"
  use eom_bg_cosmic
  use eom_pert_cosmic
  use util_mod
  
  implicit none

contains

  ! Subroutine to set initial conditions for the eom_pert_cosmic system.
  ! to do: check if dsyev stores eigenvectors as rows or columns
  ! to do: check if dpotrf puts zeros in the upper triangle
  ! to do: check that norm2 is correctly multiplied or divided
  subroutine init_pert_cosm(y, k2, norm2)
    real(dl) :: y(:)
    real(dl) :: k2, norm2

    integer :: i
    real(dl), dimension(nfld,nfld) :: l_inv, eig_vec, temp_mat, temp2_mat
    real(dl), dimension(nfld) :: eig_val
    real(dl), dimension(3*nfld-1) :: work

    ! set diagonalized \sqrt{\omega^2_{eff}} = \sqrt{k^2 + a^2D^2 - a^2(2H^2+\dot{H})}
    ! where D^2 is the diagonalized mass matrix
    eig_vec = bg_ddv_mat(y(PHI_I))
    call dsyev('V','L',nfld,eig_vec,nfld,eig_val,work,3*nfld-1,i)  ! diagonalize mass matrix
    do i=1,nfld
       eig_val(i) = sqrt(k2 + exp(2._dl*y(ALPHA_I))*(eig_val(i) - 2._dl*y(HUB_I) - d_hub_cosm(y)))
    end do
    
    ! set L by Choleski decomposition of norm^2U[\sqrt{\omega^2_{eff}}]^{-1/2}U^T/2
    m2_mat = 0._dl
    do i=1,nfld
       m2_mat(i,i) = 0.5_dl/eig_val(i)*norm2
    end do
    call dgemm('N','T',nfld,nfld,nfld,1._dl,m2_mat,nfld,eig_vec,nfld,0._dl,temp_mat,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,eig_vec,nfld,temp_mat,nfld,0._dl,temp2_mat,nfld)
    call dpotrf('L',nfld,temp2_mat,nfld,i)  
    call set_lower_ar(y(L_I), temp2_mat, nfld, 0)  ! set y(L_I)
    
    y(DL_I) = 0._dl  ! set dl = 0
    y(CORR_XDX_I) = 0._dl  ! set cxdx = 0

#ifndef VARCOND
    ! set cdxdx = norm^2/a^{2}L^{-1}U\sqrt{\omega^2_{eff}}U^T{L^{-1}}^T/2
    m2_mat = 0._dl
    do i=1,nfld
       m2_mat(i,i) = 0.5_dl*eig_val(i)*norm2/exp(2._dl*y(ALPHA_I))
    end do
    l_inv = temp2_mat
    call dtrtri('L','N',nfld,l_inv,nfld,i)  ! compute inverse of L
    call dgemm('N','T',nfld,nfld,nfld,1._dl,m2_mat,nfld,eig_vec,nfld,0._dl,temp_mat,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,eig_vec,nfld,temp_mat,nfld,0._dl,temp2_mat,nfld)   
    call dtrmm('L','L','N','N',nfld,nfld,1._dl,l_inv,nfld,temp2_mat,nfld) 
    call dtrmm('R','L','T','N',nfld,nfld,1._dl,l_inv,nfld,temp2_mat,nfld)
    call set_lower_ar(y(CORR_DXDX_I), temp2_mat, nfld, 0)
#endif
#ifdef VARCOND
    ! set cdxdx = norm^2/a^{2}L^{-1}U\sqrt{\omega^2_{eff}}U^T{L^{-1}}^T/2 - k^2/a^2
    ! calc norm^2*L^{-1}U\sqrt{\omega^2_{eff}}U^T{L^{-1}}^T/2
    m2_mat = 0._dl
    do i=1,nfld
       m2_mat(i,i) = 0.5_dl*eig_val(i)*norm2
    end do
    l_inv = temp2_mat
    call dtrtri('L','N',nfld,l_inv,nfld,i)  ! compute inverse of L
    call dgemm('N','T',nfld,nfld,nfld,1._dl,m2_mat,nfld,eig_vec,nfld,0._dl,temp_mat,nfld)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,eig_vec,nfld,temp_mat,nfld,0._dl,temp2_mat,nfld)   
    call dtrmm('L','L','N','N',nfld,nfld,1._dl,l_inv,nfld,temp2_mat,nfld) 
    call dtrmm('R','L','T','N',nfld,nfld,1._dl,l_inv,nfld,temp2_mat,nfld)
    ! calc norm^2*L^{-1}U\sqrt{\omega^2_{eff}}U^T{L^{-1}}^T/2 - k^2
    do i=1,nfld
       temp2_mat(i,i) = temp2_mat(i,i) - k2
    end do
    ! calc (norm^2*L^{-1}U\sqrt{\omega^2_{eff}}U^T{L^{-1}}^T/2 - k^2)/a^2
    temp2_mat = temp2_mat/exp(2._dl*y(ALPHA_I))
    call set_lower_ar(y(CORR_DXDX_I), temp2_mat, nfld, 0)
#endif
  end subroutine init_pert_cosm

  ! Subroutine to get the real part of the 2pt from the eom_pert_cosmic system.
  ! corr is organized with f1 in column 1, df1 in column2, f2 in column 3, ... , like wise with rows
  ! n.b. This sets the correlation for \phi_A without the conformal factor of a, and uses the
  !      conjugate momenta \Pi_A=a^3\dot{\phi_A}.
  subroutine get_2pt_cosm(y, corr, norm2, k2)
    real(dl), intent(in) :: y(:)
    real(dl) :: corr(:,:)
    real(dl), intent(in) :: norm2
    real(dl), intent(in), optional :: k2

    integer :: i,j
    real(dl), dimension(nfld,nfld) :: temp_mat, temp2_mat, temp3_mat
    
    ! set matricies for L, dL, cxdx, cdxdx
    call set_ltri_mat(y(L_I), l_mat, nfld)
    call set_ltri_mat(y(DL_I), dl_mat, nfld)
    call set_asym_mat(y(CORR_XDX_I), cxdx_mat, nfld)
    call set_sym_mat(y(CORR_DXDX_I), cdxdx_mat, nfld)

    ! set \langle\phi^A\phi_B\rangle
    call dgemm('N','T',nfld,nfld,nfld,1._dl,l_mat,nfld,l_mat,nfld,0._dl,temp_mat,nfld)
    corr(CORR_FF_I) = temp_mat/exp(2._dl*y(ALPHA_I))
    
    ! set \langle\phi_A\Pi_B\rangle
    temp_mat = l_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_mat,nfld,cxdx_mat,nfld,-y(HUB_I),temp_mat,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,temp_mat,nfld,l_mat,nfld,0._dl,temp2_mat,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,l_mat,nfld,dl_mat,nfld,1._dl,temp2_mat,nfld)
    corr(CORR_FDF_I) = temp2_mat*exp(y(ALPHA_I))

    ! set \langle\Pi_A\phi_B\rangle
    do j=2,2*nfld,2
       do i=1,2*nfld-1,2
          corr(j,i) = corr(i,j)
       end do
    end do
 
    ! set \langle\Pi^A\Pi^B\rangle
#ifndef VARCOND
    temp_mat = l_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_mat,nfld,cdxdx_mat,nfld,y(HUB_I)**2,temp_mat,nfld)
    temp2_mat = dl_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,dl_mat,nfld,cxdx_mat,nfld,y(HUB_I),temp2_mat,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,temp_mat+temp2_mat,nfld,l_mat,nfld,0._dl,temp3_mat,nfld)
    temp_mat = dl_mat + y(HUB_I)*l_mat
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,l_mat,nfld,cxdx_mat,nfld,1._dl,temp_mat,nfld)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,temp_mat,nfld,dl_mat,nfld,1._dl,temp3_mat,nfld)
    corr(CORR_DFDF_I) = temp3_mat*exp(4._dl*y(ALPHA_I))
#endif
#ifdef VARCOND
    temp_mat = cdxdx_mat
    do i=1,nfld
       temp_mat(i,i) = temp_mat(i,i) + k2*exp(-2._dl*y(ALPHA_I))  ! need to check   
    end do
    ! L*cdxdx + (H^2L - H*dL)
    temp2_mat = y(HUB_I)**2*l_mat - y(HUB_I)*dl_mat
    call dgemm('N','N',nfld,nfld,nfld,1._dl,l_mat,nfld,temp_mat,nfld,1._dl,temp2_mat,nfld)
    ! dL*cxdx + (L*cdxdx + H^2L - H*dL)
    call dgemm('N','N',nfld,nfld,nfld,1._dl,dl_mat,nfld,cxdx_mat,nfld,1._dl,temp2_mat,nfld)
    ! (dL*cxdx + (L*cdxdx + H^2L - H*dL))*L^T
    call dgemm('N','T',nfld,nfld,nfld,1._dl,temp2_mat,nfld,l_mat,nfld,0._dl,temp_mat,nfld)
    ! -L*cxdx + (dL - H*L)
    temp2_mat = dl_mat - y(HUB_I)*l_mat
    call dgemm('N','N',nfld,nfld,nfld,-1._dl,l_mat,nfld,cxdx_mat,nfld,1._dl,temp2_mat,nfld)
    ! (-L*cxdx + (dL - H*L))*dL^T + ((dL*cxdx + (L*cdxdx + H^2L - H*dL))*L^T)
    call dgemm('N','T',nfld,nfld,nfld,1._dl,temp2_mat,nfld,dl_mat,nfld,1._dl,temp_mat,nfld)
    corr(CORR_DFDF_I) = temp_mat*exp(4._dl*y(ALPHA_I))
#endif
    corr(:,:) = corr(:,:)/norm2
  end subroutine get_2pt_cosm
  
end module corr_cosmic
