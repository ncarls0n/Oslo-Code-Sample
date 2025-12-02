! lat_init_mod.f90

! Module for subroutines relating to initializing the lattice.

! to do: write subroutine that initializes fields from a supplied correlation
!        matrix and self consitantly initializes metric perturbations from
!        the same random variables
! to do: write a subroutine that will initialize zeta from the metric
!        perturbations
! to do: write a subroutine to call the evolution of the background, field
!        perturbations, and metric pertubations
! to do: probably best to move the zeta initialization stuff elsewhere
! to do: write a subroutine to initialize from a checkpoint

!#define LATINITOPT 1
! LATINITOPT 0 for homogeneous initialization
! LATINITOPT 1 for spectral initialization
! LATINITOPT 2 for convolution initialization
! LATINITOPT 3 for file initialization
! LATINITOPT 4 for checkpoint initialization
!#define ZETAINIT 0
! ZETAINIT 0 to initialize zeta to 0
! ZETAINIT 1 to initialize zeta from constraint equation
!          assuming longitudinal gauge

module lat_init_mod
#include "macros.h"
  use grv_mod
  use fftw3
  use omp_lib
  use params
  use potential_mod
  use analysis

  implicit none

  interface
     function filt_template(k_in, param_r) result(f)
       import
       real(dl) :: k_in
       real(dl) :: param_r(:)
       real(dl) :: f
     end function filt_template
  end interface
     
contains

#if 0
  ! Subroutine that wraps various options for lattice initialization
  ! to do: update inputs
  !subroutine lat_init(corr_in, k_cut, seed_in)
  subroutine lat_init(corr_in, kos, k_cut, seed_in, f0_in, df0_in, H0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
    !real(dl) :: corr_in(:,:,:)!corr_in(2*nfld,2*nfld,fp_nos)
    !real(dl), intent(in) :: k_cut
    !integer, intent(in) :: seed_in
    real(dl) :: corr_in(:,:,:)      ! input correlation !corr_in(2*nfld,2*nfld,cp_nos)
    integer, intent(in) :: kos     ! over sampling of k in corr_in
    real(dl), intent(in) :: k_cut   ! value of k cutoff
    integer, intent(in) :: seed_in  ! seed for rng
    real(dl), intent(in) :: f0_in(:)  ! mean fields inputs
    real(dl), intent(in) :: df0_in(:) ! mean momenta inputs
    real(dl), intent(in) :: H0_in   ! Hubble input
    real(dl) :: f_init(:,:,:,:)     ! array of fields to be initialized
    real(dl) :: df_init(:,:,:,:)    ! array of momenta to be initialized
    real(dl) :: zeta_init(:,:,:)    ! zeta to be initialized
    real(C_DOUBLE), pointer :: f(:,:,:) ! pointer field for lattice
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:) ! pointer for fft
    type(C_PTR) :: planb
    real(dl), intent(in) :: norm_in

#if (LATINITOPT == 0)
    call lat_init_homo()
#elif (LATINITOPT == 1)
    call lat_init_spec(corr_in, kos, k_cut, seed_in, f0_in, df0_in, H0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
#elif (LATINITOPT == 2)
    call lat_init_convolution(corr_in, k_cut, seed_in)
#endif
  end subroutine lat_init
#endif

  
  ! Subroutine to initialize homogeneous fields on the lattice from input. 
  subroutine lat_init_homo(f0_in, df0_in,f_init, df_init)
    real(dl), intent(in) :: f0_in(:), df0_in(:)
    real(dl) :: f_init(:,:,:,:), df_init(:,:,:,:)

    integer :: n

    do n=1,nfld
      f_init(FLDIND(n)) = f0_in(n)
      df_init(FLDIND(n)) = df0_in(n)
    enddo
  end subroutine lat_init_homo

  ! to do: should probably delete corr_Fk after calling
  ! to do: make sure my index ordering is consistant here and in pert_cosmo_test.f90
  ! n.b. zeta mode realization assumes yscl=1
  ! n.b. wrap fields after initializing if using a discrete stencil
  ! n.b. This should corr_in[0:nfld,0:nfld] = <phi_A,phi_B>, 
  !      corr_in[nfld+1:2*nfld, nfld+1:2*nfld] = <Pi_A,Pi_B>,
  !      corr_in[0:nfld, nfld+1:2*nfld] = <phi_A,Pi_B>,
  ! to do: Figure out what is happening with passing arrays of bounds starting at 0 into this
  !        subroutine. may be indexing from 1 internally.
  ! to do: off by a factor of sqrt(2)
  ! to do: check storing in fftw
#if 0
  subroutine lat_init_spec(corr_in, kos, k_cut, seed_in, f0_in, df0_in, H0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
    real(dl) :: corr_in(:,:,:)      ! input correlation !corr_in(2*nfld,2*nfld,cp_nos)
    integer, intent(in) :: kos     ! over sampling of k in corr_in
    real(dl), intent(in) :: k_cut   ! value of k cutoff
    integer, intent(in) :: seed_in  ! seed for rng
    real(dl), intent(in) :: f0_in(:)  ! mean fields inputs
    real(dl), intent(in) :: df0_in(:) ! mean momenta inputs
    real(dl), intent(in) :: H0_in   ! Hubble input
    real(dl) :: f_init(:,:,:,:)     ! array of fields to be initialized
    real(dl) :: df_init(:,:,:,:)    ! array of momenta to be initialized
    real(dl) :: zeta_init(:,:,:)    ! zeta to be initialized
    real(C_DOUBLE), pointer :: f(:,:,:)  ! pointer field for lattice
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)  ! pointer for fft
    type(C_PTR) :: planb  ! fft c2r plan
    real(dl), intent(in) :: norm_in

    integer :: i, j, k
    integer :: ii, jj, kk, l
    real(dl) :: rad
    integer :: n, m

    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: corr_fac
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: corr_Fk

    complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv
   
    ! initialize rng
    call init_rng(seed_in)    

    ! loop over wavenumbers
    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
      do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif 
        do i=1,nnx; ii=i-1
          rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))
          if (rad*dk .GE. k_cut) then
            corr_Fk(:,LATIND) = 0._dl*grv(:) ! cut field fluctuations
          else
            ! initialize corr_fac
            l = floor(rad*kos)
            if (l>0) then
              corr_fac(:,:) = (1._dl,0._dl)*((1._dl + l - rad*kos)*corr_in(:,:,l) & !((1._dl-l+rad*kos)*corr_in(:,:,l) &
                                             + (rad*kos - l)*corr_in(:,:,l+1))
            else
              corr_fac(:,:) = (1._dl,0._dl)*((rad-1._dl)*corr_in(:,:,2) &
                                             + (2._dl-rad)*corr_in(:,:,1))
            endif
            ! Choleski factorization of corr_fac
            call zpotrf('L',2*nfld,corr_fac,2*nfld,l)
            ! print warning if factorization failed
            if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1)) then
              print*, "Factorization Warning: l = ", l
              print*, "ii, jj, kk = ", ii, jj, kk
              print*, "corr_fac = ", corr_fac
            endif
            ! call Gaussian random variables
            if (i==1 .and. j==1 .and. k==1) then
              grv(:) = (0._dl, 0._dl)
            !elseif (i==1 .or. i==nx) then
            !  grv(:) = (1._dl, 0._dl)*get_grv_real(2*nfld)
            else
              grv(:) = get_grv_complex(2*nfld)
            endif
            ! enforce triangular form on correlation matrix
            do m=1,2*nfld; do n=1,2*nfld
              if (n<m) then
                corr_fac(n,m) = (0._dl, 0._dl)
              endif
            enddo; enddo
            ! realize field modes
            !print*, 'corr_fac(2,2): ', corr_fac(2,2)
            !print*, 'grv(:): ', grv(:)
            call ztrmv('L','N','N', 2*nfld, corr_fac, 2*nfld, grv, 1)
            corr_Fk(:,LATIND) = grv(:)/nvol*norm_in
          endif
          ! cut modes above kcut
          !if (rad .GE. k_cut) then
          !  corr_Fk(:,LATIND) = 0._dl*grv(:) ! cut field fluctuations
          !endif
#if (ZETAINIT == 1)
          ! Realize modes for zeta
          Fk(LATIND) = (0._dl,0._dl)
          do n=1,nfld
            Fk(LATIND) = Fk(LATIND) + ((bg_dv(f0_in,n) + 3._dl*H0_in*df0_in(n))/(sum(df0_in**2) - 2._dl*rad) &
                                       + (bg_dv(f0_in,n))/(3._dl*sum(df0_in**2)))*corr_Fk(2*n-1,LATIND)
            Fk(LATIND) = Fk(LATIND) + ((df0_in(n))/(sum(df0_in**2) - 2._dl*rad) &
                                       + (df0_in(n))/(3._dl*sum(df0_in**2)))*corr_Fk(2*n,LATIND)
          enddo
#endif
        enddo
      enddo
    enddo
#if (ZETAINIT == 0)
    zeta_init(:,:,:) = 0._dl
#elif (ZETAINIT == 1)
    ! Invert FFT to initialize zeta
    call fftw_execute_dft_c2r(planb, Fk, f)
    zeta_init(:,:,:) = f(:,:,:)
#endif

    ! initialize homogeneous fields
    call lat_init_homo(f0_in,df0_in,f_init, df_init)
    ! invert FFT, add fluctuations to homogeneous fields
    do m=1,nfld
      Fk(:,:,:) = corr_Fk(m,:,:,:)!Fk(:,:,:) = corr_Fk(2*m-1,:,:,:)
      call fftw_execute_dft_c2r(planb, Fk, f)
      f_init(m,IRANGE) = f_init(m,IRANGE) + f(IRANGE)
      Fk(:,:,:) = corr_Fk(nfld+m,:,:,:)!Fk(:,:,:) = corr_Fk(2*m,:,:,:)
      call fftw_execute_dft_c2r(planb, Fk, f)
      df_init(m,IRANGE) = df_init(m,IRANGE) + f(IRANGE)
    enddo
    ! set scale factor
    !yscl = 1._dl
  end subroutine lat_init_spec
#endif
  
  ! Subroutine to initialize lattice from a supplied Fourier space correlations. Initialization is done to match
  ! Fourier space correlations, not real space correlations.
  ! n.b. the homogeneous values for the fields and momenta must be set separately
  ! n.b. yscl and ysclp must be set separately
  ! n.b. corr should be supplied as the contiuum correlation matrix
  !      (modulo the factor of (2\pi)^3\delta^{(3)}(k-k')).
  ! n.b. I am initializing with corr having the ff corrlations in the upper left
  subroutine lat_init_fluc_spec(corr, filt, kos, kcut, kfilt, seed_in, f_init, df_init, f, FK, planb, norm)
    real(dl), intent(in) :: corr(:,:,:)  ! correlation matrix
    procedure(filt_template) :: filt     ! filter function
    integer, intent(in) :: kos           ! oversamping ratio of k in corr relative to dk on the lattice
    integer, intent(in) :: kcut          ! maximum k index of corr, zero all above
    real(dl), intent(in) :: kfilt(:)     ! filter parameters
    integer, intent(in) :: seed_in       ! seed for rng
    real(dl) :: f_init(:,:,:,:), df_init(:,:,:,:)
    real(C_DOUBLE), pointer :: f(:,:,:)              ! pointer field
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)  ! pointer for fft
    type(C_PTR) :: planb                             ! fft c2r plan
    real(dl), intent(in) :: norm
    
    integer :: i, j, k, n
    integer :: ii, jj, kk, l
    real(dl) :: rad
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: corr_interp
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: corr_Fk
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv

    ! initialize rng
    call init_rng(seed_in)    

    ! loop over wave numbers
    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else; kk=k-1; endif
       do j=1,ny; if (j>nny) then; jj = ny+1-j; else; jj=j-1; endif 
          do i=1,nnx; ii=i-1
             rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))  ! find radius
             l = floor(rad*kos)  ! get radius index
             if (l == 0 .or. l .gt. kcut) then  ! zero the zero/cut mode and cycle the loop
                corr_Fk(:,LATIND) = (0._dl, 0._dl)  
                cycle
             else  ! interpolate corr
                corr_interp(:,:) = (1._dl,0._dl)*((1._dl + l - rad*kos)*corr(l,:,:) + (rad*kos - l)*corr(l+1,:,:))
             endif
             call zpotrf('L',2*nfld,corr_interp,2*nfld,l)  ! Choleski decomposition
             grv(:) = get_grv_complex(2*nfld)
             call ztrmv('L','N','N', 2*nfld, corr_interp, 2*nfld, grv, 1)  ! multiply grv and Choleski
             corr_Fk(:,LATIND) = grv(:)*filt(dk*rad, kfilt)  ! apply filter
          end do
       end do
    end do

    ! Loops to set the ii=0 components to be complex conjugate
    do k=1,nz; if (k.ne.1) then; kk = nz+2-k; else; kk=k; endif
       do j=2,nny; jj = ny+2-j
          corr_FK(:,1,jj,kk) = conjg(corr_Fk(:,1,j,k))
       end do
    end do

    do k=2,nnz; kk = nz+2-k
       corr_FK(:,1,1,kk) = conjg(corr_Fk(:,1,1,k))
    end do

    corr_Fk = corr_Fk*norm/nvol  ! apply normalization
    
    ! invert FFT, add fluctuations to homogeneous fields
    do n=1,nfld
       Fk(:,:,:) = corr_Fk(n,:,:,:)
       call fftw_execute_dft_c2r(planb, Fk, f)
       f_init(FLDIND(n)) = f_init(FLDIND(n)) + f(IRANGE)
       Fk(:,:,:) = corr_Fk(nfld+n,:,:,:)
       call fftw_execute_dft_c2r(planb, Fk, f)
       df_init(FLDIND(n)) = df_init(FLDIND(n)) + f(IRANGE)
    end do
  end subroutine lat_init_fluc_spec

  ! Function for a top_hat filter
  function top_hat(k_in, param_r) result(f)
    real(dl) :: k_in
    real(dl) :: param_r(:)
    real(dl) :: f

    f = 0.5_dl + sign(0.5_dl, param_r(1)-k_in)
  end function top_hat

  ! Function for a bowler hat filter with a cubic tail
  ! n.b. must have param_r(1) .lt. param_r(2)
  function bowler_hat_cubic(k_in, param_r) result(f)
    real(dl) :: k_in
    real(dl) :: param_r(:)
    real(dl) :: f

    real(dl) :: temp
    
    temp = (2._dl*k_in - (param_r(1) + param_r(2))) / ((param_r(2) - param_r(1)))
    f = 0.5_dl + sign(0.5_dl, param_r(1)-k_in) &
         + (sign(0.5_dl, k_in-param_r(1)) + sign(0.5_dl, param_r(2)-k_in)) &
         * (0.5_dl + 0.25_dl*(temp**3 - 3._dl*temp))
  end function bowler_hat_cubic
    
#ifdef CONV_TEST
  ! Subroutine for initializing a lattice for convergence tests.
  ! This subroutine takes the lattice parameters for a smaller run and generates
  ! initial conditions with an identical set of modes.
  ! to do: test this subroutine
  ! to do: Take as input the number of lattice sites on the smaller lattice
  ! to do: When looping over modes to initialize choose only modes which are
  !        initialized on the smaller lattice.
  ! to do: Loop over modes in an identical order to the order of initializeation
  !        on the smaller lattice.
  ! to do: make sure the correct normalization is used give then corr_in calculation
  !        and lattice size
  subroutine lat_init_spec_convergence_test(corr_in, kos, k_cut, nsize_small, seed_in, f0_in, df0_in, H0_in, f_init, df_init, zeta_init, f, Fk, planb, norm_in)
    real(dl) :: corr_in(:,:,:)      ! input correlation !corr_in(2*nfld,2*nfld,cp_nos)
    integer, intent(in) :: kos     ! over sampling of k in corr_in
    real(dl), intent(in) :: k_cut   ! value of k cutoff
    integer, dimension(1:3), intent(in) :: nsize_small  ! size of the lattice which is being convergence tested
    integer, intent(in) :: seed_in  ! seed for rng
    real(dl), intent(in) :: f0_in(:)  ! mean fields inputs
    real(dl), intent(in) :: df0_in(:) ! mean momenta inputs
    real(dl), intent(in) :: H0_in   ! Hubble input
    real(dl) :: f_init(:,:,:,:)     ! array of fields to be initialized
    real(dl) :: df_init(:,:,:,:)    ! array of momenta to be initialized
    real(dl) :: zeta_init(:,:,:)    ! zeta to be initialized
    real(C_DOUBLE), pointer :: f(:,:,:) ! pointer field for lattice
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:,:,:)  ! pointer for fft
    type(C_PTR) :: planb ! fft c2r plan
    real(dl), intent(in) :: norm_in

    real(dl) :: rad           ! radius
    integer :: l              ! radial floor index
    integer :: i, j, k        ! loop indicies
    integer :: ii, jj, kk     ! wave numbers
    integer :: iii, jjj, kkk  ! FFT indicies
    integer :: n1, n2, n3     ! small lattice size
    integer :: nn1, nn2, nn3  ! small lattice Nyquist
    integer :: n, m           ! indexing integers

    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: corr_fac
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: corr_Fk
    complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv

    ! loop over wavenumbers
    n1 = nsize_small(1); n2 = nsize_small(2); n3 = nsize_small(3)
    nn1 = n1/2+1; nn2=n2/2+1; nn3=n3/2+1
    do k=1,n3; if (k>nn3) then; kk = n3+1-k; kkk = nz+1-kk; else; kk = k-1; kkk = k; endif
      do j=1,n2; if (j>nn2) then; jj = n2+1-j; jjj = ny+1-jj; else; jj=j-1; jjj = j; endif 
        do i=1,nn1; ii=i-1; iii = i
           ! get wavenumber
           rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))
           ! initialize corr_fac
           if (rad*dk .GE. k_cut) then
              corr_Fk(:,iii,jjj,kkk) = 0._dl*grv(:) ! cut field fluctuations
           else
              l = floor(rad*kos)
              if (l>0) then
                 corr_fac(:,:) = (1._dl,0._dl)*((1._dl + l - rad*kos)*corr_in(:,:,l) & !((1._dl-l+rad*kos)*corr_in(:,:,l) &
                      + (rad*kos - l)*corr_in(:,:,l+1))
              else
                 corr_fac(:,:) = (1._dl,0._dl)*((rad-1._dl)*corr_in(:,:,2) &
                      + (2._dl-rad)*corr_in(:,:,1))
              endif
              ! Choleski factorization of corr_fac
              call zpotrf('L',2*nfld,corr_fac,2*nfld,l)
              ! print warning if factorization failed
              if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1)) then
                 print*, "Factorization Warning: l = ", l
                 print*, "ii, jj, kk = ", ii, jj, kk
                 print*, "corr_fac = ", corr_fac
              endif
              ! call Gaussian random variables
              if (i==1 .and. j==1 .and. k==1) then
                 grv(:) = (0._dl, 0._dl)
              else
                 grv(:) = get_grv_complex(2*nfld)
              endif
              ! enforce triangular form on correlation matrix
              do m=1,2*nfld; do n=1,2*nfld
                 if (n<m) then
                     corr_fac(n,m) = (0._dl, 0._dl)
                 endif
              enddo; enddo
              ! realize field modes
              call ztrmv('L','N','N', 2*nfld, corr_fac, 2*nfld, grv, 1)
              corr_Fk(:,iii,jjj,kkk) = grv(:)/nvol*norm_in
           endif
        enddo
     enddo
  enddo

  ! initialize homogeneous fields
  call lat_init_homo(f0_in,df0_in,f_init, df_init, zeta_init)
  ! invert FFT, add fluctuations to homogeneous fields
  do m=1,nfld
     Fk(:,:,:) = corr_Fk(m,:,:,:)!Fk(:,:,:) = corr_Fk(2*m-1,:,:,:)
     call fftw_execute_dft_c2r(planb, Fk, f)
     f_init(m,IRANGE) = f_init(m,IRANGE) + f(IRANGE)
     Fk(:,:,:) = corr_Fk(nfld+m,:,:,:)!Fk(:,:,:) = corr_Fk(2*m,:,:,:)
     call fftw_execute_dft_c2r(planb, Fk, f)
     df_init(m,IRANGE) = df_init(m,IRANGE) + f(IRANGE)
  enddo
end subroutine lat_init_spec_convergence_test
#endif

  ! Subroutine to initialize the fields on the lattice to match Fourier space
  ! correlations using spectral initialization. Initializes zeta from metric 
  ! perturbations initialized  consitantly with the field fluctuations.
  ! to do: initialize metric perturbations in line
  ! to do: multiply mp_fp_interp and mp_fms_interp
  ! to do: multiply product by grv(1:2)
  ! to do: initialize zeta fft inline do this by filling Fk
  ! to do: check if offset is required for wavenumber from corr_in
  ! to do: matrix multiplication alpha check data type
!  subroutine lat_init_spec(corr_in, k_cut, seed_in)
!    real(dl) :: corr_in(2*nfld,2*nfld,fp_nos)
!    real(dl), intent(in) :: k_cut
!    integer, intent(in) :: seed_in
!
!    integer :: i, j, k
!    integer :: ii, jj, kk, l
!    real(dl) :: rad
!    integer :: n, m
!
!    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,2*nfld) :: corr_fac
!    complex(C_DOUBLE_COMPLEX), dimension(2*nfld,nnx,ny,nz) :: corr_Fk
!
!    complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_fp_interp
!    complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_fms_interp
!    complex(C_DOUBLE_COMPLEX), dimension(2,2) :: mp_init
!
!    complex(C_DOUBLE_COMPLEX), dimension(2*nfld) :: grv
!    complex(C_DOUBLE_COMPLEX), dimension(2) :: mp
!   
!    ! initialize rng
!    call init_rng(seed_in)    
!
!    ! loop over wavenumbers
!    do k=1,nz; if (k>nnz) then; kk = nz+1-k; else kk=k-1; endif
!      do j=1,ny; if (j>nny) then; jj = ny+1-j; else jj=j-1; endif 
!        do i=1,nnx; ii=i-1
!          rad = sqrt(dble(ii**2)+dble(jj**2)+dble(kk**2))
!          ! initialize corr_fac
!          l = floor(rad*fp_kos)
!          corr_fac(:,:)=(1._dl,0._dl)*((1._dl-l+rad*fp_kos)*corr_in(:,:,l)&
!                                       + (rad*fp_kos-l)*corr_in(:,:,l+1))
!          !do m=1,2*nfld; do n=1,2*nfld
!          !  corr_fac(m,n) = (1._dl,0._dl)*()
!          !enddo; enddo
!          l = floor(rad*mp_kos)
!          mp_fp_interp(:,:)=(1._dl,0._dl)*((1._dl-l+rad*mp_kos)*mp_fp_mat(:,:,l)&
!                                             + (rad*mp_kos-l)*mp_fp_mat(:,:,l+1))
!          mp_fms_interp(:,:)=(1._dl,0._dl)*((1._dl-l+rad*mp_kos)*mp_fms(:,:,l) &
!                                             + (rad*mp_kos-l)*mp_fms(:,:,l+1))
!          ! Multiply mp_fms by mp_fp
!          call zgemm('N','N',2,2,2,(1._dl,0._dl),mp_fms_interp,2, &
!                     mp_fp_interp,2,(0._dl,0._dl),mp_init,2)
!          ! Choleski factorization of corr_fac
!          call zpotrf('L',2*nfld,corr_fac,2*nfld,l)
!          ! print warning if factorization failed
!          if (l /= 0 .and. (i/=1 .or. j/=1 .or. k/=1))
!            print*, "Factorization Warning: l = ", l
!            print*, "ii, jj, kk = ", ii, jj, kk
!            print*, "corr_fac = ", corr_fac
!          endif
!          ! call Gaussian random variables
!          if (i==1 .and. j==1 .and. k==1) then
!            grv(:) = (0._dl, 0._dl)
!          elseif (i==1 .or. i==nx) then
!            grv(:) = (1._dl, 0._dl)*get_grv_real(2*nfld)
!          else
!            grv(:) = get_grv_complex(2*nfld)
!          endif
!          ! realize metric perturpation modes
!          call zgemv('N',2,2,(1._dl,0._dl),mp_init,2,grv(1:2),1, &
!                     (0._dl,0._dl),mp,1)
!          ! enforce triangular form on correlation matrix
!          do m=1,2*nfld; do n=1,2*nfld
!            if (n<m) then
!              corr_fac(n,m) = (0._dl, 0._dl)
!            endif
!          enddo
!          ! realize field modes
!          call ztrmv('L','N','N', 2*nfld, corr_fac, n*nfld, grv, 1)
!          corr_Fk(:,LATIND) = grv(:)/nvol
!          ! realize zeta modes
!          call set_bg_cosmo_ic(a_bg_0,hub_bg_0,f_bg_0,df_bg_0)
!          Fk(LATIND) = (2._dl*(mp(2)/hub_bg_0 + mp(1)) &
!                       /(3._dl*(1._dl+get_bg_w_conf())) + mp(1))/nvol
!          ! cut modes above kcut
!          if (rad .GE. k_cut) then
!            corr_Fk(:,LATIND) = 0._dl*grv(:) ! cut field fluctuations
!            Fk(LATIND) = (0._dl,0._dl) ! cut zeta fluctuations
!          endif
!        enddo
!      enddo
!    enddo

!    ! invert FFT to initialize zeta
!    call fftw_execute_dtf_c2r(planb, Fk, laplace)
!    zeta_lat(:,:,:) = laplace(:,:,:)
!
!    ! initialize homogeneous fields
!    call lat_init_homo()
!    ! invert FFT, add fluctuations to homogeneous fields
!    do m=1,nfld
!      Fk(:,:,:) = corr_Fk(2*m-1,:,:,:)
!      call fftw_execute_dft_c2r(planb, Fk, laplace)
!      fld(m,:,:,:) = fld(m,:,:,:) + laplace(:,:,:)
!      Fk(:,:,:) = corr_Fk(2*m,:,:,:)
!      call fftw_execute_dft_c2r(planb, Fk, laplace)
!      fldp(m,:,:,:) = fldp(m,:,:,:) + laplace(:,:,:)
!    enddo
!    ! set scale factor
!    yscl = 1._dl
!  end subroutine lat_init_spec


  ! Subroutine to initialize the fields on the lattice to match real space
  ! correlations using convolution initialization. Initializes zeta from metric
  ! perturbations initialized consitantly with the field fluctuations.
  ! to do: figure out how to 
  subroutine lat_init_conv()

  end subroutine lat_init_conv

  ! Subroutine to initialize the fields on the lattice to match those supplied in files.
  ! n.b. yscl and ysclp need to be initialized separately
  subroutine lat_init_file(file_in, n_file, f_init, df_init)
    character(len=*) :: file_in(:,:)   ! names of fields and momenta files
    integer :: n_file
    real(dl) :: f_init(:,:,:,:)     ! array of fields to be initialized
    real(dl) :: df_init(:,:,:,:)    ! array of momenta to be initialized
    integer :: i

    do i=1,nfld
       call read_fld(file_in(1,i), n_file, f_init(FLDIND(i)))
       call read_fld(file_in(2,i), n_file, df_init(FLDIND(i)))
    enddo
  end subroutine lat_init_file

  ! Subroutine to read a single field from a file.
  subroutine read_fld(file_in, n_file, f_init)
    character(len=*) :: file_in
    integer :: n_file
    real(dl) :: f_init(:,:,:)
    
    open(n_file, file=file_in, form="unformatted", access="stream",status="old")
    read(n_file) f_init(IRANGE)  ! read in field
    close(n_file)
  end subroutine read_fld

  ! Subroutine to read seed from a file
  subroutine read_seed(file_in, n_file, seed)
    character(len=*) :: file_in
    integer :: n_file
    integer :: seed

    open(n_file, file=file_in,status="old")
    read(n_file, "(I5)") seed  ! read in field
    close(n_file)
  end subroutine read_seed

  subroutine write_filt(n_file, filt, kfilt, ncut)
    integer, intent(in) :: n_file
    procedure(filt_template) :: filt     ! filter function
    real(dl), intent(in) :: kfilt(:)     ! filter parameters
    integer, intent(in) :: ncut

    integer :: i
    
    open(unit=n_file, file='filt.out')
    do i=1,ncut
       write(n_file, '(30(ES22.15, 2x))') i*dk, filt(i*dk,kfilt)
    end do
    close(n_file)
  end subroutine write_filt

end module lat_init_mod
