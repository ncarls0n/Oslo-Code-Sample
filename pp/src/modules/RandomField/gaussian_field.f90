module gaussian_field

  use grid
  use fftw_interface
  use mpivars
  use timing_diagnostics
  use globalvars
  use openmpvars

  implicit none

  !-----------------------------------------------------------------------
  ! BOOKKEEPING AND GLOBAL VARIABLES
  !-----------------------------------------------------------------------

  real xflag,yflag,zflag,rhoflag
  real kx,ky,kz,kmax,ak,dk,d3k,dq,V_k
  real z
  real D

  integer seed

contains
  
  subroutine FT_field

    implicit none
    ! --------------------------------------------------------------------
    ! FFT FROM REAL SPACE TO K SPACE
    ! --------------------------------------------------------------------
    
    if(myid==0.and.rf_report==1) call timer_begin
    
    call fftwf_mpi_execute_dft_r2c(plan, delta, deltac) ! function defined in Intel MPI module on Niagara, plan declared in src/modules/RandomField/fftw_interface.f90
    
    if(myid==0.and.rf_report==1) then
       timing_diagnostics_code='real2complex'
       call timer_end
    endif
    
    return
  end subroutine FT_field

  !=======================================================================

  subroutine iFT_field

    implicit none
    real*8 lavg,avg,lsig,sig
    ! --------------------------------------------------------------------
    ! FFT FROM K SPACE TO REAL SPACE
    ! --------------------------------------------------------------------
    
    if(myid==0.and.rf_report==1) call timer_begin
    
    call fftwf_mpi_execute_dft_c2r(iplan, deltac, delta)! function defined in Intel MPI module on Niagara, iplan decla    red in src/modules/RandomField/fftw_interface.f90
    
    if(myid==0.and.rf_report==1) then
       timing_diagnostics_code='complex2real'
       call timer_end
    endif

    ! Normalization of backwards FFT
    do k=1,local_nz
       do j=1,n
          do i=1,n
             delta(i,j,k) = delta(i,j,k) / real(n)**3
          enddo
       enddo
    enddo
    ! Report mean and variance
    lavg = sum(delta(1:n,:,:)) / real(n)**3
    call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
    lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
    call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
       
    if(myid==0.and.rf_report==1) write(*,96) avg,sqrt(sig)    
96   format(3x,'density:',       /5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)

    return


  end subroutine iFT_field

  !=======================================================================

  subroutine generate_noise

    use random ! module src/modules/RandomField/random.f90

    implicit none

    real*8 lavg,avg,lsig,sig
    real work(1)
    real(dp) :: xr 

    integer :: seeds(MaxRandNumStreams,IRandNumSize)   

    if(myid==0.and.rf_report==1) call timer_begin

    ! Create white noise by assigning pseudorandom numbers to array "seeeds" with dim=(ntasks,4(bytes)) based on a seed value "seed"
    call rans(ntasks,seed,seeds) ! subroutine in random.f90
    do k=1,local_nz !local_nz assigned in fftw_initialize() as mpi slab size, equal to n in runs with ntile=1
       do j=1,n
          do i=1,n
             call gaussdev(seeds(myid+1,:),xr) !myid=0 for peak patch runs
             delta(i,j,k)=xr !xr is a gaussian random deviate with mu=0, sigma=1
          enddo
       enddo
    enddo

    ! Report mean "avg" and variance "sig"
    lavg = sum(delta(1:n,:,:)) / real(n)**3
    call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr) !calls mpi to sum lavg components from each tile
    lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
    call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr) !calls mpi to sum lsig components from each tile

    if(myid==0.and.rf_report==1) then
       write(*,96) avg,sqrt(sig)
       timing_diagnostics_code='white noise'
       call timer_end
    endif

    return
96  format(3x,'white noise:',       /5x' mean = ',&
           1pe13.6,/5x,'sigma = ',1pe13.6)
    
  end subroutine generate_noise
  
  !=======================================================================
  
  subroutine convolve_noise

    use pktable

    implicit none    

    complex(C_FLOAT_COMPLEX) ctemp

    real, allocatable :: work(:)

    double precision lavg,avg,lsig,sig

    real tf     ! zeta transfer function
    real pk     ! CMB power spectrum
    real pkchi  ! Power preheating model chi
    real pkchix ! Power spectrum for chi_extreme
!    pksav=sqrt(pksav_in*d3k*n**3)
!    tsav=sqrt(tsav_in*d3k*n**3)
    ! --------------------------------------------------------------------
    ! FFT FROM REAL SPACE TO K SPACE
    ! --------------------------------------------------------------------

    if((code=='rho').or.(code=='zeta2delta').or.(code=='zetag').or.&
       (code=='chi').or.(code=='chi_extreme')) then
       if(myid==0.and.rf_report==1) call timer_begin

       !Assign `deltac` as complex-valued DFT of `delta`
       call fftwf_mpi_execute_dft_r2c(plan, delta, deltac)

       if(myid==0.and.rf_report==1) then
          timing_diagnostics_code='real2complex'
          call timer_end
       endif
    endif

    ! --------------------------------------------------------------------
    ! MULTIPLY delta(k1,k2,k3) by P(k)                                    
    ! --------------------------------------------------------------------
    
    if(myid==0.and.rf_report==1) call timer_begin

!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(deltac) 
    do k=1,local_nz!local_nz assigned in fftw_initialize() as mpi slab size, equal to n in runs with ntile=1
       kz=(k+local_z_start-1)*dk ! local_z_start assigned in fftw_initialize() in module fftw_interface.f90,
       if (k+local_z_start.gt.n12) kz=kz-kmax        ! index of slab start, equal to 0 for runs with ntile=1
       
       do j=1,n
          ky=(j-1)*dk
          if (j.gt.n12) ky=ky-kmax
          
          do i=1,n12+1
             kx=(i-1)*dk
             if (i.gt.n12) kx=kx-kmax
             
             ak=sqrt(kx**2+ky**2+kz**2)
           
             ! Integer values approximately equating k in spectra read in to elements in deltac(i,j,k)
             itab=int((log10(ak)-pktabminl)/dpktab)+1 ! rounds up to the nearest integer element itab of the power spectrum array pksav so that pksav(itab-1) <= P(ak) < pksav(itab)
             if(ak==0.) then
                pk    = 0.
                pkchi = 0.
                pkchix= 0.
                tf    = 0.
             !if(ak==0.0 .or. ak>n12*dk) then ! putting in tophat cut
             !   pk    = 0.
             !   pkchi = 0.
             !   tf    = 0.
             else
                dtab=(log10(ak)-(pktabminl+(itab-1)*dpktab))/dpktab
                !Interpolates values of:
                pk    = pksav(itab)*(1.-dtab) + pksav(itab+1)*dtab !square root of linear power spectrum sqrt( P(k) )
                tf    = tsav(itab) *(1.-dtab) + tsav(itab+1) *dtab !CMB transfer function sqrt( k^2 T(k) )
                pkchi = pkchisav(itab)*(1.-dtab)+pkchisav(itab+1)*dtab !chi power spectrum sqrt( P_chi(k) ) from preheating models
                pkchix= pkchixsav(itab)*(1.-dtab)+pkchixsav(itab+1)*dtab ! power spectrum for inflationary chi_extreme in intermittent, spatially-localized non-gaussianity model
                ! from the existing tabulated values read in as pksav, tsav, pkchisav
                ! by subroutine read_pktable in RandomField_init
                
             endif
             ctemp = deltac(i,j,k)
             
             if(code.eq.'rho') then
                ctemp = ctemp * pk ! delta(k) = white_noise(k) * sqrt( P(k)*(2 pi n/boxsize)^3 )
             elseif(code.eq.'chi') then
                ctemp = ctemp * pkchi
             elseif(code.eq.'zetag') then
                ctemp = ctemp * pk/tf**2 ! zeta_gaussian(k) = white_noise(k) * sqrt( P(k)*(2 pi n/boxsize)^3 ) / ( k^2 T(k) )
                !if(tf/=0) ctemp = ctemp * pk/tf**2 ! putting in tophat cut
                !if(tf==0) ctemp = ctemp * 0.0
             elseif(code.eq.'zeta2delta') then
                ctemp = ctemp * tf**2 
             elseif(code=='chi_extreme') then
                ctemp = ctemp * pkchix
             elseif(code.eq.'xlpt1') then
                ctemp = ctemp * kx / ak**2 * cmplx(0.0,1.0)
                if(i.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'ylpt1') then
                ctemp = ctemp * ky / ak**2 * cmplx(0.0,1.0)
                if(j.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'zlpt1') then
                ctemp = ctemp * kz / ak**2 * cmplx(0.0,1.0)
                if(k+local_z_start.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'lapld') then
                ctemp = ctemp * ak**2 

             !2LPT added. First get phi,ij
             elseif(code.eq.'phi11') then
                ctemp = -ctemp * kx * kx / ak**2 
                if(i.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'phi22') then
                ctemp = -ctemp * ky * ky / ak**2 
                if(j.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'phi33') then
                ctemp = -ctemp * kz * kz / ak**2 
                if(k+local_z_start.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'phi21') then
                ctemp = -ctemp * ky * kx / ak**2 
                if(i.eq.n12+1) ctemp = 0.0
                if(j.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'phi31') then
                ctemp = -ctemp * kz * kx / ak**2 
                if(i.eq.n12+1) ctemp = 0.0
                if(k+local_z_start.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'phi32') then
                ctemp = -ctemp * kz * ky / ak**2 
                if(j.eq.n12+1) ctemp = 0.0
                if(k+local_z_start.eq.n12+1) ctemp = 0.0
             !For Psi^(2)
             elseif(code.eq.'xlpt2') then
                ctemp = -ctemp * kx / ak**2 * cmplx(0.0,1.0)
                if(i.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'ylpt2') then
                ctemp = -ctemp * ky / ak**2 * cmplx(0.0,1.0)
                if(j.eq.n12+1) ctemp = 0.0
             elseif(code.eq.'zlpt2') then
                ctemp = -ctemp * kz / ak**2 * cmplx(0.0,1.0)
                if(k+local_z_start.eq.n12+1) ctemp = 0.0
             endif
             
             if(ak==0.) ctemp = 0.0

             deltac(i,j,k) = ctemp
             
          enddo
       enddo
    enddo
    
    if(myid==0.and.rf_report==1) then
       timing_diagnostics_code='convolution'
       call timer_end
    endif

    ! --------------------------------------------------------------------
    ! FFT FROM K SPACE TO REAL SPACE
    ! --------------------------------------------------------------------

    if((code.ne.'rho') .and. (code.ne.'zeta2delta')) then
       if(myid==0.and.rf_report==1) call timer_begin
    
       call fftwf_mpi_execute_dft_c2r(iplan,deltac,delta) ! delta=iFFT(deltac)

       if(myid==0.and.rf_report==1) then
          timing_diagnostics_code='complex2real'
          call timer_end
       endif

       ! Normalization of backwards FFT
       do k=1,local_nz
          do j=1,n
             do i=1,n
                delta(i,j,k) = delta(i,j,k) / real(n)**3
             enddo
          enddo
       enddo
       ! Report mean and variance
       lavg = sum(delta(1:n,:,:)) / real(n)**3
       call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                          mpi_comm_world,ierr)
       lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
       call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
                          mpi_comm_world,ierr)
       
       if(myid==0.and.rf_report==1) then
          if(code.eq.'xlpt1')       write(*,97) avg,sqrt(sig)
          if(code.eq.'ylpt1')       write(*,98) avg,sqrt(sig)
          if(code.eq.'zlpt1')       write(*,99) avg,sqrt(sig)
          if(code.eq.'zetag')       write(*,100) avg,sqrt(sig)
          if(code.eq.'zeta2delta')  write(*,101) avg,sqrt(sig)
          if(code.eq.'chi')         write(*,102) avg,sqrt(sig)
          if(code=='chi_extreme')   write(*,113) avg,sqrt(sig)
          if(code.eq.'lapld')       write(*,103) avg,sqrt(sig)
          if(code.eq.'phi11')       write(*,104) 
          if(code.eq.'phi22')       write(*,105) 
          if(code.eq.'phi33')       write(*,106) 
          if(code.eq.'phi21')       write(*,107) 
          if(code.eq.'phi31')       write(*,108) 
          if(code.eq.'phi32')       write(*,109) 
          if(code.eq.'xlpt2')       write(*,110) avg,sqrt(sig)
          if(code.eq.'ylpt2')       write(*,111) avg,sqrt(sig)
          if(code.eq.'zlpt2')       write(*,112) avg,sqrt(sig)
       endif
    endif

97   format(3x,'x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
98   format(3x,'y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
99   format(3x,'z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
100  format(3x,'zeta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
101  format(3x,'zeta2delta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
102  format(3x,'chi:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
113  format(3x,'chi_extreme:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
103  format(3x,'density laplacian:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
104  format(3x,'phi,11:')
105  format(3x,'phi,22:')
106  format(3x,'phi,33:')
107  format(3x,'phi,21:')
108  format(3x,'phi,31:')
109  format(3x,'phi,32:')
110  format(3x,'2lpt x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
111  format(3x,'2lpt y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
112  format(3x,'2lpt z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)

    call mpi_barrier(mpi_comm_world,ierr)
    
    return
    
  end subroutine convolve_noise

end module gaussian_field
