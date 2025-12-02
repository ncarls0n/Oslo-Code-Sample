program get_powerspectrum
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3-mpi.f03'
    real, parameter :: pi=4.*atan(1.)

    ! Variable declarations
    character(len=128) field_file
    integer i,j,k,ii,jj,kk,nbuff,nyquist, klen, i_low
    integer(8) file_size,n!, i,j,k, ii,jj,kk, nbuff,nyquist, klen, i_low
    real s_box, dk, mode
    real weights(2)
    real, allocatable :: kbins(:), P_k(:), P_0(:)

    ! FFTW variables
    real(c_float),            pointer :: f_x(:,:,:)
    complex(c_float_complex), pointer :: f_k(:,:,:)
    type(c_ptr)                       :: f_ptr, FT_plan, iFT_plan
    integer(c_size_t) fft_size !, total_local_sizes
    integer(c_int)    nfft
    integer           ierr

    ! Read at command line
    read(*,*) field_file ! file name for the 3D field
    read(*,*) nbuff      ! buffer thickness (lattice units)
    read(*,*) s_box      ! side length of the 3D field (Mpc)

    ! read size of file in bytes and determine resolution n^3 of the 3D field
    inquire(file=trim(field_file), size=file_size)
    if(file_size==-1) then
        write(*,*) 'INQUIRE returned -1, file '//trim(field_file)//' could not be read, exiting.'
        stop
    else
        n = (file_size/4)**(1./3.)
        if(4*n**3/=file_size) then
            write(*,*) 'File must have n^3 elements where n is an integer, exiting.'
            stop
        endif
    endif

    ! Field variables
    dk      = 2*pi/s_box
    nyquist = n/2+1 ! since nyquist is an integer, this is like int(n/2)+1

    ! Distance (lattice units) to furthest cell when you fold at Nyquist
    klen = (3.**.5)*(nyquist-1)

    ! Allocate wavenumber and power spectrum arrays
    allocate(kbins(klen+1))
    allocate(P_k(klen+1))
    allocate(P_0(klen+1))

    ! Assign wavenumber bins
    kbins(1:klen+1) = (/ (i*dk, i=1,klen+1) /)

    ! Initialize serial FFTW
    FT_plan  = fftwf_plan_dft_r2c_3d(nfft,nfft,nfft,f_x,f_k,fftw_estimate)
    iFT_plan = fftwf_plan_dft_c2r_3d(nfft,nfft,nfft,f_k,f_x,fftw_estimate)

    ! Allocate serial FFTW arrays and setup pointers
    nfft     = n
    fft_size = 2*(nfft/2+1)*nfft**2
    f_ptr    = fftwf_alloc_real(fft_size)
    call c_f_pointer( f_ptr, f_x, [2*(nfft/2+1),nfft,nfft] )
    call c_f_pointer( f_ptr, f_k, [   nfft/2+1 ,nfft,nfft] )

    open(unit=99,file=field_file,access='stream')
    read(99,pos=1) (((f_x(i,j,k),i=1,nfft),j=1,nfft),k=1,nfft)
    call fftwf_mpi_execute_dft_r2c(FT_plan, f_x, f_k)
    !call fftwf_execute_dft_r2c(FT_plan, f_x, f_k)

!    do i=1,n
!        if(i<nyquist)then; ii=i
!        else;              ii=nyquist-i ;endif
!        do j=1,n
!            if(j<nyquist)then; jj=j
!            else;              jj=nyquist-1 ;endif
!            do k=1,n
!                if(k<nyquist)then; kk=k
!                else;              kk=nyquist-k ;endif
!
!                mode    = sqrt(ii**2.0+jj**2.0+kk**2.0)
!                i_low   = sqrt(ii**2  +jj**2  +kk**2  )
!                weights = (1.-(/i_low-mode,i_low-mode+1/)**2.)**2.
!                P_k(i_low:i_low+1) = P_k(i_low:i_low+1) + weights*abs(f_k(i,j,kk))**2.
!                P_0(i_low:i_low+1) = P_0(i_low:i_low+1) + weights
!
!            enddo
!        enddo
!    enddo
!
!    do i=1,klen+1
!        if(P_0(i)/=0.:
!            P_k(i) = P_k(i) / P_0(i)
!        endif
!    enddo
!    P_k = (n/s_box)**3.*P_k
!
!    ! ! iFFT to get f_x
!    ! call fftwf_mpi_execute_dft_c2r(iFT_plan, f_k, f_x)
!    ! do k=1,n
!    !     do j=1,n
!    !         do i=1,n
!    !             f_x(i,j,k) = f_x(i,j,k) / real(n)**3
!    !         enddo
!    !     enddo
!    ! enddo
!    
end program get_powerspectrum

