program get_partial_power
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    ! Parameter declarations
    real, parameter :: pi=4.*atan(1.)
    integer, parameter :: subres=4

    ! Variable declarations
    character(len=128) field_file, out_file
    integer n,i,j,k,ii,jj,kk,l,nbuff,nyquist,klen
    integer(8) file_size, offset_bytes
    real s_box, dk, mode, weights(2)
    real, allocatable :: kbins(:), P_k(:), P_0(:)

    ! FFTW variable delcarations
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
    out_file = trim(field_file)//'_partial_power.dat'

    ! Read file size in bytes, determine resolution n^3 of the 3D field
    inquire(file=trim(field_file), size=file_size)
    if(file_size==-1) then
        write(*,*) 'Size of file '//trim(field_file)// &
            ' could not be read. Exiting.'; stop
    else
        n = (file_size/4)**(1./3.)
        if(4*n**3/=file_size) then
            write(*,*) 'File must have n^3 elements where n is an integer. Exiting.'
            stop
        endif
    endif

    ! We don't need to sample the full power spectrum for measuring sigma_8 because the tails are exponentially suppressed. The important feature is the peak, so for a large field with more k modes for k>k_peak than k modes for k<k_peak, we can drop some of the high-k modes by only sampling every 4 points in the lattice.
    n=n/subres

    ! Field variables
    dk      = 2.*pi/s_box ! Fourier space resolution (Mpc^-1)
    nyquist = n/2+1       ! the Nyquist number, an integer

    ! Distance (lattice units) to furthest cell when you fold at Nyquist
    klen = sqrt(3.)*nyquist

    ! Allocate wavenumber and power spectrum arrays
    allocate(kbins(klen))
    allocate(P_k(klen))
    allocate(P_0(klen))

    ! Assign wavenumber bins
    kbins(1:klen) = (/ (i*dk, i=0,klen-1) /)

    ! Initialize serial FFTW plans for 3D discrete Fourier transforms
    nfft     = n
    FT_plan  = fftwf_plan_dft_r2c_3d(nfft,nfft,nfft,f_x,f_k,fftw_estimate)
    iFT_plan = fftwf_plan_dft_c2r_3d(nfft,nfft,nfft,f_k,f_x,fftw_estimate)

    ! Allocate serial FFTW arrays and setup pointers
    fft_size = 2*(nfft/2+1)*nfft**2
    f_ptr    = fftwf_alloc_real(fft_size)
    call c_f_pointer( f_ptr, f_x, [2*(nfft/2+1),nfft,nfft] )
    call c_f_pointer( f_ptr, f_k, [   nfft/2+1 ,nfft,nfft] )

    ! Read real-space feild into FFTW array variable
    open(unit=99,file=field_file,access='stream')
    do i=1,n
        do j=1,n
            do k=1,n
                offset_bytes = 1 + int(4,8)*( int(subres*(i-1),8) + &
                                          int(subres**2*n*(j-1),8) + &
                                          int(subres**3*n**2*(k-1),8) )
                read(99,pos=offset_bytes) f_x(i,j,k)
            enddo
        enddo
    enddo
    close(99)

    ! Calculate 3D DFT of field F(f_x)=f_k
    call fftwf_execute_dft_r2c(FT_plan, f_x, f_k)

    ! Calculate power spectrum from Fourier space field f_k
    do i=1,n; if(i<=nyquist)then; ii=i-1; else; ii=n+1-i ;endif
        do j=1,n; if(j<=nyquist)then; jj=j-1; else; jj=n+1-j ;endif
            do k=1,n; if(k<=nyquist)then; kk=k-1; else; kk=n+1-k ;endif

                mode       = sqrt(float(ii**2+jj**2+kk**2))+1
                l          = floor(mode)
                weights    = (1.-(/l-mode,l-mode+1./)**2.)**2.
                P_k(l:l+1) = P_k(l:l+1) + weights * abs(f_k(ii+1,j,k))**2.
                P_0(l:l+1) = P_0(l:l+1) + weights

            enddo
        enddo
    enddo
    do i=1,klen
        if(P_0(i)/=0.) P_k(i) = P_k(i) / P_0(i)
    enddo

    ! Express modulus square of DFT Fourier modes P_k as samples of the of
    ! the power spectrum P(k) of a continuous field. Because of our choice
    ! of conventions, we have the relation P(k) = n^-3 (s_latt/2pi)^3 P_k
    P_k = (s_box/(2.*pi*float(n*(n-2*nbuff/subres))))**3. * P_k

    ! Write out fundamental mode
    write(*,*)
    write(*,'(a44,es10.4e2,a8)') &
        'The 0th Fourier mode is P(k=0)=<|f(k=0)|^2>=',P_k(1),' Mpc^-1.'
    write(*,*)

    ! Save power spectrum in the form used by Peak Patch initial conditions
    open(unit=42,file=out_file)
    do i=2,klen
        write(42,'(e10.5e2,a,e10.5e2)') kbins(i),' ',P_k(i)
    enddo
    close(42)
    write(*,'(a23)') &
        'P(k>0) written to file:'
    write(*,*) trim(out_file)
    write(*,*)

end program get_partial_power

