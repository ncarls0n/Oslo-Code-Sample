module halos_to_field_tools
    implicit none
contains

    ! The Einasto profile 
    ! https://en.wikipedia.org/wiki/Einasto_profile
    real function rho_einasto(r,rho_0,r_0,alpha)
        real, intent(in) :: r, rho_0, r_0, alpha
        rho_einasto = rho_0 * exp( (r/r_0)**alpha )
    end function rho_einasto

    ! The NFW profile
    ! https://en.wikipedia.org/wiki/Navarro–Frenk–White_profile
    real function rho_nfw(r,rho_0,R_s)
        real, intent(in) :: r, rho_0, R_s
        real                u
        u = r/R_s
        rho_nfw = rho_0 / ( u*(1+u)**2 )
    end function rho_nfw

    ! The generalised NFW profile (gNFW)
    ! see WebSky paper [2001.08787] figure 3.4 or [9509122]
    real function rho_gnfw(r,rho_0,r_0,a,b,c)
        real, intent(in) :: r,rho_0,a,b,c
        real                u
        u = r/r_0
        rho_0 * u**c * ( 1+u**a )**-b
    end function rho_gnfw

    subroutine halos_to_field
        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        ! Parameter declarations
        real, parameter :: pi=4.*atan(1.)
        integer, parameter :: subres=4

    ! General use variables
    integer i, N_args
    real, parameter :: pi=4.*atan(1.)

    ! Input variables
    character(len=512), intent(in) :: halo_file, field_file
    integer,            intent(in) :: nbuff
    real,               intent(in) :: s_box, Omega_m0, rho_c

    ! Halo catalogue variables
    integer N_halos, file_size, header_offset, row_offset
    real R_th_max, s_box
    real, allocatable :: x(:), y(:), z(:), R_th(:), M(:)

    ! Variable declarations
    !integer i ! n,i,j,k,ii,jj,kk,l,nbuff,nyquist,klen
    !integer(8) file_size, offset_bytes
    !real s_box, dk, mode, weights(2), R_th_max, z_obs
    !real, allocatable :: kbins(:), P_k(:), P_0(:)


    ! FFTW variable delcarations
    real(c_float),            pointer :: f_x(:,:,:)
    complex(c_float_complex), pointer :: f_k(:,:,:)
    type(c_ptr)                       :: f_ptr, FT_plan, iFT_plan
    integer(c_size_t) fft_size !, total_local_sizes
    integer(c_int)    nfft
    integer           ierr

    ! Get the number of command line arguments
    N_args = command_argument_count()

    ! Read at command line
    read(*,*) halo_file  ! file name for the 3D field
    read(*,*) nbuff      ! buffer thickness (lattice units)
    read(*,*) s_box      ! side length of the 3D field (Mpc)
    if(N_args > 3) then
        read(*,*) Omega_m0   ! Matter energy density fraction at z=0
    else
        Omega_m0 = 0.3138
    endif
    if(N_args > 4) then
        read(*,*) rho_c      ! Critical energy density (Msun/Mpc^3)
    else
        rho_c = 1.2589124e11 ! Msun/Mpc^3
    endif
    if(N_args > 5) then
        read(*,*) field_file ! file name for the 3D field
    else
        field_file = trim(halo_file)//'halo_field.dat'
    endif
    if(N_args > 6) then
        read(*,*) M_bin_min, M_bin_max, M_bin_len
    else
        M_bin_min=-1.
        M_bin_max=-1.
        M_bin_len=20

    ! Read in the header
    open(unit=99,file=halo_file,access='stream')
    read(99) N_halos
    read(99) R_th_max
    read(99) z_obs
    close(99)
    header_offset = 3*4 ! 3 4-byte floats in the header

    ! Allocate halo catalogue arrays
    allocate(   x(N_halos))
    allocate(   y(N_halos))
    allocate(   z(N_halos))
    allocate(R_th(N_halos))

    ! Read file size in bytes, determine resolution n^3 of the 3D field
    inquire(file=trim(halo_file), size=file_size)
    if(file_size==-1) then
        write(*,*) 'Size of file '//trim(halo_file)// &
            ' could not be read. Exiting.'; stop
    endif
    row_offset = (file_size-header_offset)/(4*N_halos)

    ! Read in halo catalogue
    open(unit=99,file=halo_file,access='stream')
    do i=1,N_halos
        read(99,pos=header_offset+(i-1)*row_offset   ) x(    i )
        read(99,pos=header_offset+(i-1)*row_offset+ 4) y(    i )
        read(99,pos=header_offset+(i-1)*row_offset+ 8) z(    i )
        read(99,pos=header_offset+(i-1)*row_offset+24) R_th( i )
    enddo
    M = (4.0/3.0)*pi*R_th**3*Omega_m0*rho_c

    ! Populate 


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

    end subroutine halos_to_field

end module halos_to_field_tools

program main
    use halos_to_field_tools
    implicit none
    call halos_to_field()
end program main

