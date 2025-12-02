module loc
    implicit none
    contains
    integer function findloc_sorted(array, querry) result(indx)
        ! Assumes "array" to be an array of integers sorted from least to greatest, finds the index i of the first instance of "target" in the array and returns i
        implicit none
        integer, intent(in) :: array(:)
        integer, intent(in) :: querry
        integer i_high,i_low,high,low
        logical flag
        i_low  = 1
        i_high = size(array)
        low  = array(i_low )
        high = array(i_high)
        ! If the querry value is not in the array
        if(querry<low.or.querry>high) then
            indx = -1
            flag = .true.
        ! Check boundaries
        elseif(low==querry) then
            indx = i_low
            flag = .true.
        elseif(high==querry) then
            indx = i_high
            flag = .true.
        endif
        ! Loop to search if querry within boundaries
        do while (flag.eqv..false.)
            ! This loop assumes array(high) and array(low) not equal to querry
            if(i_high<=i_low+1) then
                indx = -1
                flag = .true.
            else
                ! Starting guess assuming linear distribution
                indx = (querry-low)*(i_high-i_low)/(high-low)+i_low
                if(indx==i_high) then
                    indx=indx-1
                elseif(indx==i_low) then
                    indx=indx+1
                endif
                ! Check guess
                if(array(indx)==querry) then
                    flag = .true.
                ! If guess is too high
                elseif(array(indx)>querry) then
                    i_high = indx
                    high   = array(indx)
                ! If guess is too low
                elseif(array(indx)<querry) then
                    i_low = indx
                    low   = array(indx)
                endif
            endif
        enddo
        ! To find the first instance of querry (in case there are doubles)
        do while (flag.eqv..true.)
            if(indx>1.and.array(indx-1)==querry) then
                indx=indx-1
            else
                flag=.false.
            endif
        enddo
    end function findloc_sorted

    integer function findloc_dumb(array,querry) result(indx)
        implicit none
        integer, intent(in) :: array(:)
        integer, intent(in) :: querry
        logical flag
        flag=.false.
        indx=1
        do while((flag.eqv..false.).and.(indx<=size(array)))
            if(array(indx)==querry) then
                flag=.true.
            else
                indx=indx+1
            endif
        enddo
        if(indx>size(array)) indx=-1
    end function findloc_dumb

end module loc

program get_powerspectrum
    use, intrinsic :: iso_c_binding
    use loc
    implicit none
    include 'fftw3f-mpi.f03'
    real, parameter :: pi=4.*atan(1.)
    character(len=128) field_file
    integer file_size, n, i,j,k,l, cells, radii, start, finish
    real s_box
    integer, allocatable :: radial(:,:,:), key(:,:), counts(:), temp(:)
    logical, allocatable :: mask(:)
    real,    allocatable :: wavenum(:), P(:)

    real(c_float),            pointer :: f_x(:,:,:)
    complex(c_float_complex), pointer :: f_k(:,:,:)
    type(c_ptr)                       :: f_ptr, FT_plan, iFT_plan
    integer(c_size_t)                    fft_size !, total_local_sizes
    integer(c_int)                       nfft
    integer                              ierr
!    ! Global
!    real(c_float),            pointer :: f_x_global(:,:,:)
!    complex(c_float_complex), pointer :: f_k_global(:,:,:)
!    type(c_ptr)                       :: f_global_ptr

!integer q(5,5), thing(5)
!thing = (/ 0,1,2,3,4 /)
!do i=-1,5
!write(*,*) i, findloc_sorted( thing,i ), findloc_dumb( thing,i )
!enddo
!goto 1

    read(*,*) field_file

    ! read size of file in bytes
    inquire(file=trim(field_file), size=file_size)
    
    ! check that file size was read correctly and determine resolution n^3 of the box
    if(file_size==-1) then
        write(*,*) 'INQUIRE returned -1, file '//trim(field_file)//' could not be read. Exiting.'
        goto 1
    else
        n = (file_size/4)**(1./3.)
        if(4*n**3/=file_size) then
            write(*,*) 'File must have n^3 elements where n is an integer, exiting.'
            goto 1
        endif
    endif

    ! read s_box
    ! field files take the form *_###Mpc_n###_nb##_nt##
    do i=126,1,-1
        if(field_file(i:i)=='M'.and.field_file(i+1:i+1)=='p'.and.&
           field_file(i+2:i+2)=='c') finish=i-1
    enddo
    do i=finish,1,-1
        if(field_file(i:i)=='_') then
            start=i+1
            goto 2
        endif
    enddo
    2 continue

    ! Convert string to real
    read(field_file(start:finish),*) s_box

    ! Count the number of cells within the radius 2*pi/s_box of cell i,j,k=0,0,0 (in Fourier space) and set that to "cells" to allocate arrays
    cells=0
    do i=1,n
        do j=1,min(i,int((n**2-i**2)**.5))
            do k=1,min(j,int((n**2-i**2-j**2)**.5),n/2+1)
                cells=cells+1
            enddo
        enddo
    enddo
    allocate( counts(cells) ) ! number of cells with a given radius
    allocate(   temp(cells) ) 
    allocate(   mask(cells) ) ! logical array for sorting counts(:)

    ! Sort "counts" from least to greatest radius and save it as "temp"
    write(*,*) 'Making k field'
    cells = 0
    mask  = .true.
    do i=1,n
        do j=1,min(i,int((n**2-i**2)**.5))
            do k=1,min(j,int((n**2-i**2-j**2)**.5),n/2+1)
                cells=cells+1
                counts(cells) = (i-1)**2 + (j-1)**2 + (k-1)**2
            enddo
        enddo
    enddo
    write(*,*) 'Sorting by k'
    do i=1,cells
        j=minval(counts,dim=1,mask=mask)
        k=minloc(counts,dim=1,mask=mask)
        temp(i)=j ! minval(counts,mask)
        mask(k)=.false. !minloc(counts,mask))=.false.
    enddo
    deallocate(mask)

    write(*,*) 'Binning by k'
    radii=1
    counts=1
    do i=1,cells-1
        if(temp(i)==temp(i+1)) then
            counts(radii)=counts(radii)+1
        else
            radii=radii+1
        endif
    enddo
    deallocate(counts)

    ! Make a key to go back and forth between P(k) and bins
    write(*,*) 'making key'

    allocate(key(radii,2))
    key(:,2)=1
    radii=1
    do i=1,cells-1
        if(temp(i)==temp(i+1)) then
            key(radii,2)=key(radii,2)+1
        else
            key(radii,1)=temp(i)
            radii=radii+1
        endif
    enddo

    deallocate(temp)
    write(*,*) 'making key but actually this time'
    allocate( radial(n,n,n/2+1)  )
    do i=1,n
        do j=1,min(i,int((n**2-i**2)**.5))
            do k=1,min(j,int((n**2-i**2-j**2)**.5),n/2+1)
                radial(i,j,k)=findloc_dumb( key(:,1), (i-1)**2+(j-1)**2+(k-1)**2)
                if(i/=j)              radial(j,i,k)=radial(i,j,k)
                if(i/=k.and.i<=n/2+1) radial(k,j,i)=radial(i,j,k)
                if(j/=k.and.j<=n/2+1) radial(i,k,j)=radial(i,j,k)
                if(i/=j.and.i/=k.and.j/=k) then
                    if(j<=n/2+1) radial(k,i,j)=radial(i,j,k)
                    if(i<=n/2+1) radial(j,k,i)=radial(i,j,k)
                endif
            enddo
        enddo
    enddo

    ! ! Initialize parallel FFTs
    ! ierr = fftwf_init_threads()
    ! call fftwf_plan_with_nthreads(omp_get_max_threads())
    ! total_local_sizes = fftwf_mpi_local_size_3d(nmpifft,nmpifft,nmpifft,mpi_com_world,local_nz,local_z_start)
    ! mpiFT_plan  = fftwf_mpi_plan_dft_r2c_3d(nmpifft,nmpifft,nmpifft,&
    !     f_x_global,f_k_global,mpi_comm_world,fftw_estimate)
    ! mpiiFT_plan = fftwf_mpi_plan_dft_c2r_3d(nmpifft,nmpifft,nmpifft,&
    !     f_k_global,f_x_global,mpi_comm_world,fftw_estimate)

    ! Initialize serial FFTs
    FT_plan  = fftwf_plan_dft_r2c_3d(nfft,nfft,nfft,f_x,f_k,fftw_estimate)
    iFT_plan = fftwf_plan_dft_c2r_3d(nfft,nfft,nfft,f_k,f_x,fftw_estimate)
!
!    ! Allocate serial FFTs
!    nfft     = n
!    fft_size = (nfft+2)*nfft**2
!    f_ptr    = fftwf_alloc_real(fft_size)
!    call c_f_pointer( f_ptr, f_x, [2*(nfft/2+1),nfft,nfft] )
!    call c_f_pointer( f_ptr, f_k, [   nfft/2+1 ,nfft,nfft] )

    ! ! Allocate parallel FFTs
    ! fft_size = 2*(nmpifft/2+1)*nmpifft*local_nz
    ! f_ptr_global = fftf_alloc_real(fft_size)
    ! call c_f_pointer(f_ptr_global, f_x_global, [2*(nmpifft/2+1),nmpifft,local_nz])
    ! call c_f_pointer(f_ptr_global, f_k_global, [   nmpifft/2+1 ,nmpifft,local_nz])

    ! For a Fourier-space field f_k(i,j,k) the wavenumber corresponding to this point is key(radial(i,j,k),1) and the power spectrum is P(radial(i,j,k))
!    allocate(wavenum(radii))
!    allocate(counts( radii))
!    allocate(P(      radii))
!    wavenum = key(:,1)*2.*pi/s_box
!    counts  = key(:,2)
!    deallocate(key)
!    P       = 0.0

    ! Read in field and FFT to get Fourier-space field
!    open(unit=99,file=field_file,access='stream')
!    read(99,pos=0) (((f_x(i,j,k),i=1,n),j=1,n),k=1,n)
!    call fftwf_mpi_execute_dft_r2c(FT_plan, f_x, f_k)
!    do i=1,n
!        do j=1,n
!            do k=1,n//2+1
!                !wavenum(radial(i,j,k)) = key(radial(i,j,k),1) * 2.0*pi/s_box
!                !wavenum(radial(i,j,k)) = real((i**2+j**2+k**2)**.5) * 2.0*pi/s_box
!                !P(radial(i,j,k)) = P(radial(i,j,k))+field(i,j,k)/key(radial(i,j,k),2)
!                P(radial(i,j,k)) = P(radial(i,j,k)) + f_k(i,j,k)/counts(radial(i,j,k))
!            enddo
!        enddo
!    enddo

    ! ! iFFT to get f_x
    ! call fftwf_mpi_execute_dft_c2r(iFT_plan, f_k, f_x)
    ! do k=1,n
    !     do j=1,n
    !         do i=1,n
    !             f_x(i,j,k) = f_x(i,j,k) / real(n)**3
    !         enddo
    !     enddo
    ! enddo
    
    1 continue
end program get_powerspectrum

