program make_zeta_runner
    use, intrinsic :: iso_c_binding
    use intreal_types      ! from src/modules/External, defines integer and real types used in peak patch
    use arrays             ! from src/hpkvd, subroutines: allocate_halos, allocate_boxes
    use params             ! from src/modules/GlobalVariables
    use cosmoparams        ! from src/modules/GlobalVariables, declares cosmological variables
    use input_parameters   ! from src/modules/GlobalVariables, declares input parameters
!   use cosmology          ! from src/modules/RandomField
!   use textlib            ! from src/modules/External
!   use random             ! from src/modules/RandomField
    use globalvars         ! from src/modules/RandomField
!   use pktable            ! from src/modules/RandomField
    use grid               ! from src/modules/RandomField
    use mpivars            ! from src/modules/External
    use timing_diagnostics ! from src/modules/External
    use TabInterp          ! from src/modules/TabInterp/TabInterp.f90
!   use openmpvars         ! from src/modules/External
    use fftw_interface     ! from src/modules/RandomField
    use gaussian_field     ! from src/modules/RandomField
!   use tiles              ! from src/modules/RandomField
!   use chi2zeta           ! from src/modules/RandomField
    use RandomField        ! from src/modules/RandomField
    implicit none
    real, parameter :: pi = 4.0*atan(1.0)
    character(len=512) run_dir,hpkvd_params!,fielddir
    integer            nmesh

    ! Reads in file name for hpkvd_inputs from command line
    read(*,*) run_dir

    ! Reads in parameters from the text file hpkvd_params.txt (which is
    ! just a fortran-readable reformatting of the parameter file
    ! param.params that was made by the python runner make_zeta.py)
    hpkvd_params = trim(run_dir)//'hpkvd_params.txt'
    call read_parameters

    ! Initialize MPI process
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,ntasks,ierr)

    ! Run Peak Patch initial conditions
    call make_zeta

    ! Exit MPI process
    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_finalize(ierr)

contains

    subroutine make_zeta
        ! Random Field variables
        character(len=512) rf_filepktab,rf_outcode,rf_incode,rf_fielddir
        real               rf_boxsize,rf_ainit
        integer            rf_nmesh,rf_nbuff,rf_ntile,rf_seed,rf_ireadfield
        ! Array variables
        real(c_float),            pointer :: delta_(:,:,:),  deltag_(:,:,:),&
                                             delta_w(:,:,:), deltag_w(:,:,:)
        complex(c_float_complex), pointer :: deltac_(:,:,:), deltagc_(:,:,:),&
                                             deltac_w(:,:,:),deltagc_w(:,:,:)
        type(c_ptr)                       :: deltap_,        deltagp_,&
                                             deltap_w,       deltagp_w
        ! Smoothing parameters
        integer half_n1, n_nyquist
        real    k_min,   f_nyquist
        real(c_float) sigma_8 !,sigma_8_serial 
        !complex speq(:,:)

        ! Initialize Random Field
        rf_filepktab  = pkfile ! file w/ tabulated powerspectra & transfer functions
        rf_outcode    = filein ! file for field output
        rf_incode     = densfilein ! file for field input
        rf_fielddir   = trim(run_dir)//trim(fielddir) ! director of <...>/fields/
        rf_boxsize    = dcore_box*(nlx-1) + dL_box 
        ! Note here that dcore_box is the sidelength of a simulation tile
        ! in Mpc/h excluding buffers, and dL_box is the sidelength of a
        ! simulation tile in Mpc/h including buffers, so rf_boxsize is the
        ! sidelength of the whole simulation including exterior buffers
        ! (i.e. ntile*dcore_box+dL_box-dcore_box). This is DIFFERENT from
        ! the definition of "boxsize" used in the python script parameter
        ! file, where it is defined as the sidelength of the whole
        ! simulaiton excluding exterior buffers (i.e. ntile*dcore_box).

        ! For a parallel run divided up into nlx*nly*nlz boxes (which we
        ! call "tiles", and for a cubic run nlx=nly=nlz=ntile), dcore_box
        ! is the sidelength in Mpc/h of a tile excluding buffers, and
        ! dL_box is the sidelength in Mpc/h of a tile including buffers, so
        ! rf_boxsize is the total sidelength of the simmulation volume
        ! including the external buffers only.
        rf_nmesh      = nmesh ! simulation box side length in array elements
        rf_nbuff      = nbuff ! simulation box buffer thickness in array elements
        rf_ntile      = nlx ! simulation box side length in parallelization tiles
        rf_ainit      = 1.0 ! 
        rf_seed       = iseedFFT ! seed number in random field making
        rf_ireadfield = ireadfield ! 0 if make field, 1 if read field instead

        ! subroutine in src/modules/RandomField/RandomField.f90
        myid=1
        call RandomField_init(rf_filepktab,rf_ireadfield,rf_outcode,rf_incode,&
            rf_fielddir,rf_boxsize,rf_nmesh,rf_nbuff,rf_ntile,rf_ainit,&
            NonGauss,fNL,A_nG,B_nG,R_nG)
        myid=0

        ! Allocates arrays delta and deltac respectively as real and k-space
        ! Fourier conjugates
!        call fftw_allocate_serial(   delta_w, deltac_w, deltap_w )
        call fftw_allocate_parallel( deltag_w,deltagc_w,deltagp_w)
!        call fftw_allocate_serial(   delta_,  deltac_,  deltap_  )
        call fftw_allocate_parallel( deltag_, deltagc_, deltagp_ )
        ! subroutines in src/modules/RandomField/fftw_interface.f90

        ! Assigns the seed number used in subsequent RandomField_make calls
        call RandomField_seed(rf_seed)
        call RandomField_make(-1, deltag_, deltagc_, deltagc_)
        ! subroutines in src/modules/RandomField/RandomField.f90
        ! -1 is the ID for making an overdensity field, so deltag_ is
        ! populated with Gaussian random deviate, deltagc_ is set to the FT
        ! deltag_, convolved with the powerspectrum and approprate transfer
        ! functions for the type of NonGaussianity used

!        ! Smooth k-space field with real-space top-hat kernel with
!        !  R_TH = 8 Mpc/h 
!        k_min = 2*pi/dL_box ! Smallest wavenumber allowed by tile of this size
!        call smooth_field(deltagc_,deltagc_w,k_min,8.0)
!
!        ! Calculate standard deviation (sigma_8 since R_TH = 8 Mpc/h)
!        sigma_8=sqrt( sum(deltag_w(:,:,:)**2) / (n1*nlx*n2*nly*n3*nlz) )
!        sigma_8_serial=sqrt(sum(delta_w(1:n1,:,:)**2)/n1/n2/n3)
!        call mpi_reduce(sigma_8_serial,sigma_8,1,mpi_real,mpi_sum,0,&
!            mpi_comm_world,ierr)
!        sigma_8=sigma_8/nlx/nly/nlz
!        write(*,*) 'sigma_8 = ',sigma_8
!
!        ! Output the final smoothed overdensity field
!        call RandomField_FT(deltag_w,deltagc_w,-1)
!        myid=1
!        call RandomField_Output(-99,deltag_w)
!        myid=0

        ! Get field in real space deltag_ as iFT of deltagc_
        call RandomField_FT(deltag_,deltagc_,-1)

!        open(unit=99,file='sigma_8.dat',access='stream')
!        write(99,pos=1) sigma_8
!        close(99)
!        sigma_8=sqrt( sum(deltag_(:,:,:)) / (n1*nlx*n2*nly*n3*nlz) )
!        write(*,*) 'sigma = ',sigma_8, n1,nlx
!
        ! Output the final overdensity field
        myid=1
        call RandomField_Output(-1,deltag_)
        myid=0

!        call test_random_number_stuff

    end subroutine make_zeta

    subroutine read_parameters
        ! Note the changes from the read_parameters subroutine used in a
        ! peak patch run:
        !
        !     1) The subroutine accepts an argument hpkvd_params, which is
        !        the name (including directory) of the text file printed by
        !        make_zeta.py after having read the parameter file.
        !
        !     2) Rather than passing the seed number at command line, we
        !        just read it in as the last value (see the line ending in
        !        !* below). ntasks, nmesh also passed here.
        open(unit=1,file=trim(hpkvd_params))
        read(1,*) ireadfield,ioutshear,&
          global_redshift,maximum_redshift,num_redshifts,&
          Omx,OmB,Omvac,h,nlx,nly,nlz,&
          dcore_box,dL_box,cenx,ceny,cenz,nbuff,next,ievol,&
          ivir_strat,fcoll_3,fcoll_2,fcoll_1,dcrit,iforce_strat,&
          TabInterpNx,TabInterpNy,TabInterpNz,&
          TabInterpX1,TabInterpX2,TabInterpY1,TabInterpY2,TabInterpZ1,TabInterpZ2,&
          wsmooth,rmax2rs,ioutfield, NonGauss, fNL, ilpt, iwant_field_part,&
          largerun,iseedFFT,ntasks,nmesh
        read(1,'(a)') fielddir
        read(1,'(a)') densfilein
        read(1,'(a)') filein
        read(1,'(a)') pkfile
        read(1,'(a)') filterfile
        read(1,'(a)') fileout
        read(1,'(a)') TabInterpFile
        close(1)
    end subroutine read_parameters

    subroutine smooth_field(field_in,field_out,k_min,R_smooth)
        complex(c_float_complex), pointer :: field_in(:,:,:),&
                                             field_out(:,:,:)
        integer ii,jj,kk,n1_nyquist,n2_nyquist,n3_nyquist
        real    k_min,R_smooth,k_x,k_y,k_z,kR,w

        ! Nyquist number for each axis
        n1_nyquist=n1/2+1
        n2_nyquist=n2/2+1
        n3_nyquist=n3/2+1

        ! OpenMP commands for parallel computing
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(n1_nyquist,n2_nyquist,n3_nyquist,k_min,R_smooth,field_in,field_out)
        do k=1,n3 ! iterate over z-axis
            kk=k
            ! fold the lattice in the z-axis at the Nyquist number
            ! kk={1,2,...,nyq1-1,nyq1,nyq1-1,...}
            if(k>n3_nyquist) kk=2*n3_nyquist-k
            k_z=(kk-1)*k_min ! z-component of wavevector

            do j=1,n2 ! iterate over y-axis
                jj=j 
                ! fold the lattice in the y-axis at the Nyquist number
                ! jj={1,2,...,nyq1-1,nyq1,nyq1-1,...}
                if(j>n2_nyquist) jj=2*n2_nyquist-j
                k_y=(jj-1)*k_min ! y-component of wavevector

                do i=1,n1_nyquist ! iterate over x-axis up to Nyquist number
                    ii=i
                    ! fold the lattice in the x-axis at the Nyquist number
!                    if(i>n1_nyquist) ii=2*n1_nyquist-i
                    k_x=(ii-1)*k_min ! x-component of wavevector

                    ! Wavenumber times smoothing scale
                    kR=sqrt(k_x**2+k_y**2+k_z**2)*R_smooth

                    ! Sets window function
                    if(wsmooth==0) then 
                        kR = kR/2
                        w = wgauss(kR) ! Gaussian window function
                    elseif(wsmooth==1) then 
                        w = wth(kR) ! real-space top-hat window function 
                    elseif(wsmooth==2) then
                        w = kR/R_smooth * wth(kR) ! sigma_1
                    elseif(wsmooth==3) then
                        w = (kR/R_smooth)**2 * wth(kR) ! sigma_2
                    elseif(wsmooth==4) then
                        w = k_x * wgauss(kR) ! to debug gradients
                    endif
                    if(kR==0) w=1.0 ! keeps DFT from blowing up at origin

                    ! Convolves input field with window function
                    if(wsmooth==4) then 
                        field_out(i,j,k) = cmplx(0.0,1.0)*field_in(i,j,k) * w
                    else 
                        field_out(i,j,k) = field_in(i,j,k) * w
                    endif
                enddo
            enddo
        enddo
        return
    end subroutine smooth_field

    ! Fourier-space Gaussian smoothing window function
    real function wgauss(x)
        implicit none
        real x
        wgauss = exp(-(x**2)/2.)
        return
    end function wgauss

    ! Fourier-space top-hat smoothing window function
    real function wth(x)
        implicit none
        real x
        wth = 3.*(sin(x)-x*cos(x))/x**3
        return
    end function wth

    subroutine test_random_number_stuff
        use random ! module src/modules/RandomField/random.f90
        implicit none
        integer q
        integer seed
        integer :: seeds(MaxRandNumStreams,IRandNumSize),&
                   seed1(MaxRandNumStreams,IRandNumSize),& 
                   seed2(MaxRandNumStreams,IRandNumSize) ! from random.f90
        real(dp) :: gauss_dev
        integer  :: poiss_dev
        real(kind=8) :: ranf_num,thing

        thing=50.

        call rans(ntasks,12345,seeds)
        call rans(ntasks,54321,seed1)
        call rans(ntasks,31442,seed2)

        open(unit=42,file=trim(run_dir)//'gaussdev.txt')
        open(unit=43,file=trim(run_dir)//'poissdev.txt')
        open(unit=44,file=trim(run_dir)//'ranf.txt')
        do q=1,10000
            call gaussdev(seeds(1,:), gauss_dev)
            write(42,*) gauss_dev
            call poissdev(seed1(1,:), thing, poiss_dev)
            write(43,*) poiss_dev
            call ranf(seed2(1,:), ranf_num)
            write(44,*) ranf_num
        enddo
        close(42)
        close(43)
        close(44)
    end subroutine test_random_number_stuff

end program make_zeta_runner
