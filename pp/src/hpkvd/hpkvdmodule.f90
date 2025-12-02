module hpkvdmodule

  ! INCLUDE MODULES
  use intreal_types      ! from src/modules/External, declares: dp, sp, i4b, i8b, dpc, spc, lgt
  use mpivars            ! from src/modules/External, declares: ierr, myid, ntasks, request, mpi_status, requests
  use timing_diagnostics ! from src/modules/External, subroutines: timer_begin, timer_end
  use input_parameters   ! from src/modules/GlobalVariables, declares input parameters
  use fftw_interface     ! from src/modules/RandomField

  use params             ! from src/modules/GlobalVariables
  use arrays             ! from src/hpkvd, subroutines: allocate_halos, allocate_boxes
  use RandomField        ! from src/modules/RandomField, generates unigrid initial conditions for simulation
  use SlabToCube         ! from src/modules/SlabToCube, routines for interfacing slab and cube domain decomps
  use TabInterp          ! from src/modules/TabInterp, tables for interpolation
  use memory_management  ! from src/modules/RandomField, subroutines: report_allocation,report_deallocation

  use io                 ! from src/hpkvd, subroutine: myfwrite

  use config_reader      ! from src/modules/ini_reader

  public :: hpkvd

contains 

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine hpkvd(parameters_filename, input_field, maketable, iseedFFT)

  ! ======================================================
  ! DECLARATIONS
  ! ======================================================

  ! PARAMETERS
  real,    parameter :: pi = 4.0 * atan(1.0)
  real,    parameter :: f3p = 4.0/3.0*pi
  integer :: maketable, iseedFFT

  ! LATTICE 
  integer  nhunt
  real     alattv(3),Sbar(3),Ebar(3,3),Sbar2(3),gradpk(3),gradpkrf(3),gradpkf(3)

  ! FILTERBANK 
  integer  idclus(nclmax)
  real     fcritv(nclmax),Rfclv(nclmax),sigrho_list(nclmax,3)
  integer  nic, nicz
  ! HALOS
  integer, allocatable :: Npk_eachtask(:), Npk_eachtaskl(:), Npk_begin(:)
  integer, allocatable :: ndump_eachtask(:), ndump_eachtaskl(:), ndump_begin(:)
  
  ! REDSHIFTS
  real, allocatable :: redshifts(:)  

  ! DATE AND TIME 
  character(8)          :: dnt_date
  character(10)         :: dnt_time
  character(5)          :: dnt_zone
  integer, dimension(8) :: dnt_values
  integer initial_time, final_time, inter_time, count_rate
  integer elapsed_hours, elapsed_minutes, elapsed_seconds
  double precision elapsed_time

  ! MEMORY USAGE
  real grid_mem, fftw_mem, s2c_mem, cat_mem, lat_mem, total_mem, n_fields, n_tot
  real local_grid_size, fftw_grid_size
  ! June 2024 Edit: need to define ngrid
  integer ngrid

  ! HALO STATISTICS
  real xpkmin,   xpkmax,   ypkmin,   ypkmax,   zpkmin,   zpkmax
  real xpkminl,  xpkmaxl,  ypkminl,  ypkmaxl,  zpkminl,  zpkmaxl
  real vxpkmin,  vxpkmax,  vypkmin,  vypkmax,  vzpkmin,  vzpkmax
  real vxpkminl, vxpkmaxl, vypkminl, vypkmaxl, vzpkminl, vzpkmaxl
  real vx2pkmin, vx2pkmax, vy2pkmin, vy2pkmax, vz2pkmin, vz2pkmax
  real vx2pkminl,vx2pkmaxl,vy2pkminl,vy2pkmaxl,vz2pkminl,vz2pkmaxl
  real RTHLmin
  real RTHLminl, RTHLmaxl
  real Fcollvmin, Fcollvmax, F_evmin, F_evmax, F_pvmin, F_pvmax
  real Fcollvminl,Fcollvmaxl,F_evminl,F_evmaxl,F_pvminl,F_pvmaxl
  real F_d2min, F_d2max, F_d2minl, F_d2maxl
  real sigrho, sigrho1, sigrho2
  real sigrhol, sigrho1l, sigrho2l
  ! RANDOM FIELD PARAMETERS
  character(len=512) rf_filepktab,rf_fileTktab,rf_outcode,rf_incode,rf_fielddir
  real               rf_boxsize,rf_ainit
  integer            rf_nmesh,rf_nbuff,rf_ntile,rf_seed,rf_ireadfield
  real               mi, fField
  integer*8          npart_eull, npart_eul, nlatt_total
  ! BOXES
  logical findpeaks, report_sigma
  double precision D, abar
  integer debug
  
  ! I/O
  character(len=512)      :: outputfilename, jboxout
  character(len=256)      :: parameters_filename
  character(len=8)        :: parameters_extension
  integer                 :: outnum
  integer(C_LONG)         :: pos_offset
  integer(C_LONG)         :: bytesout
  real(C_FLOAT),  pointer :: outbuffer(:)
  real(C_FLOAT),  target  :: RTHLmax, m1, global_redshiftl
  integer(C_INT), target  :: npk_tot
  logical                 :: file_exists
  integer(kind=mpi_offset_kind) offset_bytes  
  real(C_FLOAT), pointer, optional :: input_field(:,:,:) ! Random field if passed as a pointer
  ! EXTERNAL FUNCTIONS
  real evolve_ellipse_function
  external evolve_ellipse_function

  ! ======================================================
  ! EXECUTABLE STATMENTS
  ! ======================================================

  ! Set debug flag to 1 to have info printed
  debug = 0

  if (debug == 1) then
    ! Print the arguments
    write(*,*) 'maketable: ', maketable
    write(*,*) 'iseedFFT: ', iseedFFT
  end if
  ! HARD CODED PARAMETERS ARE IN MODULE params AND SET HERE
  verbose        = .true.
  rf_report      = 1
  report_memory  = .true.
  report_sigma   = .false.  !If true calls mpi all reduce to get exact sigma aross all processors
                            !Not needed unless exact calculation requires an exact value 
  ! HARD CODED TABLE PARAMETERS ARE IN MODULE TabInterp AND SET HERE
  TabInterpWrite = .True.  
  TabInterpOut   = -1

  ! REPORT DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     write(*,11) dnt_values(3),dnt_values(2),&
          dnt_values(1),dnt_values(5),dnt_values(6),&
          dnt_values(7),dnt_zone
     call system_clock(initial_time,count_rate)
  endif
  ! from Intel mpi module on Niagara
  if(report_memory) &
       call ReportMemory(C_CHAR_"Initial"//C_NULL_CHAR)

  ! READ PARAMETERS AND SET COSMOLOGY
! Get the file extension
  parameters_extension = adjustl(trim(get_file_extension(parameters_filename)))

  ! Check the file extension and call the appropriate subroutine
  if ((parameters_extension == "ini") .or. (parameters_extension == "conf")) then
      call read_parameters(parameters_filename) ! An ini-based approach to reading parameters, or conf-based (MUSIC)
  else if (parameters_extension == "bin") then
      call read_parameters_from_bin(parameters_filename) ! reads hpkvd_params.params (created by peak-patch.py)
  else
      print *, 'Error: Unsupported file extension: ', trim(parameters_extension)
      error stop 'Unsupported file extension'
  end if
  ! Print the parameters if needed
  !call verbose_writing
  if (debug == 1) then
    call print_legacy
  end if

  call Dlinear_cosmology(Omx,OmB,Omvac,h,iamcurved,dcurv) ! subroutine in src/cosmology/psubs_Dlinear.f90, makes tables for linear-order approximation to Zeldovich D(t)
  
  ! HOMOGENEOUS ELLIPSOID TABLES  
  call TabInterpInit() ! subroutine in src/modules/TabInterp.f90, allocates TabInterpArray
  TabInterpX1    = log10(TabInterpX1) ! Frhomin, TabInterpX1 read from parameter file with default 1.5
  TabInterpX2    = log10(TabInterpX2) ! Frhomax, TabInterpX2 read from parameter file with default 8.0

  if(maketable==1) then 
     timing_diagnostics_code='Making collapse table'
     call timer_begin ! from src/modules/External/timing_diagnostics.f90
     if(myid==0) then ! Checks that myid is processor 0 as assigned by MPI
        call TabInterpMake(evolve_ellipse_function) ! passes evolve_ellipse_function defined in src/hpkvd/peakvoidsubs.f90 to subroutine TabInterpMake from src/modules/TabInterp/TabInterp.f90
        call timer_end ! from src/modules/External/timing_diagnostics.f90, displays timing_diagnostics_code and cpu time for TabInterpMake() to run to completion
     endif
     call mpi_finalize(ierr) ! from Intel mpi module on Niagara
     stop
  else
     timing_diagnostics_code='Loading collapse table'
     call timer_begin     ! from src/modules/External/timing_diagnostics.f90
     call TabInterpLoad() ! from src/modules/TabInterp.f90
     if(myid==0) call timer_end ! from src/modules/External/timing_diagnostics.f90, displays timing_diagnostics_code and cpu time for TabInterpLoad() to run to completion
  endif

  ! INTIALIZE RANDOM FIELD
  rf_filepktab      = pkfile
  rf_fileTktab      = Tkfile
  rf_outcode        = filein
  rf_incode         = densfilein
  rf_fielddir       = fielddir
  rf_boxsize        = dcore_box*(nlx-1) + dL_box 
  rf_nmesh          = n1
  rf_nbuff          = nbuff
  rf_ntile          = nlx
  rf_ainit          = 1.0
  rf_seed           = iseedFFT
  rf_ireadfield     = ireadfield
call RandomField_init(rf_filepktab, rf_ireadfield, rf_outcode, rf_incode, rf_fielddir, rf_boxsize, &
                      rf_nmesh, rf_nbuff, rf_ntile, rf_ainit, NonGauss, fNL, A_nG, B_nG, R_nG, H_e, &
                      ng_seed, rf_fileTktab, h) ! subroutine in src/modules/RandomField/RandomField.f90

  ! INTIALIZE SLAB TO CUBE
  call SlabToCube_init(nlx,nly,nlz) ! subroutine in src/modules/SlabToCube/SlabToCube.f90

  ! MISCELLANEOUS PARAMETERS 

  n1_2     = n1/2   ! number of lattice cells in a tile (including buffers) divided by two
  nyq1     = n1/2+1 ! Nyquist number for tile
  n1xn2    = n1*n2  ! for cubic simulations, nmesh^2
  nF       = nFmax 
  alatt    = dL_box / n1   ! Lagrangian space lattice spacing
  boxlen   = alatt * n1    ! Lagrangian space tile side length (including buffers)
  alatt_2  = (1 / alatt)**2
  ak_min   = 2 * pi / boxlen
  alattv   = alatt
  vTHvir0  = 100 * h * sqrt(Omnr)
  vTHvir02 = vTHvir0*vTHvir0
  outnum   = 0
  ncore    = n1-2*nbuff
  n_tot    = ncore*nlx + 2*nbuff
  nlatt_total = (ncore*nlx)**3
  ! READ FILTERBANK AND MAKE BOXES
  call read_filterbank(nic,idclus,fcritv,Rfclv,Rfclmax) ! reads cluster ID, mass fcritv, filter scale Rfclv(ic)

  ! INITIALIZE LATTICE HUNT ARRAY
  !nhunt = 48.0/alatt
  nhunt = Rfclmax * 1.75 / alatt ! M.A.A. 30.11.16 
  nhunt = min(nbuff-1,nhunt)
  call icloud(npartmax,nhunt) ! subroutine in src/hpkvd/peakvoidsubs.f90, generates lattice positions
  ! defines npartmax=npart, the number of lattice sites within a radius nhunt (in lattice units) of
  ! a given cell of the lattice (npart is not declared so has implicit type integer as it starts with
  ! the letter n). The subroutine also defines ixsvec(1:npartmax,1:3) which defines 1 <= i <= npartmax
  ! 3-vectors ixsvec(i,:) defining the position in lattice units of the ith cell within radius nhunt,
  ! and defines irs2(1:npartmax), which is the radius squared of the corresponding point ixsvec.
  ! npartmax is used in hpkvd, and ixsvec, irs2 are used in a subsequent call to the subroutine
  ! get_homel in src/hpkvd/peakvoidsubs.f90

  ! REPORT MEMORY BUDGET
  fftw_grid_size  = total_local_sizes ! assigned in subroutine fftw_initialize() of src/modules/RandomField/fftw_interface.f90
  local_grid_size = (n1+2)*n1**2

  if(ilpt<2)  n_fields=4 ! Determines whether first-order (1LPT) or second-order (2LPT) Lagrangian perturbations to be used
  if(ilpt==2) n_fields=7 ! second-order perturbation theory is to be used in peak displacement

  ! June 2024 Edit: need to define ngrid. It is defined nowhere 
  ! in the code, so please change it if you need it.
  ngrid = 0

  fftw_mem = n_fields*4.*fftw_grid_size/1024**3              ! slab(4)
  s2c_mem  = 4. *(ngrid**3+local_nz*n1**2)/1024**3           ! sendbuff(4) + recbuff(4)
  grid_mem = (n_fields+1+1./4)*4.*local_grid_size/1024**3    ! rho(4), eta(12), mask(1), recvbuff(4)
  cat_mem  = 92.*(2*npkmax+npkmaxl)/1024**3                  ! catalogs(184)
  lat_mem  = 96.*npartmax/1024**3                            ! nshell(4),Sshell(12),SRshell(36),
                                                             ! Fbar(4),rad(4),indx(4),ixsvec(6)
                                                             ! irs2(4)

  total_mem   = fftw_mem + grid_mem + &     ! FFTW + GRID + CATALOG + LATTICE +
                cat_mem + lat_mem + &       ! S2C BUFFER
                s2c_mem

  if(myid==0) then
        write(*,*)
        write(*,*) '  FFTW memory:        ',fftw_mem,   ' GB'
        write(*,*) '  Grid memory:        ',grid_mem,   ' GB'
        write(*,*) '  S2C memory:         ',s2c_mem,    ' GB'
        !write(*,*) '     S2C, n1:         ',n1
        !write(*,*) '     S2C, ngrid:      ',ngrid
        !write(*,*) '     S2C, local_nz:   ',local_nz
        write(*,*) '  Catalog memory:     ',cat_mem,    ' GB'
        write(*,*) '  Lattice hunt memory:',lat_mem,    ' GB'
        write(*,*) '  Total memory:       ',total_mem  ,' GB'
        write(*,*) 
  endif

  ! ALLOCATE MASK AND FFT ARRAYS
  allocate(mask(n1,n2,n3))

  ! LOCAL, using subroutine from src/modules/RandomField/fftw_interface
  ! fftw_allocate_serial(f1,f2,c) assigns the target of a C pointer c to Fortran pointers f1 and f2 for serial operations (i.e. on one simulation tile of a parallel run)
  ! fftw_allocate_parallel(f1,f2,c) assigns the target of a C pointer c to Fortran pointers f1 and f2 for parallel operations (i.e. across tiles in a parallel run)
  call fftw_allocate_serial(   delta,   deltac,   deltap )
  call fftw_allocate_serial( delta_u, deltac_u, delta_up ) ! to save unsmoothed
  call fftw_allocate_serial(    etax,    etaxc,    etaxp )
  call fftw_allocate_serial(    etay,    etayc,    etayp )
  call fftw_allocate_serial(    etaz,    etazc,    etazp )

  if(ioutshear>=1) then
     call fftw_allocate_serial(     lapd,   lapdc, lapdp   )
     call fftw_allocate_serial(   lapd_u, lapdc_u, lapd_up )
     call fftw_allocate_parallel(  lapdg,  lapdgc, lapdgp  )
  endif
  ! GLOBAL, using subroutine from src/modules/RandomField/fftw_interface
  call fftw_allocate_parallel( deltag, deltagc, deltagp )
  call fftw_allocate_parallel(  etaxg,  etaxgc,  etaxgp )
  call fftw_allocate_parallel(  etayg,  etaygc,  etaygp )
  call fftw_allocate_parallel(  etazg,  etazgc,  etazgp )
  if (ilpt==2) then
     call fftw_allocate_serial( eta2x, eta2xc, eta2xp )
     call fftw_allocate_serial( eta2y, eta2yc, eta2yp )
     call fftw_allocate_serial( eta2z, eta2zc, eta2zp )
     !Use etaxg as F(q), see comment 1111 below 
     call fftw_allocate_parallel( eta2xg,  eta2xgc,  eta2xgp )
     call fftw_allocate_parallel( eta2yg,  eta2ygc,  eta2ygp )
     call fftw_allocate_parallel( eta2zg,  eta2zgc,  eta2zgp )
  endif

  if(report_memory) & ! calls subroutine in src/modules/External/memory_management.f90
       call ReportMemory(C_CHAR_"After grid allocations"//C_NULL_CHAR)

  ! ALLOCATE REDSHIFTS
  allocate(redshifts(num_redshifts))  

  ! MAKE BOXES AND BOX TO PROCESS MAPPING
  m = 0

  ! Set arrays xbox(i), ybox(i), zbox(i) as the coordinates of centres of each
  ! parallelisation tile (so i runs from 1 to ntile**3)
  call allocate_boxes ! subroutine in src/hpkvd/arrays_gen.f90
  call make_boxes     ! populates boxes allocated above
  do i=1,nboxes       ! nboxes = total simulation tiles (i.e. for square volume nboxes=ntile^3)

     ! xbx, ybx, zbx are set to the distances along each Cartesian axis from the origin to the near
     ! corner of the cubic simulation tile, such that rbx is the absolute distance from the origin
     ! to the near corner of the simulation tile
     xbx = abs(xbox(i)) - dcore_box / 2 ! dcore_box is the comoving sidelength
     ybx = abs(ybox(i)) - dcore_box / 2 ! of a simulation tile in Mpc so rbx
     zbx = abs(zbox(i)) - dcore_box / 2 ! is the comoving distance to the near

     ! Note: you may encounter an issue here where zbox(i) can get zeroed if the number of tiles is
     ! too large. This may be a result of MPI not being allocated enough RAM/node.

     rbx = sqrt(xbx**2+ybx**2+zbx**2)
     abx = afn(rbx) ! scale factor at rbx in light cone runs
     zbx = 1/abx-1  ! minimum redshift of tile i in light cone runs

     if(zbx > maximum_redshift.and.ievol==1) cycle ! evolve along line of sight in lightcone runs

     m = m + 1      ! For lightcone runs, lattice sites for which the redshift is greater than the max are dropped
     boxlist(m) = i
     
  enddo
 
  nboxes = m ! this is equivalent to no change nboxes = nboxes for non-lightcone runs or light cone runs where maximum_redshift is set to greater than the largest redshift of any halo in the box, which is typically the case
  neach  = ceiling((nboxes+0.0)/ntasks)

  ! GET GRID DENSITY ASSIGNMENT KERNEL TABLE
  call atab4 ! subroutine in src/hpkvd/peakvoidsubs.f90

  ! INITIALIZE REDSHIFTS
  dz = (maximum_redshift-global_redshift)/num_redshifts
  do i=1,num_redshifts
    redshifts(i) = global_redshift + (i-1)*dz
  enddo
  if(num_redshifts>1) ievol=0

  ! RANDOM FIELD REALIZATION
  if(ireadfield==0) then
     !Create Gaussian random field
     rf_seed = iseedFFT
     call RandomField_seed(rf_seed) ! subroutine in src/modules/RandomField/RandomField.f90, assigns the seed number used in subsequent RandomField_make calls
     call RandomField_make(-1, deltag, deltagc, deltagc) ! subroutine in src/modules/RandomField/RandomField.f90
     ! -1 is the ID for making an overdensity field, so deltag is populated with Gaussian random deviate, deltagc is set to FT(deltag), convolved with the powerspectrum and approprate transfer functions for the type of NonGaussianity used 

  elseif(ireadfield==1) then
     !Read in linear density field from external file
     call RandomField_Input(-1,deltag) ! subroutine in src/modules/RandomField/RandomField.f90
     if(ireadfield==1) call RandomField_FT(deltag,deltagc,1) ! subroutine in src/modules/RandomField/RandomField.f90, computes deltagc = FT(deltag)
  elseif(ireadfield > 1) then
     !Read in linear density field from a field, passed as a pointer
     call RandomField_Input_pointer(-1,deltag,input_field) ! subroutine in src/modules/RandomField/RandomField.f90
     if(ireadfield > 0) call RandomField_FT(deltag,deltagc,1) ! subroutine in src/modules/RandomField/RandomField.f90, computes deltagc = FT(deltag)
  endif

  ! 2LPT BEFORE 1LPT IN ORDER TO STORE INTERMEDIATE FIELD IN etaxg TO SAVE MEMORY
  ! 1111
  ! To calculate 2LPT need d^2(phi)/dx^i/dx^j \def phi,ij for 
  ! (i,j) = (1,1),(2,2),(3,3),(3,2),(3,1),(2,1)
  ! where phi,ij = -k_i*k_j/k^2 * delta(k)
  ! F(q) = (phi,33 * phi,22) + (phi,22 * phi,11) + (phi,33 * phi,11)
  !       -(phi,32)^2 - (phi,31)^2 - (phi,21)^2
  ! Finally, Psi^(2)_i = D_2 * i*k_i/k^2 * F(k)
  if (ilpt==2) then
     if(ireadfield<2) then
        call RandomField_make(-5,  eta2xg,  eta2xgc, deltagc)  ! phi,11
        call RandomField_make(-6,  eta2yg,  eta2ygc, deltagc)  ! phi,22
        call RandomField_make(-7,  eta2zg,  eta2zgc, deltagc)  ! phi,33
        do k=1,local_nz
           do j=1,n_tot
              do i=1,n_tot
                 etaxg(i,j,k) = eta2xg(i,j,k)*eta2yg(i,j,k) + &
                                eta2xg(i,j,k)*eta2zg(i,j,k) +eta2yg(i,j,k)*eta2zg(i,j,k) 
              enddo ! Determines order term in 2LPT F
           enddo    ! F = phi,11 * phi,22 + phi,11 * phi,33 + phi,22 * phi,33
        enddo       ! 
        ! the third field in these call statements, "deltagc" is the field to convolve with
        call RandomField_make(-8,  eta2xg,  eta2xgc, deltagc)  ! phi,21
        call RandomField_make(-9,  eta2yg,  eta2ygc, deltagc)  ! phi,31
        call RandomField_make(-10, eta2zg,  eta2zgc, deltagc)  ! phi,32
        do k=1,local_nz
           do j=1,n_tot
              do i=1,n_tot
                 etaxg(i,j,k) = etaxg(i,j,k) - eta2xg(i,j,k)**2 - &
                                eta2yg(i,j,k)**2 - eta2zg(i,j,k)**2 
              enddo ! Determines second order term in 2LPT F(q)
           enddo    ! F = F - phi,21^2 - phi,31^2 - phi,32^2
        enddo       ! 

        call RandomField_FT(etaxg,etaxgc,1)

        call RandomField_make(-11,  eta2xg,  eta2xgc, etaxgc) ! determines 2LPT x displacement
        if(ioutfield>1) call RandomField_Output(-11,eta2xg)
        call RandomField_make(-12,  eta2yg,  eta2ygc, etaxgc) ! determines 2LPT y displacement
        if(ioutfield>1) call RandomField_Output(-12,eta2yg)
        call RandomField_make(-13,  eta2zg,  eta2zgc, etaxgc) ! determines 2LPT z displacement
        if(ioutfield>1) call RandomField_Output(-13,eta2zg)
        
     elseif(ireadfield==2) then
        call RandomField_Input(-11,eta2xg)
        call RandomField_Input(-12,eta2yg)
        call RandomField_Input(-13,eta2zg)

     endif
  endif

  ! Next, whether or not 2LPT computation is done, the 1LPT fields are created
  if(ireadfield<2) then
     call RandomField_make(-2,  etaxg,  etaxgc, deltagc) ! determines 1LPT x displacement
     if(ioutfield>1) call RandomField_Output(-2,etaxg)
     call RandomField_make(-3,  etayg,  etaygc, deltagc) ! determines 1LPT y displacement
     if(ioutfield>1) call RandomField_Output(-3,etayg)  
     call RandomField_make(-4,  etazg,  etazgc, deltagc) ! determines 1LPT z displacement
     if(ioutfield>1) call RandomField_Output(-4,etazg)
     if(ioutshear>=1) call RandomField_make(-14,  lapdg,  lapdgc, deltagc) ! determines laplacian of overdensity delta

  elseif(ireadfield==2) then
     call RandomField_Input(-2,etaxg)
     call RandomField_Input(-3,etayg)
     call RandomField_Input(-4,etazg)

  endif
  
  if(ioutshear>=1) then
     if(ireadfield==2) call RandomField_Input(-14,lapdg)
     if(ioutfield>1)   call RandomField_Output(-14,lapdg)
  endif

  if(ireadfield<2) call RandomField_FT(deltag,deltagc,-1) ! compute deltag as iFT of deltagc
  if((ireadfield==0) .and. (ioutfield>0))  call RandomField_Output(-1,deltag) ! outputs deltag as the field Fvec_###Mpc...

  if(ioutfield==3) return ! Just creates 1LPT and 2LPT fields, but doesn't run peak patch

  ! INITIALIZE NUMBER OF PEAKS
  Npk        = 0
  ndump      = 0
  npart_eull = 0
  
  ! ALLOCATE HALOS
  call allocate_halos(num_redshifts,ioutshear) ! from src/hpkvd/arrays.f90, allocates arrays for halo search 

  ! IC Timing
  if(myid==0) then
     call system_clock(inter_time,count_rate)

     elapsed_seconds = int((inter_time - initial_time) / count_rate)
     write(*,501) elapsed_seconds
  endif

  ! ======================================================
  ! BEGIN MAIN PEAK PATCH CALCULATION

  ! LOOP OVER BOXES: myid is the CPU number running from 0 to ntasks-1 (inclusive), neach is the
  ! number of tasks run in series on each CPU, so ibox runs over all the tiles run as a separate
  ! task on each CPU
  do ibox=myid+1,ntasks*neach,ntasks
     findpeaks = .true.

     ! reduce in get_pks fails if every proccessor not involved
     if(ibox+ntasks > ntasks*neach) report_sigma = .false.
 
     ! Not all CPUs will necessarily be running the same number of tasks in series. To keep the
     ! parallelisation from failing, we introduce pseudo-periodic boundary conditions where the same
     ! tile is attributed to two different box IDs, but the second time peaks are not found.
     if(ibox>nboxes) then
       jbox = boxlist(ibox-nboxes)
       findpeaks = .false.
     else
       jbox = boxlist(ibox)
     endif

     ! Set (xbx,ybx,zbx) to the coordinates of the centre of the jbox^{th} cubic simulation tile in
     ! Mpc relative to the observer at the origin
     xbx = xbox(jbox)
     ybx = ybox(jbox)
     zbx = zbox(jbox)

     if(ievol==1) then ! for lightcone runs (redshift is a function of distance from centre of simulation box)
        rbx = sqrt(xbx**2+ybx**2+zbx**2) ! distance to the centre of the cubic simulation tile
        abx = afn(rbx)                   ! scale factor a(t), function of radial coordinate for lightcone runs
        Dbx = Dlinear(abx,rc,HD_a,D_a)   ! linear-order approximation of Zeldovich D(t)
     else ! for runs where the entire simulation volume is at constant redshift
        abx = 1./(1+global_redshift)   ! scale factor a(t), constant over entire box
        Dbx = Dlinear(abx,rc,HD_a,D_a) ! linear-order approximation of Zeldovich D(t)
     endif
     rdbx = 1./abx - 1 ! redshift of box (either a function of radial coord for lightcone runs or constant)

     delta = 0
     etax  = 0
     etay  = 0
     etaz  = 0

     ! SLAB TO CUBE EXCHANGE
     call timer_begin ! from src/modules/External/timing_diagnostics.f90
     timing_diagnostics_code='Slab to cube exchange'
     call get_cube_indices(xbx,ybx,zbx,nx,ny,nz) ! converts xbx,ybx,zbx from spatial coordinates (in Mpc) to lattice indices (integers)
     call SlabToCube_map(nx,ny,nz,nlx,nly,nlz) ! <SlabToCube_...> subroutines from:
     call SlabToCube_exchange(delta, deltag)   ! src/modules/SlabToCube/SlabToCube.f90
     call SlabToCube_exchange( etax,  etaxg)   ! These subroutines are used in parallelization as the simulation volume
     call SlabToCube_exchange( etay,  etayg)   ! is broken up into a series of smaller regions to optimize speed. Reads
     call SlabToCube_exchange( etaz,  etazg)   ! etajg from slab and writes etaj for tile

     if (ilpt==2) then
        eta2x  = 0
        eta2y  = 0 
        eta2z  = 0
        call SlabToCube_exchange( eta2x,  eta2xg) ! Reads in slabp=eta2xg, writes out tilep=eta2x
        call SlabToCube_exchange( eta2y,  eta2yg) ! Reads in slabp=eta2yg, writes out tilep=eta2y
        call SlabToCube_exchange( eta2z,  eta2zg) ! Reads in slabp=eta2zg, writes out tilep=eta2z
     endif

     call fftwf_execute_dft_r2c(plan_s, delta, deltac_u)  ! from FFTW module on Niagara, returns complex-valued deltac_u 
     if(ioutshear>=1) then                                ! as DFT of real-valued delta
        call SlabToCube_exchange( lapd,  lapdg) ! constructs aggregate lapd from lapdg on each parallelization tile
        call fftwf_execute_dft_r2c(plan_s, lapd, lapdc_u) ! from FFTW module on Niagara, returns complex-valued lapdc_u
     endif                                                ! as DFT of real-valued lapd

     ! PRE-SMOOTH TO DAMP HIGH-FREQUENCY NOISE
     ! call smooth_field(nyq1,ak_min,n1_2,0.5*alatt*2,fknyq,speq)
     ! deltac_u = deltac
     if(findpeaks) then ! IF FINDPEAKS {

     ! LOOP OVER REDSHIFTS {
     do ired=1,num_redshifts
     mask = 0
     if(ievol==1) then ! for lightcone runs
        cen = 0.5*(n1+1)                 ! centre of simulaiton volume
        rc  = sqrt(xbx**2+ybx**2+zbx**2) ! radial coordinate
        a   = afn(rc)                    ! scale factor a(t), function of radial coordinate for lightcone runs
        z   = 1./a-1                     ! redshift z also a function of radial coordinate
        D   = Dlinear(a,rc,HD_a,D_a)     ! linear-order approximation of Zeldovich D(t)
     elseif(ievol==0) then ! for nonlightcone runs
        z   = redshifts(ired)        ! redshift z
        a   = 1 / (1+z)              ! scale factor a(t), constant over entire box
        rc  = 0.0                    ! D(t) does not change with radial coordinate so we set rc to zero
        D   = Dlinear(a,rc,HD_a,D_a) ! linear-order approximation of Zeldovich D(t)
     endif
     nicz = nic ! nic is the number of entries in the filterbank, read in previously by subroutine read_filterbank()

     call timer_end ! from src/modules/External/timing_diagnostics.f90, displays timing_diagnostics_code and cpu time for subroutine calls to run to completion

     if(myid==0) write(*,101) ibox,nboxes,D,z,Npk ! writes formatted string to screen reporting box ID, D(t), z and number of peaks it contains

     timing_diagnostics_code='Smoothing'
     call timer_begin ! from src/modules/External/timing_diagnostics.f90

     ! LOOP OVER SMOOTHING SCALES
     if(myid==0) write(*,102) ! writes formatted string "Smoothing:" to screen
     Nprev=Npk
     nbad=0
     do ic=nicz,1,-1 ! start from largest smoothing scale Rfclv(nicz) and count down the smallest Rfclv(1)
        if(ioutshear>=1) then ! CALCULATE SIGMA_1 AND SIGMA_2 ON GRID BEFORE FINDING PEAKS
           ! ONLY CALCULATES TRUE MEAN FOR NPROC=NBOXES DUE TO REDUCE HERE 
           ! AND NOT AT VERY END. CLOSE ENOUGH TO TRUE MEAN FOR LARGE BOXES

           wsmooth = 2 !SIGMA_1
           call smooth_field(nyq1,ak_min,n1_2,Rfclv(ic),fknyq,speq) ! smooth on filter/q scale Rfclv(ic)
           call fftwf_execute_dft_c2r(iplan_s, deltac, delta) ! creates real-valued delta as DFT of complex-valued deltac
           delta = delta / real(n1)**3 ! normalisation

           if(report_sigma) then ! report_sigma=0 unless exact value needed
           sigrho1l=sqrt(sum(delta(1:n1,:,:)**2)/n1/n2/n3)
           call mpi_reduce(sigrho1l,sigrho1,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr) ! sigrho1 = mpi_sum( sigrho1l ), effectively copies sigrho1l on one of Peak Patch's parallel processors to sigrho1 on the root processor so it can be reported
           sigrho1 = sigrho1/nboxes ! normalisation
           endif

           wsmooth = 3 !SIGMA_2
           call smooth_field(nyq1,ak_min,n1_2,Rfclv(ic),fknyq,speq) ! smooth on filter scale Rfclv(ic)
           call fftwf_execute_dft_c2r(iplan_s, deltac, lapd) ! creates real-valued lapd as DFT of complex-valued deltac
           lapd = lapd / real(n1)**3 ! normalisation

           if(report_sigma) then ! report_sigma=0 unless exact value needed
           sigrho2l=sqrt(sum(lapd(1:n1,:,:)**2)/n1/n2/n3)
           call mpi_reduce(sigrho2l,sigrho2,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr) ! sigrho2 = mpi_sum( sigrho2l ), effectively copies sigrho1l on one of Peak Patch's parallel processors to sigrho1 on the root processor so it can be reported
           sigrho2 = sigrho2/nboxes ! normalisation
           sigrho_list(ic,2) = sigrho1
           sigrho_list(ic,3) = sigrho2
           endif

           wsmooth = 1 ! set to 1 so that peaks will be found after shear determined in this if statement
        endif

        ! FIND PEAKS
        call smooth_field(nyq1,ak_min,n1_2,Rfclv(ic),fknyq,speq) ! smooth on filter scale Rfclv(ic)

        ! FFT TO REAL SPACE
        call fftwf_execute_dft_c2r(iplan_s, deltac, delta) ! creates real-valued delta as DFT of complex-valued deltac
        delta = delta / real(n1)**3 ! normalisation

        ! CALCULATE VARIANCE ON SCALE Rfclv(ic)
        if(report_sigma) then ! report_sigma=0 unless exact value needed
           sigrhol=sqrt(sum(delta(1:n1,:,:)**2)/n1/n2/n3)
           call mpi_reduce(sigrhol,sigrho,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr) ! sigrho = mpi_sum( sigrhol )
           sigrho = sigrho/nboxes ! normalisation
           if(ioutshear>=1) sigrho_list(ic,1) = sigrho
        else
           sigrhol = sqrt(sum(delta(1:n1,:,:)**2)/n1/n2/n3) ! variance of delta
        endif

        ! FIND PEAKS
        Npk_in=Npk
        call get_pks(z,Npk,xbx,ybx,zbx,alattv,ired, Rfclv(ic)) ! subroutine in src/hpkvd/peakvoidsubs.f90, finds peaks in current simulation tile at filter scale Rfclv(ic) using nearest neighbours
        if (Npk.ge.Npkmaxl) then ! if Npk >= 200000, parameter Npkmaxl is assigned in src/hpkvd/arrays.f90
           write(*,*) 'Reached cluster limit ',Npkmaxl," myid=",myid
           call mpi_finalize(ierr)
           stop ! stops looking for new peaks if the max allowed value of 200000 peaks in the catalogue is exceded
        endif

        Npk_new=Npk-Npk_in ! number of peaks found by call to get_pks at current filter scale Rfclv(ic)
        if(myid==0) write(*,103) Rfclv(ic),sigrhol, Npk_new
        if(myid==0.and.ioutshear>=1) write(*,109) sigrho,sigrho1,sigrho2

        if(Npk_new.lt.1) cycle ! goes to the next value of ic if no peaks found at filter scale Rfclv(ic)

        do jpp=Npk_in+1,Npk  ! Preserves relationship between peak ID (jpp)
           iwas(jpp,ired)=ic ! redshift at which it virializes (ired) and
        enddo                ! element of the filter bank (ic)

     enddo! END LOOP OVER SMOOTHING SCALES

     if(myid==0) then ! Checks that myid is processor 0 as assigned by MPI
        if(ioutshear>=1) then !WRITE SIGMA_i list to file
           open(1,file="sigrho_out"//"."//str(iseedFFT),access="stream")
           write(1) nicz
           do ic=1,nicz
              write(1) Rfclv(ic),sigrho_list(ic,1),sigrho_list(ic,2),sigrho_list(ic,3)
           enddo
           close(1)
        endif
     endif

     if(myid==0) write(*,107) Npk-Nprev

     ! SET DELTA BACK TO UNSMOOTHED VALUE
     call fftwf_execute_dft_c2r(iplan_s, deltac_u, delta) ! creates real-valued delta as DFT of complex-valued deltac_u
     delta = delta / real(n1)**3 ! unsmoothed overdensity normalisation
     if(ioutshear>=1) then
        call fftwf_execute_dft_c2r(iplan_s, lapdc_u, lapd) ! creates real-valued lapd as DFT of complex-valued lapdc_u
        lapd = lapd / real(n1)**3 ! normalization, lapd is the laplacian of unsmoothed overdensity
     endif
     mask = 0
     if(myid==0) call timer_end ! from src/modules/External/timing_diagnostics.f90, displays timing_diagnostics_code and cpu time for subroutine calls to run to completion
     
     ! ELLIPSOIDAL COLLAPSE ON PEAKS
     timing_diagnostics_code='Ellipsoidal collapse'
     call timer_begin
     srms = 0.0
     grms = 0.0
     arms = 0.0
     ndump_i = 0
     do jpp=Nprev+1,Npk ! cycle over peaks found at current filter scale Rfclv(ic)
        if(mod((jpp-Nprev+1),5000)==0.and.myid==0) then ! Reports when ellipsoidal collapse completes every
           write(*,108) jpp-Nprev+1                     ! 5000 peaks to give an idea of how long collapse
        endif                                           ! process will take.

        ic=iwas(jpp,ired)     ! sets current cluster index ic to that corresponding for jpp
        jp=lagrange(jpp,ired) ! sets jp to the location of the cluster with respect to the lattice, the 3D array
        ! lagrange is assigned in subroutine get_pks of src/hpkvd/peakvoidsubs.f90 and uniquely maps elements
        ! (i,j,k) of the n1 x n2 x n3 lattice to an element of an array of length n1*n2*n3, i.e. if a peak
        ! (jpp,ired) is at lattice site (i,j,k), then lagrange(jpp,ired) = i+(j-1)*n1+(k-1)*n1*n2.

        if (ievol.eq.1) then ! adjusts peak redshifts based on radial position in light-cone runs
           radpk=sqrt(xpk(jpp,ired)**2+ypk(jpp,ired)**2+zpk(jpp,ired)**2)
           chipk=radpk
           apk=afn(chipk)
           redshiftpk=1/apk - 1
        else
           redshiftpk=z ! = 0.
        endif

        if(redshiftpk > maximum_redshift.and.ievol==1) then ! discards peaks that will not have collapsed
           ndump_i = ndump_i+1                              ! by z=maximum_redshift in light-cone runs
           cycle
        endif

        Fnupk=0. ! Note: Rfclv(ic) is the filter radius at which the current peak (index ic) was found,
        Fevpk=0. ! alatt is the sidelength of a lattice cell in Mpc/h, so Rfclv(ic)/alatt is the current
        Fpvpk=0. ! filter radius expressed in lattice units.
        Fd2pk=0. ! 
        ir2min=min( int((1.75*Rfclv(ic)/alatt)**2), int((40.0/alatt-1)**2)) !GS hardcode 45Mpc instead of 
        if(rmax2rs .gt. 0.0) then ! rmax2rs=0. in parameter file            !Rfclv(nic) so don't 
           npart = int(f3p*(rmax2rs*Rfclv(ic)/alatt)**3)                    !waste smoothing
           npart = min(npart,npartmax) ! npartmax is declared implicitly as an integer
        else                           ! because it starts with an n, and is assigned
           npart = npartmax            ! in a call to the subroutine "icloud" above
        endif ! npartmax is the number of lattice sites within a radius nhunt of the centre of the simulation
              ! tile where nhunt=min(nbuff-1, Rfclmax*1.75/alatt), or ~7/4 the max filter radius 

        ! HOMOGENEOUS ELLIPSOID CALCULATION
        call get_homel(npart,jp,alatt,ir2min,& ! subroutine in src/modules/hpkvd/peakvoidsubs.f90
             Sbar,RTHL,Srb,1+redshiftpk,Fnupk,Fevpk,&
             Fpvpk,Fd2pk,Ebar,Sbar2,Rfclv(ic),gradpkrf,gradpk,gradpkf,zformpk)
        ! This subroutine performs a radial integration of F profile to find Fbar=f_v crossing point
        ! closest to R^2 = ir2min * alatt. 

        ! get_homel returns the top-hat filter radius RTHL in lattice units, so here RTHLv is defined
        ! in terms of real-space units (Mpc/h) by multiplying RTHL by the real-space width of a lattice
        ! cell, alatt
        report_memory=.false.
        RTHLv(jpp,ired)=RTHL*alatt

        ! The integer ioutshear is assigned when a value is read from the parameter file
        ! (e.g. <run_directory>/param/param.params). If it is not set to zero, the components
        ! of the shear tensor for each peak will be appended to the main output. 
        if(ioutshear>0) zform(jpp,ired) = zformpk

        if(RTHLv(jpp,ired)<=0.) then ! If a peak doesn't fully collapse, or in a lightcone run where
           RTHLv(jpp,ired)=0.0       ! the peak lies beyond the horizon, get_homel sets RTHL=-1.0, so
           vTHvir(jpp,ired)=0.0      ! here we wipe those failed peaks

           vxpk(jpp,ired)=0.0
           vypk(jpp,ired)=0.0
           vzpk(jpp,ired)=0.0

           vx2pk(jpp,ired)  = 0.0
           vy2pk(jpp,ired)  = 0.0
           vz2pk(jpp,ired)  = 0.0

           Fcollv(jpp,ired) = 0.0

           if(ioutshear>0) F_ev(jpp,ired)   = 0.0 ! ellipticity
           if(ioutshear>0) F_pv(jpp,ired)   = 0.0 ! prolateness
           if(ioutshear>0) F_d2(jpp,ired)   = 0.0 ! laplacian of overdensity

           if(ioutshear>0) strain_bar(1,jpp,ired)=0.0 ! 
           if(ioutshear>0) strain_bar(2,jpp,ired)=0.0 ! Components of the strain tensor
           if(ioutshear>0) strain_bar(3,jpp,ired)=0.0 ! defined as in equation (2.16) of
           if(ioutshear>0) strain_bar(4,jpp,ired)=0.0 ! BM1.
           if(ioutshear>0) strain_bar(5,jpp,ired)=0.0 ! 
           if(ioutshear>0) strain_bar(6,jpp,ired)=0.0 ! 

           if(ioutshear>0) grad(1,jpp,ired)   = 0.0
           if(ioutshear>0) grad(2,jpp,ired)   = 0.0
           if(ioutshear>0) grad(3,jpp,ired)   = 0.0
           if(ioutshear>0) gradrf(1,jpp,ired) = gradpkrf(1)
           if(ioutshear>0) gradrf(2,jpp,ired) = gradpkrf(2)
           if(ioutshear>0) gradrf(3,jpp,ired) = gradpkrf(3)
           if(ioutshear>0) gradf(1,jpp,ired)  = 0.0
           if(ioutshear>0) gradf(2,jpp,ired)  = 0.0
           if(ioutshear>0) gradf(3,jpp,ired)  = 0.0

           ! iwas(jpp,ired)=-10000
           ndump_i = ndump_i+1
        else
           a  = 1./(1+redshiftpk) ! scale factor a(t)
           D  = Dlinear(a,rc,HD_a,D_a) ! linear-order Zeldovich time dependence D(t)
           
           !Changed from BM x=q-psi to standard x=q+psi LPT once and for all
           Sbar  = Sbar*D
           ! D_2 = -3/7*OmegaM**(-1/143)*D_1**2 Bouchet et al. 1995
           Sbar2 = -Sbar2 * 3./7*(Omnr/a**3 / (Omnr/a**3+Omvac))**(-1./143)*D**2 
           if(ilpt<2) Sbar2 = 0.0

           vxpk(jpp,ired)   = Sbar(1)
           vypk(jpp,ired)   = Sbar(2)
           vzpk(jpp,ired)   = Sbar(3)

           vx2pk(jpp,ired)  = Sbar2(1)
           vy2pk(jpp,ired)  = Sbar2(2)
           vz2pk(jpp,ired)  = Sbar2(3)

           Fcollv(jpp,ired) = Fnupk

           if(ioutshear>0) F_ev(jpp,ired)   = Fevpk ! ellipticity
           if(ioutshear>0) F_pv(jpp,ired)   = Fpvpk ! prolateness
           if(ioutshear>0) F_d2(jpp,ired)   = Fd2pk ! laplacian of overdensity
           
           if(ioutshear>0) strain_bar(1,jpp,ired)=Ebar(1,1) ! 
           if(ioutshear>0) strain_bar(2,jpp,ired)=Ebar(2,2) ! Components of the strain tensor
           if(ioutshear>0) strain_bar(3,jpp,ired)=Ebar(3,3) ! defined as in equation (2.16) of
           if(ioutshear>0) strain_bar(4,jpp,ired)=Ebar(2,3) ! BM1.
           if(ioutshear>0) strain_bar(5,jpp,ired)=Ebar(1,3) ! 
           if(ioutshear>0) strain_bar(6,jpp,ired)=Ebar(1,2) ! 

           if(ioutshear>0) grad(1,jpp,ired)   = gradpk(1) 
           if(ioutshear>0) grad(2,jpp,ired)   = gradpk(2)
           if(ioutshear>0) grad(3,jpp,ired)   = gradpk(3)
           if(ioutshear>0) gradrf(1,jpp,ired) = gradpkrf(1)
           if(ioutshear>0) gradrf(2,jpp,ired) = gradpkrf(2)
           if(ioutshear>0) gradrf(3,jpp,ired) = gradpkrf(3)
           if(ioutshear>0) gradf(1,jpp,ired)  = gradpkf(1)
           if(ioutshear>0) gradf(2,jpp,ired)  = gradpkf(2)
           if(ioutshear>0) gradf(3,jpp,ired)  = gradpkf(3)

           vE2 = vTHvir02 * Fnupk * Srb * RTHLv(jpp,ired)**2
           vTHvir(jpp,ired)=sqrt(abs(vE2))
           if(vE2.lt.0.0) vTHvir(jpp,ired)=-vTHvir(jpp,ired)              

        endif

        if (iwas(jpp,ired).gt.-10000) iwas(jpp,ired)=idclus(ic)
        if (idoboxf.eq.1) lagrange(jpp,ired)=jbox


     enddo

     ndump = ndump + ndump_i
     if(myid==0) call timer_end ! from src/modules/External/timing_diagnostics.f90, displays timing_diagnostics_code and cpu time for subroutine calls to run to completion

     if (largerun==1) then

        ! Write halos found by each tile to a seperate file
        ! This is needed for large runs with > 100 million halos
        ! in order to drastically speed up merge,
        ! where now only halos from tile neighbours are considered
        
        write(jboxout,'(I6.6)') jbox ! constructs 6 character entity name for output
        open(unit=3,file=trim(fileout)//'.'//trim(str(iseedFFT))//'_'//trim(jboxout),access='stream')
        
        write(3) Npk-Nprev-ndump_i, 36., global_redshift ! hardcoded max halo size of 36Mpc 
        do jpp=Nprev+1,Npk
           if(ioutshear==0) then
              if(RTHLv(jpp,ired) > 0.) then
                 write(3) xpk(jpp,ired),ypk(jpp,ired),zpk(jpp,ired),&
                      vxpk(jpp,ired),vypk(jpp,ired),vzpk(jpp,ired), RTHLv(jpp,ired),&
                      vx2pk(jpp,ired),vy2pk(jpp,ired),vz2pk(jpp,ired), Fcollv(jpp,ired)
              endif

           elseif(ioutshear>=1) then
              write(3) xpk(jpp,ired),ypk(jpp,ired),zpk(jpp,ired),&
                   vxpk(jpp,ired),vypk(jpp,ired),vzpk(jpp,ired), RTHLv(jpp,ired),&
                   vx2pk(jpp,ired),vy2pk(jpp,ired),vz2pk(jpp,ired), Fcollv(jpp,ired),&
                   F_ev(jpp,ired),F_pv(jpp,ired),strain_bar(1,jpp,ired),strain_bar(2,jpp,ired),&
                   strain_bar(3,jpp,ired),strain_bar(4,jpp,ired),strain_bar(5,jpp,ired),&
                   strain_bar(6,jpp,ired),f_d2(jpp,ired),zform(jpp,ired),grad(1,jpp,ired),&
                   grad(2,jpp,ired),grad(3,jpp,ired),gradf(1,jpp,ired),gradf(2,jpp,ired),&
                   gradf(3,jpp,ired),Rfclv(iwas(jpp,ired)),FcollvRf(jpp,ired),F_d2Rf(jpp,ired),&
                   gradrf(1,jpp,ired),gradrf(2,jpp,ired),gradrf(3,jpp,ired)

           endif
        enddo
        
        close(3)
     endif

     !Write out field particles depending on mask
     !Need to change to MPI_write for improved speed
     if(iwant_field_part == 1) then
        open(unit=4,file="field_particles.bin",access='stream')
        outnum       = 4
        offset_bytes = 3*int(4,8) + outnum*(ncore**3)*int(4,8)*(jbox-1) + 1

        cen1=0.5*(n1+1)
        cen2=0.5*(n2+1)
        cen3=0.5*(n3+1)

        do k=nbuff+1,n3-nbuff
           zc=zbx+alatt*(k-cen3)
           do j=nbuff+1,n2-nbuff
              yc=ybx+alatt*(j-cen2)
              do i=nbuff+1,n1-nbuff
                 xc=xbx+alatt*(i-cen1)
                 
                 npart_eull = npart_eull + (1-mask(i,j,k))
                 mi         = 1.-float(mask(i,j,k))        !if masked give mass of 0
                 
                 if(ievol==1) then
                    rc = sqrt(xc**2+yc**2+zc**2)
                    a   = afn(rc)
                    D   = Dlinear(a,rc,HD_a,D_a)
                 else
                    D = Dbx
                 endif
                 xpart_eul =  xc + etax(i,j,k)*D - eta2x(i,j,k) * 3./7*(Omnr/a**3 / (Omnr/a**3+Omvac))**(-1./143)*D**2  
                 ypart_eul =  yc + etay(i,j,k)*D - eta2y(i,j,k) * 3./7*(Omnr/a**3 / (Omnr/a**3+Omvac))**(-1./143)*D**2  
                 zpart_eul =  zc + etaz(i,j,k)*D - eta2z(i,j,k) * 3./7*(Omnr/a**3 / (Omnr/a**3+Omvac))**(-1./143)*D**2  
                 
                 write(4,pos=offset_bytes) mi, xpart_eul, ypart_eul, zpart_eul
                 offset_bytes = offset_bytes + outnum*4
              enddo
           enddo
        enddo
        close(4)
     endif
  enddo
  ! END LOOP OVER REDSHIFTS }

  endif ! ENDIF FINDPEAKS }

  enddo
  ! END LOOP OVER BOXES

  ! END MAIN PEAK PATCH CALCULATION
  ! ======================================================
  timing_diagnostics_code='Writing to disk'
  call timer_begin
  call mpi_barrier(mpi_comm_world,ierr)
  report_memory=.true.

  ! FREE ALL GRID MEMORY
  if(report_memory) &
       call ReportMemory(C_CHAR_"Before grid deallocations"//C_NULL_CHAR)

  ! LOCAL 
  call fftwf_free(deltap)
  call fftwf_free(etaxp)
  call fftwf_free(etayp)
  call fftwf_free(etazp)
  
  ! GLOBAL
  call fftwf_free(deltagp)
  call fftwf_free(etaxgp)
  call fftwf_free(etaygp)
  call fftwf_free(etazgp)
  if(ilpt==2) then
     call fftwf_free(eta2xgp)
     call fftwf_free(eta2ygp)
     call fftwf_free(eta2zgp)
  endif

  if(report_memory) &
       call ReportMemory(C_CHAR_"After grid deallocations"//C_NULL_CHAR)

  ! FIND OFFSET AND TOTAL NUMBER OF PEAKS
  allocate(Npk_eachtaskl(0:ntasks-1))
  allocate(Npk_eachtask (0:ntasks-1))
  allocate(Npk_begin    (0:ntasks-1))
  allocate(ndump_eachtaskl(0:ntasks-1))
  allocate(ndump_eachtask (0:ntasks-1))
  allocate(ndump_begin    (0:ntasks-1))


  Npk_begin           = 0
  Npk_eachtaskl       = 0
  if(ioutshear==0) then
     Npk_eachtaskl(myid) = Npk-ndump
  else
     Npk_eachtaskl(myid) = Npk
  endif
  npk_begin             = 0
  ndump_eachtaskl       = 0
  ndump_eachtaskl(myid) = ndump

  call mpi_allreduce(Npk_eachtaskl,Npk_eachtask,ntasks,mpi_integer,mpi_sum,&
                     mpi_comm_world,ierr)
  call mpi_allreduce(ndump_eachtaskl,ndump_eachtask,ntasks,mpi_integer,mpi_sum,&
                     mpi_comm_world,ierr)

  Npk_begin(0)   = 0
  ndump_begin(0) = 0
  do i=1,ntasks-1
    Npk_begin(i)=Npk_eachtask(i-1)+Npk_begin(i-1)
    ndump_begin(i)=ndump_eachtask(i-1)+ndump_begin(i-1)
  enddo

  if(ioutshear==0) then
     outnum = 11 ! position + 1lpt_vel + radius + 2lpt_vel + F
  elseif(ioutshear>=1) then
     outnum = 33 ! position + 1lpt_vel + radius + 2lpt_vel + F + e + p + strain + d2F + zform + grad(1:3) + gradf(1:3) + Rf + FRf + d2FRf  + gradRf(1:3) 
  endif

  Npk_tot      = Npk_begin(ntasks-1) + Npk_eachtask(ntasks-1)
  ndump_tot    = ndump_begin(ntasks-1) + ndump_eachtask(ntasks-1)

  ! LOAD OUTPUT BUFFERS
  fileout=trim(fileout)//'.'//str(iseedFFT)

  if(report_memory) &
       call ReportMemory(C_CHAR_"Before output buffer allocation"//C_NULL_CHAR)

  if(ioutshear==0) then
     outbufferp = fftwf_alloc_real(int(outnum*(Npk-ndump), C_SIZE_T))
     call c_f_pointer(outbufferp,outbuffer,[outnum*(Npk-ndump)])
  else
     outbufferp = fftwf_alloc_real(int(outnum*(Npk), C_SIZE_T))
     call c_f_pointer(outbufferp,outbuffer,[outnum*(Npk)])
  endif

  if(report_memory) &
       call ReportMemory(C_CHAR_"After output buffer allocation"//C_NULL_CHAR)

  m=0
  do ired=1,num_redshifts
     do i=1,Npk
        !G.S 24/09/2015 only write out good peaks (Npk-ndump)
        if(ioutshear==0) then
           if (RTHLv(i,ired) > 0.0) then  
              outbuffer(m+1) =   xpk(i,ired)
              outbuffer(m+2) =   ypk(i,ired)
              outbuffer(m+3) =   zpk(i,ired)
              outbuffer(m+4) =  vxpk(i,ired)
              outbuffer(m+5) =  vypk(i,ired)
              outbuffer(m+6) =  vzpk(i,ired)
              outbuffer(m+7) = RTHLv(i,ired)
              outbuffer(m+8)  =  vx2pk(i,ired)
              outbuffer(m+9)  =  vy2pk(i,ired)
              outbuffer(m+10) =  vz2pk(i,ired)
              outbuffer(m+11) =  Fcollv(i,ired)
              m = m + outnum
           endif
        else !Write out bad peaks for ioutshear>0
           outbuffer(m+1)  =   xpk(i,ired)
           outbuffer(m+2)  =   ypk(i,ired)
           outbuffer(m+3)  =   zpk(i,ired)
           outbuffer(m+4)  =  vxpk(i,ired)
           outbuffer(m+5)  =  vypk(i,ired)
           outbuffer(m+6)  =  vzpk(i,ired)
           outbuffer(m+7)  = RTHLv(i,ired)
           outbuffer(m+8)  =  vx2pk(i,ired)
           outbuffer(m+9)  =  vy2pk(i,ired)
           outbuffer(m+10) =  vz2pk(i,ired)
           outbuffer(m+11) = Fcollv(i,ired)
           if(ioutshear>=1) then 
              ! position + 1lpt_vel + radius + 2lpt_vel + F + e + p + strain + d2F + zform 
              ! + grad(1:3) + gradf(1:3) + Rf + FRf + d2FRf  + gradRf(1:3) 
              outbuffer(m+12) = F_ev(i,ired)
              outbuffer(m+13) = F_pv(i,ired)
              outbuffer(m+14) = strain_bar(1,i,ired)
              outbuffer(m+15) = strain_bar(2,i,ired)
              outbuffer(m+16) = strain_bar(3,i,ired)
              outbuffer(m+17) = strain_bar(4,i,ired)
              outbuffer(m+18) = strain_bar(5,i,ired)
              outbuffer(m+19) = strain_bar(6,i,ired)
              outbuffer(m+20) = f_d2(i,ired)
              outbuffer(m+21) = zform(i,ired)
              outbuffer(m+22) = grad(1,i,ired)
              outbuffer(m+23) = grad(2,i,ired)
              outbuffer(m+24) = grad(3,i,ired)
              outbuffer(m+25) = gradf(1,i,ired)
              outbuffer(m+26) = gradf(2,i,ired)
              outbuffer(m+27) = gradf(3,i,ired)
              outbuffer(m+28) = Rfclv(iwas(i,ired))
              outbuffer(m+29) = FcollvRf(i,ired)
              outbuffer(m+30) = F_d2Rf(i,ired)
              outbuffer(m+31) = gradrf(1,i,ired)
              outbuffer(m+32) = gradrf(2,i,ired)
              outbuffer(m+33) = gradrf(3,i,ired)

           endif
           m = m + outnum
        endif
     enddo
  enddo
  ! GATHER PEAK INFORMATION
  xpkminl   = minval(xpk)
  xpkmaxl   = maxval(xpk)
  ypkminl   = minval(ypk)
  ypkmaxl   = maxval(ypk)
  zpkminl   = minval(zpk)
  zpkmaxl   = maxval(zpk)

  vxpkminl  = minval(vxpk)
  vxpkmaxl  = maxval(vxpk)
  vypkminl  = minval(vypk)
  vypkmaxl  = maxval(vypk)
  vzpkminl  = minval(vzpk)
  vzpkmaxl  = maxval(vzpk)

  vx2pkminl  = minval(vx2pk)
  vx2pkmaxl  = maxval(vx2pk)
  vy2pkminl  = minval(vy2pk)
  vy2pkmaxl  = maxval(vy2pk)
  vz2pkminl  = minval(vz2pk)
  vz2pkmaxl  = maxval(vz2pk)

  Fcollvminl = minval(Fcollv)
  Fcollvmaxl = maxval(Fcollv)

  if(ioutshear>0) then
    F_evminl   = minval(F_ev)
    F_evmaxl   = maxval(F_ev)
    F_pvminl   = minval(F_pv)
    F_pvmaxl   = maxval(F_pv)
  endif

  RTHLminl = minval(RTHLv)
  RTHLmaxl = maxval(RTHLv)

  call mpi_reduce(xpkminl,xpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(xpkmaxl,xpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(ypkminl,ypkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(ypkmaxl,ypkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(zpkminl,zpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(zpkmaxl,zpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  call mpi_reduce(vxpkminl,vxpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vxpkmaxl,vxpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(vypkminl,vypkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vypkmaxl,vypkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(vzpkminl,vzpkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vzpkmaxl,vzpkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  call mpi_reduce(vx2pkminl,vx2pkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vx2pkmaxl,vx2pkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(vy2pkminl,vy2pkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vy2pkmaxl,vy2pkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  call mpi_reduce(vz2pkminl,vz2pkmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(vz2pkmaxl,vz2pkmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  call mpi_reduce(Fcollvminl,Fcollvmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)  
  call mpi_reduce(Fcollvmaxl,Fcollvmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)  
  call mpi_reduce(F_evminl,F_evmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)  
  call mpi_reduce(F_evmaxl,F_evmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)  
  call mpi_reduce(F_pvminl,F_pvmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)  
  call mpi_reduce(F_pvmaxl,F_pvmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)  
  call mpi_reduce(RTHLminl,RTHLmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(RTHLmaxl,RTHLmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  ! WRITE HEADER FOR FIELD PARTICLES
  call mpi_reduce(npart_eull, npart_eul, 1, mpi_integer8, mpi_sum, 0, mpi_comm_world, ierr)

  if(iwant_field_part == 1) then
     if(myid==0) then
        fField = real(npart_eul)/ real(nlatt_total) 
        open(unit=4,file="field_particles.bin",access='stream')
        offset_bytes = 1
        write(4,pos=offset_bytes) nlatt_total, fField
        close(4)
     endif
  endif
  ! WRITE HEADER IF myid == 0
  bytesout   = 4
  if(myid==0) then

    ! WRITE npk_tot
    pos_offset = 0
    call myfwrite(trim(fileout)//c_null_char,&
                  c_loc(npk_tot), bytesout, pos_offset) ! calls subroutine myfwrite in src/hpkvd/io.f90, which is a wrapper for the C function myfwrite in src/modules/External/myio.C

    ! WRITE RTHLmax
    pos_offset = pos_offset + bytesout
    call myfwrite(trim(fileout)//c_null_char,&
                  c_loc(RTHLmax), bytesout, pos_offset) ! calls subroutine myfwrite in src/hpkvd/io.f90, which is a wrapper for the C function myfwrite in src/modules/External/myio.C

    ! WRITE BOXREDSHIFT 
    pos_offset = pos_offset + bytesout
    if(ievol==1) then
       m1=-1.0
       call myfwrite(trim(fileout)//c_null_char,&
            c_loc(m1), bytesout, pos_offset) ! calls subroutine myfwrite in src/hpkvd/io.f90, which is a wrapper for the C function myfwrite in src/modules/External/myio.C
    elseif(ievol==0) then
       global_redshiftl=global_redshift
       call myfwrite(trim(fileout)//c_null_char,&
            c_loc(global_redshiftl), bytesout, pos_offset) ! calls subroutine myfwrite in src/hpkvd/io.f90, which is a wrapper for the C function myfwrite in src/modules/External/myio.C
    endif

  endif 

  ! WRITE PEAKS
  pos_offset = int(4,8)*(3+outnum*(Npk_begin(myid))) 

  if(ioutshear==0) then
     bytesout   = int(4,8)*outnum*(Npk-ndump)
  else
     bytesout   = int(4,8)*outnum*(Npk)
  endif
  if(Npk  >0) call myfwrite(trim(fileout)//c_null_char,&
       outbufferp, bytesout, pos_offset) ! calls subroutine myfwrite in src/hpkvd/io.f90, which is a wrapper for the C function myfwrite in src/modules/External/myio.C

  ! REPORT PEAK INFORMATION
  if(myid==0) then
     write(*,301) Npk_tot
     write(*,302) ndump_tot
     write(*,303) RTHLmin,RTHLmax
     write(*,304) xpkmin,xpkmax,ypkmin,ypkmax,zpkmin,zpkmax
     write(*,305) vxpkmin,vxpkmax,vypkmin,vypkmax,vzpkmin,vzpkmax
     write(*,306) vx2pkmin,vx2pkmax,vy2pkmin,vy2pkmax,vz2pkmin,vz2pkmax
     if(ioutshear>0) write(*,307) Fcollvmin,Fcollvmax,F_evmin,F_evmax,F_pvmin,F_pvmax
  endif

  call timer_end ! from src/modules/External/timing_diagnostics.f90, displays timing_diagnostics_code and cpu time for subroutine calls to run to completion

  ! DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     call system_clock(final_time,count_rate)

     elapsed_time    = (final_time - initial_time) / count_rate
     elapsed_hours   = int(elapsed_time / 3600)
     elapsed_minutes = int((elapsed_time - elapsed_hours * 3600)/60)
     elapsed_seconds = int(elapsed_time - elapsed_hours * 3600 &
                                        - elapsed_minutes * 60)
     
     write(*,401) dnt_values(3),dnt_values(2),dnt_values(1),&
       dnt_values(5),dnt_values(6),dnt_values(7),&
       dnt_zone,elapsed_hours,elapsed_minutes,elapsed_seconds

  endif

  return

 11 format(/,3x,61('-'),/,3x,'Peak-patch HPKVD running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)
101 format(/,3x,'Box ',i0,' of ',i0,/,3x,&
         'D = ',f5.3,/,3x,'z = ',f5.2,/,3x,'Npks = ',i0)
102 format(/,5x,'Smoothing:')
103 format(  5x,'R = ',f7.3,' Mpc: sigma = ',f5.3,2x,' Peaks found = ',i0)
105 format(  10x, ' max delta(R), R = ',f7.3,f7.3)
107 format(/,5x,'Homogeneous Ellipsoid running on ',i0,&
         ' peaks:')
108 format(  5x,'Peak',1x,i0,1x,'done')
109 format(  5x,'sigma_0 = ',f5.3,2x,'sigma_1 = ',f5.3,2x,'sigma_2 = ',f5.3)

201 format(3x,/, 'Gathering peak information on process 0')
301 format(/,3x, 'Total number of peaks kept:   ',11x,i0)
302 format(3x, 'Total number of peaks dumped: ',11x,i0)
303 format(/,3x,   'Minimum and maximum of peak radii:',&
           /,5x,1('(',1pe13.6',',1pe13.6,')',2x,/))
304 format(3x,   'Minimum and maximum of coordinates:',&
           /,5x,3('(',1pe13.6',',1pe13.6,')',/,5x))
305 format(3x,   'Minimum and maximum of 1LPT displacements:',&
           /,5x,3('(',1pe13.6',',1pe13.6,')',/,5x))
306 format(3x,   'Minimum and maximum of 2LPT displacements:',&
           /,5x,3('(',1pe13.6',',1pe13.6,')',/,5x))
307 format(3x,   'Minimum and maximum of F_v, e_v, p_v, and d2_v:',&
           /,5x,4('(',1pe13.6',',1pe13.6,')',/,5x))

401 format(/,3x,61('-'),/,3x,&
         'Peak-patch HPKVD finished running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         40x,'Total run time: ',i2.2,':',i2.2,':',i2.2,/&
         3x,61('-'),/)

501 format(/,3x,61('-'),/,3x,&
         'Initial Conditions Timing:',3x,i0,1x,'sec'/,&
         3x,61('-'),/)
                                    
end subroutine hpkvd


! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!SMOOTHES COMPLEX FIELD ON SCALE Rfclv AND RETURNS        

subroutine smooth_field(nyq1,ak_min,n1_2,Rfclv,fknyq,speq)

  use input_parameters
  use mpivars
  use arrays
  use fftw_interface

  sigrhol=0.0
  fknyq=(nyq1-1)*ak_min
  
  ! OpenMP commands for parallel computing
!$OMP PARALLEL DO &
!$OMP DEFAULT(FIRSTPRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(nyq1,ak_min,n1_2,Rfclv,delta,deltac,fknyq,speq)

  do k=1,n3 ! cycle over every lattice site along the 3 axis of a given parallelization tile
     kk=k
     if(k.ge.nyq1) kk=n3+2-k ! folds the lattice in the 3 axis at the Nyquist frequency kk={1,2,...,nyq1-1,nyq1,nyq1-1,...}
     fkz=(kk-1)*ak_min       ! coorinate in the Fourier conjugate 3 axis

     do j=1,n2 ! cycle over every lattice site along the 2 axis of a given parallelization tile
        jj=j
        if(j.ge.nyq1) jj=n2+2-j ! folds the lattice in the 2 axis at the Nyquist frequency kk={1,2,...,nyq1-1,nyq1,nyq1-1,...}
        fky=(jj-1)*ak_min       ! coordinate in the Fourier conjugate 2 axis

        do i=1,n1_2+1 ! cycle over first half of lattice sites along 1 axis of a given parallelization tile
           ii=i
           if(i>n1_2) ii=n1+2-i ! folds the lattice in the 1 axisy
           fkx=(ii-1)*ak_min    ! coordinate in the truncated Fourier conjugate 1 axis

           ! CONVOLVE WITH GAUSSIAN OR TOP-HAT FILTER                                                
           fkR=sqrt(fkx**2+fky**2+fkz**2)*Rfclv

           if(wsmooth==0) then
              fkR = fkR/2
              w = wgauss(fkR) ! k-space Gaussian window function
           elseif(wsmooth==1) then
              w = wth(fkR) ! k-space top-hat window function 
           elseif(wsmooth==2) then !to calculate sigma_1
              w = fkR/Rfclv * wth(fkR)
           elseif(wsmooth==3) then !to calculate sigma_2
              w = (fkR/Rfclv)**2 * wth(fkR)
           elseif(wsmooth==4) then !to debug gradients
              w = fkx * wgauss(fkR)
           endif

           if(fkR==0) w=1 ! keeps DFT from blowing up at origin

           ! Calls fftw to do the convolution integral
           if(wsmooth==4) then
              deltac(i,j,k) = cmplx(0.0,1.0)*deltac_u(i,j,k) * w
           else
              deltac(i,j,k) = deltac_u(i,j,k) * w
           endif
        enddo

     enddo
  enddo

  return
end subroutine smooth_field


! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine get_cube_indices(xbx,ybx,zbx,i,j,k)

  USE intreal_types
  use input_parameters
  use arrays
  use mpivars
  
  real xmin,ymin,zmin

  ! Minimum values of x,y,z coordinates for the entire simulation volume (in Mpc)
  xmin = cenx - nlx/2.*dcore_box ! Note: dcore_box read from parameter file, defined as boxsize/ntile, or
  ymin = ceny - nly/2.*dcore_box ! the sidelength of simulation cube in Mpc divided by the number of tiles
  zmin = cenz - nlz/2.*dcore_box ! used for parallel run. Also, the default setting is nlx=nly=nlz=ntile.

  ! (xbx,ybx,zbx) are the coordinates of the centre of a given cubic parallelisation tile (in Mpc),
  ! so starting with (i,j,k)=(1,1,1) for the tile with the smallest values for each of xbx, ybx and
  ! zbx, (i,j,k) are integers defining the tile relative to that tile.
  i = int((xbx - xmin)/dcore_box)+1
  j = int((ybx - ymin)/dcore_box)+1
  k = int((zbx - zmin)/dcore_box)+1

  return

end subroutine get_cube_indices

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

real function wgauss(x)

  implicit none
  real x
    
  wgauss = exp(-(x**2)/2.)

  return 

end function wgauss

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

real function wth(x)

  implicit none
  real x

  wth = 3.*(sin(x)-x*cos(x))/x**3

  return 

end function wth

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine read_external_field(xbx,ybx,zbx,&
     ibox,speq,speqx,speqy,speqz)

  USE intreal_types
  use input_parameters
  use arrays
  use mpivars
  
  character*512 filename,filebox
  
  integer ibox
  real    xbx,ybx,zbx
  
  complex speq(  :,:)
  
  complex speqx(  :,:),speqy(  :,:),speqz(  :,:)

  integer nxc,nyc,nzc,ntot
  real    boxsize, boxext
  
  ! FIND WHERE THIS BOX IS IN INPUT FIELD
  ntot    = nlx * ( n1 - 2 * nbuff )
  boxsize = dcore_box * nlx 
  boxext  = boxsize * next / ntot
  dx      = boxext / next 
  deltax  = xbx + boxext / 2
  deltay  = ybx + boxext / 2
  deltaz  = zbx + boxext / 2
  nxc     = int( deltax / dx + 1.1)
  nyc     = int( deltay / dx + 1.1)
  nzc     = int( deltaz / dx + 1.1)
  
  ! READ EXTERNAL FIELD 
  
  filename='Fvec_'//filein
  F = delta
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  filename='etax_'//filein
  F = etax
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  filename='etay_'//filein
  F = etay
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  filename='etaz_'//filein
  F = etaz
  call readsubbox(filename,nxc,nyc,nzc,ibox)
  
  ! eta = -displacement, so need to multiply by -1
  etax = -etax
  etay = -etay
  etaz = -etaz
  
  return
  
end subroutine read_external_field

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine readsubbox(filename,nxc,nyc,nzc,ibox)
  
  use mpivars
  use omp_lib

  use input_parameters
  use arrays
  use endian_utility

  character*512 filename

  integer nxc,nyc,nzc
  integer nx1,nx2,ny1,ny2,nz1,nz2
  integer ii,jj,kk
  real dum
  integer*8 m,ioffset,joffset,koffset,ne
  integer, dimension(3) :: sizes,sizese,start
  integer parallel_read
  ! Fix from April 2024
  integer offset

  ne=next
  ns=n1

  nx1 = nxc - int(ns/2) 
  nx2 = nxc + int(ns/2) - 1
  ny1 = nyc - int(ns/2) 
  ny2 = nyc + int(ns/2) - 1 
  nz1 = nzc - int(ns/2) 
  nz2 = nzc + int(ns/2) - 1 

  ! READ INPUT DATA
  parallel_read = 0
  if(parallel_read == 0) then
     
     offset = int(4,8)*(ibox-1)*n1**3+1
     open(unit=30,file=filename,access='stream')
     read(30,pos=int(offset)) F
     close(30)
     
  elseif(parallel_read == 1) then

     sizese=ne
     sizes=n1
     start(1)=nx1-1
     start(2)=ny1-1
     start(3)=nz1-1     

     call mpi_type_create_subarray(3,sizese,sizes,start,&
          mpi_order_fortran,mpi_real,newtype,ierr)
     call mpi_type_commit(newtype,ierr)
     call mpi_file_open(mpi_comm_world,filename,mpi_mode_rdonly,&
          mpi_info_null,ifile,ierr)  
     call mpi_file_set_view(ifile,0,mpi_real,newtype,"native",&
          mpi_info_null,ierr)
     call mpi_file_read_all(ifile,F,n1*n2*n3,mpi_real,&
          mpi_status_ignore,ierr)
     call mpi_file_close(ifile,ierr)

  else

     !offset = int(4,8)*(ibox-1)*n1**3
     offset = INT((ibox-1)*n1**3, 8)
     call mpi_file_open(mpi_comm_world,filename,mpi_mode_rdonly,&
          mpi_info_null,ifile,ierr)  
     call mpi_file_set_view(ifile,offset,mpi_real,newtype,"native",&
          mpi_info_null,ierr)
     call mpi_file_read(ifile,F,n1*n2*n3,mpi_real,&
          mpi_status_ignore,ierr)
     call mpi_file_close(ifile,ierr)

  endif

  if(big_endian()) then ! ASSUMES INPUT DATA IS LITTLE ENDIAN
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(F) 
    do k=1,n3
      do j=1,n2
        do i=1,n1
          F(i,j,k)=swap_endian(F(i,j,k))
        enddo
      enddo
    enddo
  endif

  return

end subroutine readsubbox
 
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine make_boxes

 use mpivars
 use arrays
 USE intreal_types
 use input_parameters
 use mpivars

 nboxes=nlx*nly*nlz

 i=0
 cen1=0.5*float(nlx+1)
 cen2=0.5*float(nly+1)
 cen3=0.5*float(nlz+1)
 do k3=1,nlz
    zBh=cenz+(k3-cen3)*dcore_box
    do k2=1,nly
       yBh=ceny+(k2-cen2)*dcore_box
       do k1=1,nlx
          xBh=cenx+(k1-cen1)*dcore_box

          i=i+1
          xbox(i)=xBh
          ybox(i)=yBh
          zbox(i)=zBh

       enddo
    enddo
 enddo

 return

end subroutine make_boxes

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine read_filterbank(nic,idclus,fcritv,Rfclv,Rfclmax)
  USE intreal_types
  use input_parameters

  dimension fcritv(*),Rfclv(*),idclus(*)

  open(1,file=filterfile,status='old',form='formatted')

  read(1,*) nic

  Rfclmax=0.0
  do ic=1,nic
     read(1,*) idclus(ic),fcritv(ic),Rfclv(ic)
     Rfclmax=max(Rfclmax,Rfclv(ic)) ! typically Rfclmax=Rfclv(nic)=34.0 [Mpc/h]
  enddo
  close(1)

  return

end subroutine read_filterbank

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

subroutine read_parameters_from_bin(parameters_filename)

  use input_parameters ! from src/modules/GlobalVariables, declares input parameters
  use arrays           ! from src/hpkvd, subroutines: allocate_halos, allocate_boxes
  use params           ! from src/modules/GlobalVariables
  use TabInterp        ! from src/modules/TabInterp/TabInterp.f90
  integer str_len,offset
  character(len=256)      :: parameters_filename
  ! READ INPUT PARAMETERS

  open(unit=1,file=trim(parameters_filename),access='stream')
  read(1,pos=1) ireadfield       , ioutshear        , global_redshift , maximum_redshift , &
                num_redshifts    , Omx              , OmB             , Omvac            , &
                h                , nlx              , nly             , nlz              , &
                dcore_box        , dL_box           , cenx            , ceny             , &
                cenz             , nbuff            , next            , ievol            , &
                ivir_strat       , fcoll_3          , fcoll_2         , fcoll_1          , &
                dcrit            , iforce_strat     , TabInterpNx     , TabInterpNy      , &
                TabInterpNz      , TabInterpX1      , TabInterpX2     , TabInterpY1      , &
                TabInterpY2      , TabInterpZ1      , TabInterpZ2     , wsmooth          , &
                rmax2rs          , ioutfield        , NonGauss        , fNL              , &
                A_nG             , B_nG             , R_nG            , H_e              , &
                ng_seed          , ilpt             , iwant_field_part, largerun

  offset=4*48
  str_len=0
  call read_string_from_bin(offset,str_len,fielddir     )
  call read_string_from_bin(offset,str_len,densfilein   )
  call read_string_from_bin(offset,str_len,filein       )
  call read_string_from_bin(offset,str_len,pkfile       )
  call read_string_from_bin(offset,str_len,Tkfile       )
  call read_string_from_bin(offset,str_len,filterfile   )
  call read_string_from_bin(offset,str_len,fileout      )
  call read_string_from_bin(offset,str_len,TabInterpFile)
  close(1)
  nstepmax      = 10000
  iwant_evmap   = 3     ! 1: turn vir z vs e_v, 2: z vs p_v for fixed e_v, 3: table
  iwant_rd      = 1
  zinit_fac_ell = 20
  tfac          = 0.01  ! fraction of local 1-axis Hubble time for dt
  ihard         = 1
  ! WRITE INPUT PARAMETERS IF VERBOSE=TRUE
end subroutine read_parameters_from_bin

subroutine verbose_writing()

  if(verbose.and.myid==0) then
    call print_parameter_char('Boolean Switches:', '')
    call print_parameter_int('  hpkvd_params:', hpkvd_params)
    call print_parameter_int('  compile_hpkvd:', compile_hpkvd)
    call print_parameter_int('  create_filterbank:', create_filterbank)
    call print_parameter_int('  merge_params:', merge_params)
    call print_parameter_int('  compile_merge:', compile_merge)
    call print_parameter_int('  map_params:', map_params)
    call print_parameter_int('  compile_maps:', compile_maps)
    call print_parameter_int('  batch:', batch)
    call print_parameter_int('  submit:', submit)

    call print_parameter_char('Host Computing Cluster:', '')
    call print_parameter_char('  machine:', machine)
    call print_parameter_char('  submit_command:', submit_command)

    call print_parameter_char('Peak Patch Main:', '')
    call print_parameter_int('  seed:', seed)
    call print_parameter_char('  run_name:', run_name)
    call print_parameter_char('  short_name:', short_name)
    call print_parameter_char('  runtype:', runtype)

    call print_parameter_char('Pixellation and Parallelization:', '')
    call print_parameter_real('  boxsize:', boxsize)
    call print_parameter_int('  nmesh:', nmesh)
    call print_parameter_int('  nbuff:', nbuff)
    call print_parameter_int('  ntile:', ntile)
    call print_parameter_int('  largerun:', largerun)

    call print_parameter_char('Cluster Parameters for Peak Patch:', '')
    call print_parameter_char('  tlimit:', tlimit)
    call print_parameter_int('  nnodes:', nnodes)
    call print_parameter_int('  tpnode:', tpnode)
    call print_parameter_int('  ntasks:', ntasks)
    call print_parameter_int('  ncpus:', ncpus)
    call print_parameter_int('  nompth:', nompth)

    call print_parameter_char('Cluster Parameters for Websky:', '')
    call print_parameter_char('  tlimit_map:', tlimit_map)
    call print_parameter_int('  nnodes_map:', nnodes_map)
    call print_parameter_int('  ppn_map:', ppn_map)
    call print_parameter_int('  ntasks_map:', ntasks_map)
    call print_parameter_int('  np_map:', np_map)
    call print_parameter_int('  nompth_map:', nompth_map)

    call print_parameter_char('Type of Peak Patch:', '')
    call print_parameter_int('  ievol:', ievol)
    call print_parameter_int('  num_redshifts:', num_redshifts)
    call print_parameter_real('  maximum_redshift:', maximum_redshift)
    call print_parameter_real('  global_redshift:', global_redshift)

    call print_parameter_char('Peak Displacement:', '')
    call print_parameter_int('  ilpt:', ilpt)
    call print_parameter_int('  ioutfield:', ioutfield)
    call print_parameter_int('  ireadfield:', ireadfield)
    call print_parameter_int('  iwant_field_part:', iwant_field_part)
    call print_parameter_char('  fielddir:', fielddir)
    call print_parameter_char('  densfilein:', densfilein)
    call print_parameter_char('  densfileout:', densfileout)

    call print_parameter_char('Non-Gaussianities:', '')
    call print_parameter_int('  NonGauss:', NonGauss)
    call print_parameter_real('  fNL:', fNL)
    call print_parameter_real('  A_nG:', A_nG)
    call print_parameter_real('  B_nG:', B_nG)
    call print_parameter_real('  R_nG:', R_nG)
    call print_parameter_real('  m_phi:', m_phi)
    call print_parameter_real('  m_chi:', m_chi)
    call print_parameter_real('  phi_w:', phi_w)
    call print_parameter_real('  phi_p:', phi_p)
    call print_parameter_real('  vev:', vev)
    call print_parameter_real('  m_tach:', m_tach)
    call print_parameter_real('  a_e:', a_e)
    call print_parameter_int('  ng_seed:', ng_seed)

    call print_parameter_char('Merging Algorithm:', '')
    call print_parameter_int('  iZeld:', iZeld)
    call print_parameter_int('  ntilemerge:', ntilemerge)
    call print_parameter_int('  ntasksmerge:', ntasksmerge)
    call print_parameter_int('  iwrap:', iwrap)

    call print_parameter_char('Websky Parameters:', '')
    call print_parameter_char('  maps:', maps)
    call print_parameter_int('  nside_map:', nside_map)
    call print_parameter_int('  npix_map:', npix_map)
    call print_parameter_real('  fov_map:', fov_map)
    call print_parameter_real('  zmin_map:', zmin_map)
    call print_parameter_real('  zmax_map:', zmax_map)
    call print_parameter_char('  tabfile_map:', tabfile_map)
    call print_parameter_char('  tabfile_sfr:', tabfile_sfr)
    call print_parameter_int('  model_map:', model_map)
    call print_parameter_int('  scramble_map:', scramble_map)
    call print_parameter_int('  center_map:', center_map)
    call print_parameter_int('  chihview_map:', chihview_map)
    call print_parameter_int('  PSZcut_map:', PSZcut_map)
    call print_parameter_int('  ellmax:', ellmax)

    call print_parameter_char('Cosmology Parameters:', '')
    call print_parameter_real('  Omx:', Omx)
    call print_parameter_real('  OmB:', OmB)
    call print_parameter_real('  Omvac:', Omvac)
    call print_parameter_real('  h:', h)
    call print_parameter_real('  ns:', ns)
    call print_parameter_real('  As:', As)
    call print_parameter_real('  sigma8:', sigma8)
    call print_parameter_real('  tau:', tau)
    call print_parameter_real('  mnu:', mnu)
    call print_parameter_char('  pkfile:', pkfile)

    call print_parameter_char('Ellipsoidal Collapse:', '')
    call print_parameter_int('  ioutshear:', ioutshear)
    call print_parameter_real('  rmax2rs:', rmax2rs)
    call print_parameter_int('  wsmooth:', wsmooth)
    call print_parameter_real('  Rsmooth_max:', Rsmooth_max)
    call print_parameter_char('  rapi:', rapi)
    call print_parameter_char('  TabInterpFile:', TabInterpFile)
    call print_parameter_int('  TabInterpNx:', TabInterpNx)
    call print_parameter_int('  TabInterpNy:', TabInterpNy)
    call print_parameter_int('  TabInterpNz:', TabInterpNz)
    call print_parameter_real('  TabInterpX1:', TabInterpX1)
    call print_parameter_real('  TabInterpX2:', TabInterpX2)
    call print_parameter_real('  TabInterpY1:', TabInterpY1)
    call print_parameter_real('  TabInterpY2:', TabInterpY2)
    call print_parameter_real('  TabInterpZ1:', TabInterpZ1)
    call print_parameter_real('  TabInterpZ2:', TabInterpZ2)
    call print_parameter_char('  filterfile:', filterfile)
    call print_parameter_int('  iforce_strat:', iforce_strat)
    call print_parameter_int('  ivir_strat:', ivir_strat)
    call print_parameter_real('  dcrit:', dcrit)
    call print_parameter_real('  fcoll_1:', fcoll_1)
    call print_parameter_real('  fcoll_2:', fcoll_2)
    call print_parameter_real('  fcoll_3:', fcoll_3)

    call print_parameter_char('Lattice Parameters for hpkvd:', '')
    call print_parameter_int('  nsub:', nsub)
    call print_parameter_int('  next:', next)
    call print_parameter_real('  dcore_box:', dcore_box)
    call print_parameter_real('  cellsize:', cellsize)
    call print_parameter_real('  buffersize:', buffersize)
    call print_parameter_real('  dL_box:', dL_box)
    call print_parameter_real('  mlatt:', mlatt)
    call print_parameter_real('  cenx:', cenx)
    call print_parameter_real('  ceny:', ceny)
    call print_parameter_real('  cenz:', cenz)
    call print_parameter_int('  nlx:', nlx)
    call print_parameter_int('  nly:', nly)
    call print_parameter_int('  nlz:', nlz)
    call print_parameter_int('  n1:', n1)
    call print_parameter_int('  n2:', n2)
    call print_parameter_int('  n3:', n3)

    call print_parameter_char('Advanced Merge Parameters:', '')
    call print_parameter_int('  iLexc:', iLexc)
    call print_parameter_int('  iLmrg:', iLmrg)
    call print_parameter_int('  iFexc:', iFexc)
    call print_parameter_int('  iFmrg:', iFmrg)
  endif
end subroutine verbose_writing

subroutine print_legacy()
  if(verbose.and.myid==0) then
     write(*,*) 'ireadfield  ', ireadfield
     write(*,*) 'ioutshear  ', ioutshear
     write(*,*) 'global_redshift  ', global_redshift
     write(*,*) 'maximum_redshift  ', maximum_redshift
     write(*,*) 'num_redshifts  ', num_redshifts
     write(*,*) 'Omx  ', Omx
     write(*,*) 'OmB  ', OmB
     write(*,*) 'Omvac  ', Omvac
     write(*,*) 'h  ', h
     write(*,*) 'nlx  ', nlx
     write(*,*) 'nly  ', nly
     write(*,*) 'nlz  ', nlz
     write(*,*) 'dcore_box  ', dcore_box
     write(*,*) 'dL_box  ', dL_box
     write(*,*) 'cenx  ', cenx
     write(*,*) 'ceny  ', ceny
     write(*,*) 'cenz  ', cenz
     write(*,*) 'nbuff  ', nbuff
     write(*,*) 'next  ', next
     write(*,*) 'ievol  ', ievol
     write(*,*) 'ivir_strat  ', ivir_strat
     write(*,*) 'fcoll_3  ', fcoll_3
     write(*,*) 'fcoll_2  ', fcoll_2
     write(*,*) 'fcoll_1  ', fcoll_1
     write(*,*) 'dcrit  ', dcrit
     write(*,*) 'iforce_strat  ', iforce_strat
     write(*,*) 'TabInterpNx  ', TabInterpNx
     write(*,*) 'TabInterpNy  ', TabInterpNy
     write(*,*) 'TabInterpNz  ', TabInterpNz
     write(*,*) 'TabInterpX1  ', TabInterpX1
     write(*,*) 'TabInterpX2  ', TabInterpX2
     write(*,*) 'TabInterpY1  ', TabInterpY1
     write(*,*) 'TabInterpY2  ', TabInterpY2
     write(*,*) 'TabInterpZ1  ', TabInterpZ1
     write(*,*) 'TabInterpZ2  ', TabInterpZ2
     write(*,*) 'wsmooth  ', wsmooth
     write(*,*) 'rmax2rs  ', rmax2rs
     write(*,*) 'ioutfield  ', ioutfield
     write(*,*) 'NonGauss  ', NonGauss
     write(*,*) 'fNL  ', fNL
     write(*,*) 'A_nG  ', A_nG
     write(*,*) 'B_nG  ', B_nG
     write(*,*) 'R_nG  ', R_nG
     write(*,*) 'H_e  ', H_e
     write(*,*) 'ng_seed  ', ng_seed
     write(*,*) 'ilpt  ', ilpt
     write(*,*) 'iwant_field_part  ', iwant_field_part
     write(*,*) 'largerun  ', largerun
     write(*,*) 'fielddir  ', trim(fielddir)
     write(*,*) 'densfilein  ', trim(densfilein)
     write(*,*) 'filein  ', trim(filein)
     write(*,*) 'pkfile  ', trim(pkfile)
     write(*,*) 'Tkfile  ', trim(Tkfile)
     write(*,*) 'filterfile  ', trim(filterfile)
     write(*,*) 'fileout  ', trim(fileout)
     write(*,*) 'TabInterpFile  ', trim(TabInterpFile)
  endif
end subroutine print_legacy


subroutine print_parameter_int(name, value)
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  character(len=32) :: str

  write(str, '(I0)') value
  call print_formatted(name, str)
end subroutine print_parameter_int

subroutine print_parameter_real(name, value)
  character(len=*), intent(in) :: name
  real, intent(in) :: value
  character(len=32) :: str

  write(str, '(F0.7)') value
  call print_formatted(name, str)
end subroutine print_parameter_real

subroutine print_parameter_char(name, value)
  character(len=*), intent(in) :: name, value

  call print_formatted(name, value)
end subroutine print_parameter_char

subroutine print_formatted(name, value)
  character(len=*), intent(in) :: name, value
  integer :: name_length, padding_length, total_length

  ! Define the total length of the formatted line
  total_length = 35

  ! Calculate the length of the name and the padding
  name_length = len_trim(name)
  padding_length = total_length - name_length

  ! Print the formatted line
  write(*, '(A, A)', advance='no') trim(name), repeat(' ', padding_length)
  print *, trim(value)
end subroutine print_formatted

subroutine read_string_from_bin(offset,str_len,str)
  integer offset, str_len, i
  character(len=512) str
  character(len=1) letter
  offset=offset+str_len
  read(1,pos=1+offset) str_len
  offset=offset+4
  str=' '
  do i=1,str_len
      read(1,pos=offset+i) letter
      str=trim(str)//letter
  enddo
end subroutine read_string_from_bin

function get_file_extension(filepath) result(extension)
    implicit none
    character(len=*), intent(in) :: filepath
    character(len=32) :: extension
    integer :: i

    extension = ""

    ! Find the position of the last dot in the filepath
    do i = len(filepath), 1, -1
        if (filepath(i:i) == '.') then
            extension = filepath(i+1:)
            return
        end if
    end do

end function get_file_extension

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

end module hpkvdmodule
