module config_reader

  use input_parameters
  use cosmoparams
  use params
  use mpivars
  use TabInterp
  use arrays

  implicit none

  ! Declare variables for boolean switches
  integer :: hpkvd_params, compile_hpkvd, create_filterbank, merge_params
  integer :: compile_merge, map_params, compile_maps, batch, submit
  
  ! Declare variables for host computing cluster
  character(len=32) :: machine, submit_command
  character(len=128) :: run_dir
  
  ! Declare variables for peak patch main parameters
  integer :: seed
  character(len=32) :: run_name, short_name, runtype 
  character(len=128) :: default_fileout
  character(len=128) :: default_filein
  
  ! Declare variables for pixellation and parallelization
  real :: boxsize
  integer :: nmesh, ntile
  
  ! Declare variables for cluster parameters for peak patch
  character(len=8) :: tlimit
  integer :: nnodes, tpnode, ncpus, nompth
  
  ! Declare variables for cluster parameters for websky
  character(len=8) :: tlimit_map
  integer :: nnodes_map, ppn_map, ntasks_map, np_map, nompth_map
  
  ! Declare variables for type of peak patch
  ! integer :: ng_seed

  ! Declare variables for peak displacement
  character(len=64) :: densfileout
  
  ! Declare variables for non-Gaussianities
  real :: m_phi, m_chi, phi_w, phi_p, vev, m_tach, a_e
  
  ! Declare variables for merging algorithm
  integer :: iZeld, ntilemerge, ntasksmerge
  
  ! Declare variables for websky parameters
  character(len=64) :: maps
  integer :: nside_map, npix_map
  real :: fov_map, zmin_map, zmax_map
  character(len=64) :: tabfile_map, tabfile_sfr
  integer :: model_map, scramble_map, center_map, chihview_map, PSZcut_map, ellmax

  ! Declare variables for cosmology parameters
  real :: ns, As, sigma8, tau, mnu
  
  ! Declare variables for ellipsoidal collapse parameters
  real :: Rsmooth_max 
  character(len=32) :: rapi
  
  ! Declare variables for lattice parameters for hpkvd
  integer :: nsub
  !integer :: n1, n2, n3
  real :: cellsize, buffersize, mlatt
  
  ! Declare variables for advanced merge parameters
  integer :: iLexc, iLmrg, iFexc, iFmrg

  ! MUSIC parameters
  integer :: levelmin, levelmax

contains

subroutine read_parameters(file_name)
  character(len=*), intent(in) :: file_name
  character(len=256) :: line, section, key, value
  integer :: ios
  logical :: in_main, in_cosmology, in_boolean_switches, in_machine_params
  logical :: in_box_params, in_parallelization_params
  logical :: in_parallelization_params_websky, in_redshifts
  logical :: in_peak_displacement, in_nongaussianities
  logical :: in_merging_algorithm, in_websky_parameters
  logical :: in_ellipsoidal_collapse, in_lattice_parameters_hpkvd
  logical :: in_advanced_merge_parameters, in_seed
  logical :: in_music_setup

  call set_default_params

  ! Initialize sections
  in_main = .false.
  in_cosmology = .false.
  in_seed = .false.
  in_boolean_switches = .false.
  in_machine_params = .false.
  in_box_params = .false.
  in_parallelization_params = .false.
  in_parallelization_params_websky = .false.
  in_redshifts = .false.
  in_peak_displacement = .false.
  in_nongaussianities = .false.
  in_merging_algorithm = .false.
  in_websky_parameters = .false.
  in_ellipsoidal_collapse = .false.
  in_lattice_parameters_hpkvd = .false.
  in_advanced_merge_parameters = .false.
  in_music_setup = .false.

  ! Open the configuration file
  open(unit=10, file=file_name, status='old', action='read', iostat=ios)
  if (ios /= 0) then
    print *, 'Error: Unable to open the file.'
    stop
  end if

  ! Read the file line by line
  do
    read(10, '(A)', iostat=ios) line
    if (ios /= 0) exit

    ! Trim leading and trailing whitespaces
    line = adjustl(line)
    line = trim(line)

    ! Skip comments and empty lines
    if (line(1:1) == '#' .or. len_trim(line) == 0) cycle

    ! Detect sections
    if (line(1:1) == '[' .and. line(len_trim(line):len_trim(line)) == ']') then
      section = line(2:len_trim(line)-1)

      ! Reset all section flags
      in_main = .false.
      in_cosmology = .false.
      in_seed = .false.
      in_boolean_switches = .false.
      in_machine_params = .false.
      in_box_params = .false.
      in_parallelization_params = .false.
      in_parallelization_params_websky = .false.
      in_redshifts = .false.
      in_peak_displacement = .false.
      in_nongaussianities = .false.
      in_merging_algorithm = .false.
      in_websky_parameters = .false.
      in_ellipsoidal_collapse = .false.
      in_lattice_parameters_hpkvd = .false.
      in_advanced_merge_parameters = .false.
      in_music_setup = .false.

      select case (trim(section))
        case ('boolean_switches')
          in_boolean_switches = .true.
        case ('peak_patch_main')
          in_main = .true.
        case ('cosmology')
          in_cosmology = .true.
        case ('seed')
          in_seed = .true.
        case ('machine_params')
          in_machine_params = .true.
        case ('box_params')
          in_box_params = .true.
        case ('parallelization_params')
          in_parallelization_params = .true.
        case ('parallelization_params_websky')
          in_parallelization_params_websky = .true.
        case ('redshifts')
          in_redshifts = .true.
        case ('peak_displacement')
          in_peak_displacement = .true.
        case ('nongaussianities')
          in_nongaussianities = .true.
        case ('merging_algorithm')
          in_merging_algorithm = .true.
        case ('websky_parameters')
          in_websky_parameters = .true.
        case ('ellipsoidal_collapse')
          in_ellipsoidal_collapse = .true.
        case ('lattice_parameters_hpkvd')
          in_lattice_parameters_hpkvd = .true.
        case ('advanced_merge_parameters')
          in_advanced_merge_parameters = .true.
        case ('setup')
          in_music_setup = .true.
      end select
      cycle
    end if

    ! Parse key-value pairs
    call parse_key_value(line, key, value)

    ! Assign values based on the current section
    if (in_main) then
      call assign_peak_patch_main(key, value)
    elseif (in_cosmology) then
      call assign_cosmology_parameters(key, value)
    elseif (in_seed) then
      call assign_seed_parameters(key, value)
    elseif (in_boolean_switches) then
      call assign_boolean_switches(key, value)
    elseif (in_machine_params) then
      call assign_machine_params(key, value)
    elseif (in_box_params) then
      call assign_box_params(key, value)
    elseif (in_parallelization_params) then
      call assign_parallelization_params(key, value)
    elseif (in_parallelization_params_websky) then
      call assign_parallelization_params_websky(key, value)
    elseif (in_redshifts) then
      call assign_redshifts(key, value)
    elseif (in_peak_displacement) then
      call assign_peak_displacement(key, value)
    elseif (in_nongaussianities) then
      call assign_nongaussianities(key, value)
    elseif (in_merging_algorithm) then
      call assign_merging_algorithm(key, value)
    elseif (in_websky_parameters) then
      call assign_websky_parameters(key, value)
    elseif (in_ellipsoidal_collapse) then
      call assign_ellipsoidal_collapse(key, value)
    elseif (in_lattice_parameters_hpkvd) then
      call assign_lattice_parameters_hpkvd(key, value)
    elseif (in_advanced_merge_parameters) then
      call assign_advanced_merge_parameters(key, value)
    elseif (in_music_setup) then
      call music_assign_setup_params(key, value)
    end if
  end do

  ! Close the file
  close(10)

  ! If -1 was passed for some parameters, evaluate them according to the formulas provided
  call evaluate_parameters

  nstepmax      = 10000
  iwant_evmap   = 3     ! 1: turn vir z vs e_v, 2: z vs p_v for fixed e_v, 3: table
  iwant_rd      = 1
  zinit_fac_ell = 20
  tfac          = 0.01  ! fraction of local 1-axis Hubble time for dt
  ihard         = 1
  H_e           = 0.0   ! Would have to fix this for non-gaussian runs
  Tkfile = 'tables/Dphi2Dzeta_Tk_LSS_units.dat'

  call add_rundir_prefix

end subroutine read_parameters

  subroutine parse_key_value(line, key, value)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key, value
    integer :: pos

    ! Find the position of the '=' character
    pos = index(line, '=')

    ! Extract key and value
    if (pos > 0) then
      key = adjustl(line(1:pos-1))
      value = adjustl(line(pos+1:len_trim(line)))
    else
      key = ''
      value = ''
    end if
  end subroutine parse_key_value

  subroutine assign_boolean_switches(key, value)
  character(len=*), intent(in) :: key, value
  integer :: int_value

  read(value, *) int_value

  select case (trim(key))
    case ('hpkvd_params')
      hpkvd_params = int_value
    case ('compile_hpkvd')
      compile_hpkvd = int_value
    case ('create_filterbank')
      create_filterbank = int_value
    case ('merge_params')
      merge_params = int_value
    case ('compile_merge')
      compile_merge = int_value
    case ('map_params')
      map_params = int_value
    case ('compile_maps')
      compile_maps = int_value
    case ('batch')
      batch = int_value
    case ('submit')
      submit = int_value
  end select
end subroutine assign_boolean_switches

subroutine assign_machine_params(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('machine')
      machine = value
    case ('submit_command')
      submit_command = value
  end select
end subroutine assign_machine_params

subroutine assign_peak_patch_main(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('seed')
      read(value, *) seed
    case ('run_name')
      run_name = value
      default_filein = run_name ! Following Python, supposed to be densefileout 
      default_fileout = 'output/' // trim(run_name) // '_raw.pksc' ! Following Python, supposed to be mergepkout
    case ('short_name')
      short_name = value
    case ('runtype')
      runtype = value
  end select
end subroutine assign_peak_patch_main

subroutine assign_box_params(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('boxsize')
      read(value, *) boxsize
    case ('nmesh')
      read(value, *) nmesh
    case ('nbuff')
      read(value, *) nbuff
    case ('ntile')
      read(value, *) ntile
    case ('largerun')
      read(value, *) largerun
  end select
end subroutine assign_box_params

subroutine assign_parallelization_params(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('tlimit')
      tlimit = value
    case ('nnodes')
      read(value, *) nnodes
    case ('tpnode')
      read(value, *) tpnode
    case ('ntasks')
      read(value, *) ntasks
    case ('ncpus')
      read(value, *) ncpus
    case ('nompth')
      read(value, *) nompth
  end select
end subroutine assign_parallelization_params

subroutine assign_parallelization_params_websky(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('tlimit_map')
      tlimit_map = value
    case ('nnodes_map')
      read(value, *) nnodes_map
    case ('ppn_map')
      read(value, *) ppn_map
    case ('ntasks_map')
      read(value, *) ntasks_map
    case ('np_map')
      read(value, *) np_map
    case ('nompth_map')
      read(value, *) nompth_map
  end select
end subroutine assign_parallelization_params_websky

subroutine assign_redshifts(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('ievol')
      read(value, *) ievol
    case ('num_redshifts')
      read(value, *) num_redshifts
    case ('maximum_redshift')
      read(value, *) maximum_redshift
    case ('global_redshift')
      read(value, *) global_redshift
  end select
end subroutine assign_redshifts

subroutine assign_peak_displacement(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('ilpt')
      read(value, *) ilpt
    case ('ioutfield')
      read(value, *) ioutfield
    case ('ireadfield')
      read(value, *) ireadfield
    case ('iwant_field_part')
      read(value, *) iwant_field_part
    case ('fielddir')
      fielddir = value
    case ('densfilein')
      densfilein = value
    case ('densfileout')
      densfileout = value
  end select
end subroutine assign_peak_displacement

subroutine assign_nongaussianities(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('NonGauss')
      read(value, *) NonGauss
    case ('fNL')
      read(value, *) fNL
    case ('A_nG')
      read(value, *) A_nG
    case ('B_nG')
      read(value, *) B_nG
    case ('R_nG')
      read(value, *) R_nG
    case ('m_phi')
      read(value, *) m_phi
    case ('m_chi')
      read(value, *) m_chi
    case ('phi_w')
      read(value, *) phi_w
    case ('phi_p')
      read(value, *) phi_p
    case ('vev')
      read(value, *) vev
    case ('m_tach')
      read(value, *) m_tach
    case ('a_e')
      read(value, *) a_e
    case ('H_e')
      read(value, *) H_e
    case ('chi_seed') ! Keep in mind the difference in naming conventions
      read(value, *) ng_seed
  end select
end subroutine assign_nongaussianities

subroutine assign_merging_algorithm(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('iZeld')
      read(value, *) iZeld
    case ('ntilemerge')
      read(value, *) ntilemerge
    case ('ntasksmerge')
      read(value, *) ntasksmerge
    case ('iwrap')
      read(value, *) iwrap
  end select
end subroutine assign_merging_algorithm

subroutine assign_websky_parameters(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('maps')
      maps = value
    case ('nside_map')
      read(value, *) nside_map
    case ('npix_map')
      read(value, *) npix_map
    case ('fov_map')
      read(value, *) fov_map
    case ('zmin_map')
      read(value, *) zmin_map
    case ('zmax_map')
      read(value, *) zmax_map
    case ('tabfile_map')
      tabfile_map = value
    case ('tabfile_sfr')
      tabfile_sfr = value
    case ('model_map')
      read(value, *) model_map
    case ('scramble_map')
      read(value, *) scramble_map
    case ('center_map')
      read(value, *) center_map
    case ('chihview_map')
      read(value, *) chihview_map
    case ('PSZcut_map')
      read(value, *) PSZcut_map
    case ('ellmax')
      read(value, *) ellmax
  end select
end subroutine assign_websky_parameters

subroutine assign_cosmology_parameters(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('Omx')
      read(value, *) Omx
    case ('OmB')
      read(value, *) OmB
    case ('Omvac')
      read(value, *) Omvac
    case ('h')
      read(value, *) h
    case ('ns')
      read(value, *) ns
    case ('As')
      read(value, *) As
    case ('sigma8')
      read(value, *) sigma8
    case ('tau')
      read(value, *) tau
    case ('mnu')
      read(value, *) mnu
    case ('pkfile')
      pkfile = value
    ! MUSIC parameters
    case ('nspec')
      read(value, *) ns
    case ('H0')
      read(value, *) h
      h = h / 100
  end select
end subroutine assign_cosmology_parameters

subroutine assign_seed_parameters(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('seed')
      read(value, *) seed
    case ('chi_seed')
      read(value, *) ng_seed
  end select
end subroutine assign_seed_parameters

subroutine assign_ellipsoidal_collapse(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('ioutshear')
      read(value, *) ioutshear
    case ('rmax2rs')
      read(value, *) rmax2rs
    case ('wsmooth')
      read(value, *) wsmooth
    case ('Rsmooth_max')
      read(value, *) Rsmooth_max
    case ('rapi')
      rapi = value
    case ('TabInterpFile')
      TabInterpFile = value
    case ('TabInterpNx')
      read(value, *) TabInterpNx
    case ('TabInterpNy')
      read(value, *) TabInterpNy
    case ('TabInterpNz')
      read(value, *) TabInterpNz
    case ('TabInterpX1')
      read(value, *) TabInterpX1
    case ('TabInterpX2')
      read(value, *) TabInterpX2
    case ('TabInterpY1')
      read(value, *) TabInterpY1
    case ('TabInterpY2')
      read(value, *) TabInterpY2
    case ('TabInterpZ1')
      read(value, *) TabInterpZ1
    case ('TabInterpZ2')
      read(value, *) TabInterpZ2
    case ('filterfile')
      filterfile = value
    case ('iforce_strat')
      read(value, *) iforce_strat
    case ('ivir_strat')
      read(value, *) ivir_strat
    case ('dcrit')
      read(value, *) dcrit
    case ('fcoll_3')
      read(value, *) fcoll_3
    case ('fcoll_2')
      read(value, *) fcoll_2
    case ('fcoll_1')
      read(value, *) fcoll_1
  end select
end subroutine assign_ellipsoidal_collapse

subroutine assign_lattice_parameters_hpkvd(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('nsub')
      read(value, *) nsub
    case ('next')
      read(value, *) next
    case ('dcore_box')
      read(value, *) dcore_box
    case ('cellsize')
      read(value, *) cellsize
    case ('buffersize')
      read(value, *) buffersize
    case ('dL_box')
      read(value, *) dL_box
    case ('mlatt')
      read(value, *) mlatt
    case ('cenx')
      read(value, *) cenx
    case ('ceny')
      read(value, *) ceny
    case ('cenz')
      read(value, *) cenz
    case ('nlx')
      read(value, *) nlx
    case ('nly')
      read(value, *) nly
    case ('nlz')
      read(value, *) nlz
    !case ('n1')
    !  read(value, *) n1
    !case ('n2')
    !  read(value, *) n2
    !case ('n3')
    !  read(value, *) n3
  end select
end subroutine assign_lattice_parameters_hpkvd

subroutine assign_advanced_merge_parameters(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('iLexc')
      read(value, *) iLexc
    case ('iLmrg')
      read(value, *) iLmrg
    case ('iFexc')
      read(value, *) iFexc
    case ('iFmrg')
      read(value, *) iFmrg
  end select
end subroutine assign_advanced_merge_parameters

subroutine music_assign_setup_params(key, value)
  character(len=*), intent(in) :: key, value

  select case (trim(key))
    case ('boxlength')
      read(value, *) boxsize
    case ('zstart')
      read(value, *) global_redshift
    case ('levelmin')
      read(value, *) levelmin
      nmesh = 2**levelmin
    case ('levelmax')
      read(value, *) levelmax
      nmesh = 2**levelmax
  end select
end subroutine music_assign_setup_params

subroutine evaluate_parameters()

  implicit none
  
  real, parameter :: epsilon = 1.0e-12
  
  ! Handle algebraic operations
  if (ntasks == -1) ntasks = nnodes * tpnode
  if (ntasks_map == -1) ntasks_map = ppn_map * nnodes_map
  if (trim(densfilein) == '-1') densfilein = run_name
  if (trim(densfileout) == '-1' ) densfileout = run_name
  if (trim(filein) == '-1' ) filein = default_filein
  if (trim(fileout) == '-1' ) fileout = default_fileout
  if (iZeld == -1) iZeld = ilpt
  if (ntilemerge == -1) ntilemerge = ntile
  if (ntasksmerge == -1) ntasksmerge = ntasks
  if (nsub == -1) nsub = nmesh - 2 * nbuff
  if (next == -1) next = nsub * ntile + 2 * nbuff
  if (abs(Omvac + 1.0) < epsilon) Omvac = 1.0 - Omx - OmB
  !if (abs(Omvac + Omx - 1.0) > epsilon) OmB = 1.0 - Omx - Omvac ! Required for MUSIC
  if (abs(dcore_box + 1.0) < epsilon) dcore_box = boxsize / ntile
  if (abs(cellsize + 1.0) < epsilon) cellsize = dcore_box / nsub
  if (abs(buffersize + 1.0) < epsilon) buffersize = cellsize * nbuff
  if (abs(dL_box + 1.0) < epsilon) dL_box = dcore_box + 2 * buffersize
  if (abs(mlatt + 1.0) < epsilon) mlatt = 2.775e11 * (Omx + OmB) * h**2 * (cellsize)**3
  if (nlx == -1) nlx = ntile
  if (nly == -1) nly = ntile
  if (nlz == -1) nlz = ntile
  !if (n1 == -1) n1 = nmesh
  !if (n2 == -1) n2 = nmesh
  !if (n3 == -1) n3 = nmesh
end subroutine evaluate_parameters

subroutine add_rundir_prefix

  run_dir = 'RUNDIR'
  fielddir = trim(run_dir) // '/' // trim(fielddir)
  pkfile = trim(run_dir) // '/tables/' // trim(pkfile)
  Tkfile = trim(run_dir) // '/' // trim(Tkfile)
  filterfile = trim(run_dir) // '/' // trim(filterfile)
  fileout = trim(run_dir) // '/' // trim(fileout)
  TabInterpFile = trim(run_dir) // '/' // trim(TabInterpFile)

end subroutine add_rundir_prefix

subroutine set_default_params
  
  ! By default, these parameters will be calculated in evaluate_parameters subroutine
  ntasks = -1
  ntasks_map = -1
  densfilein = '-1'
  densfileout = '-1'
  filein = '-1'
  fileout = '-1'
  iZeld = -1
  ntilemerge = -1
  ntasksmerge = -1
  nsub = -1
  next = -1
  Omvac = -1.0
  dcore_box = -1.0
  cellsize = -1.0
  buffersize = -1.0
  dL_box = -1.0
  mlatt = -1.0
  nlx = -1
  nly = -1
  nlz = -1

  ! Set default parameters if they aren't declared
  seed = 13579
  ng_seed = 90101
  
end subroutine set_default_params

end module config_reader
