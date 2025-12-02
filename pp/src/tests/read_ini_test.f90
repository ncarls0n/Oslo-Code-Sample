program main
  use config_reader
  !use print_module
  implicit none

  ! Call the subroutine to read parameters
  call read_parameters('parameters.ini')

  ! Print values to verify using the print_parameter subroutine
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
  call print_parameter_int('  chi_seed:', chi_seed)

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

contains
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

    write(str, '(F0.3)') value
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
end program main
