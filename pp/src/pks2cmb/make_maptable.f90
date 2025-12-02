program make_maptable

  use maptable
  use cosmology
  use textlib

  implicit none

  ! So far 3 models implemented
  !    model = 1 --> BPPS2 AGN Feedback Delta = 200
  !    model = 2 --> BPPS2 AGN Feedback Delta = 500
  !    model = 3 --> BPPS2    Adiabatic Delta = 500

  integer model

  ! Output filename
  character *512 fileout

  ! Usage check
  if(command_argument_count()/=2) then
     write(*,11) 
11   format('Usage: make_maptable <fileout> <model>')
     stop
  endif

  ! Get commandline
  call get_command_argument(1,fileout)
  model = i4arg(2,1)

  ! Set background cosmology
  omegam = 0.25
  omegab = 0.043
  omegal = 0.75
  h      = 0.7
  sigma8 = 0.8
  ns     = 0.96
  w      = -1
  fb     = omegab/omegam

  rho_0     = 2.775e11*omegam*h**2
  rhocrit_0 = 2.775e11*h**2  

  ! Make table 
  call makemaptable(fileout,model)

  stop
  
end program make_maptable
