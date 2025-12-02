program hpkvdmain
  ! You probably expected to see a whole hpkvd logic in here?
  ! For modularity, all the main hpkvd code was moved to src/hpkvd/hpkvdmodule.f90
  ! This program is now just a wrapper around hpkvdmodule, which produces an executable

  use mpivars       ! from src/modules/External
  use hpkvdmodule   ! from src/hpkvd/hpkvdmodule
  use textlib       ! from src/modules/External, procedures: i4arg, s256arg, r4arg, r8arg, indlnb, getu, onedot, s2i4, s2r8
  
  implicit none
  
  integer :: maketable_flag, seedFFT
  character(len=512) :: argument  ! Declare a variable to hold the arguments
  character(len=256) :: parameters_filename  ! Declare a variable to hold the name of parameters file (.ini)
  real, pointer :: input_field(:)  ! Declare input_field as a pointer

  ! Initialize the MPI environment
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call mpi_comm_size(mpi_comm_world, ntasks, ierr)

  ! Run hpkvd module with passing command line arguments into it
  call run_hpkvd()
  
  ! Synchronize MPI processes and finalize MPI environment
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)

contains

subroutine run_hpkvd

  ! Get command line arguments using i4arg (4-byte integer argument)
  maketable_flag = i4arg(1, 0)
  seedFFT  = i4arg(2, 13579)
  parameters_filename  = s256arg(3, "hpkvd_params.bin")
  
  ! ADDING THESE LINES ALLOWS FOR THE CLI ARGS TO PASS CORRECTLY.
  ! DELETING THE LINES BREAKS THIS. I HAVE NO IDEA WHY THIS IS HAPPENING
  call get_command_argument(1, argument)
  call get_command_argument(2, argument)
  call get_command_argument(3, argument)

  ! Pass input_field as null(), maketable_flag, and seedFFT to the hpkvd subroutine
  call hpkvd(parameters_filename, null(), maketable_flag, seedFFT)

end subroutine run_hpkvd

end program hpkvdmain
