module hpkvd_c_wrapper
    use, intrinsic :: iso_c_binding
    use mpivars
    use hpkvdmodule
    use textlib
    implicit none

contains

    subroutine c_run_hpkvd(parameters_filename, input_field, maketable_flag, seedFFT) bind(C, name="c_run_hpkvd")
        implicit none

        integer(c_int), intent(in) :: maketable_flag, seedFFT
        real(c_float), pointer, intent(in), optional :: input_field(:,:,:) ! Random field if passed as a pointer
        character(kind=c_char), intent(in) :: parameters_filename

        ! Local variables
        integer :: ierr, myid, ntasks

        ! Initialize the MPI environment
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world, myid, ierr)
        call mpi_comm_size(mpi_comm_world, ntasks, ierr)

        ! Debug print statements
        print *, "MPI initialized with rank ", myid, " out of ", ntasks
        print *, "maketable_flag = ", maketable_flag, "seedFFT = ", seedFFT

        ! Pass input_field (or null if not present), maketable_flag, and seedFFT to the hpkvd subroutine
        if (present(input_field)) then
            print *, "input_field is present"
            call hpkvd(parameters_filename, input_field, maketable_flag, seedFFT)
        else
            print *, "input_field is not present"
            call hpkvd(parameters_filename, null(), maketable_flag, seedFFT)
        end if

        ! Synchronize MPI processes and finalize MPI environment
        call mpi_barrier(mpi_comm_world, ierr)
        call mpi_finalize(ierr)
        print *, "MPI finalized"

    end subroutine c_run_hpkvd

end module hpkvd_c_wrapper
