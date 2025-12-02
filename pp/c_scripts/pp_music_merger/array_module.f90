module array_module
  use iso_c_binding
  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: n = 3
  real(dp), dimension(n,n,n), target :: array

contains

  subroutine initialize_array() bind(C)
    integer :: i, j, k
    do k = 1, n
      do j = 1, n
        do i = 1, n
          array(i, j, k) = i + (j - 1) * n + (k - 1) * n * n
        end do
      end do
    end do
  end subroutine initialize_array

  subroutine get_array_pointer(ptr) bind(C)
    use iso_c_binding, only: c_ptr
    type(c_ptr), intent(out) :: ptr
    ptr = c_loc(array)
  end subroutine get_array_pointer

end module array_module

