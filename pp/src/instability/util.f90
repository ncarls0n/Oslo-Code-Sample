! util.f90

! Module for some utility subroutines when dealing with conversion between matricies and arrays.

module util_mod
  use params
  
  implicit none

contains
  
  subroutine set_ltri_mat(y, mat, dim)
    real(dl), intent(in) :: y(:)
    real(dl), intent(out) :: mat(:,:)
    integer, intent(in) :: dim

    integer :: i,j,k

    mat(:,:) = 0._dl
    k = 0
    do j=1,dim
       do i=j,dim
          k = k+1
          mat(i,j) = y(k)
       end do
    end do
  end subroutine set_ltri_mat

  subroutine set_sym_mat(y, mat, dim)
    real(dl), intent(in) :: y(:)
    real(dl), intent(out) :: mat(:,:)
    integer, intent(in) :: dim

    integer :: i,j,k

    mat(:,:) = 0._dl
    k = 0
    do j=1,dim
       do i=j,dim
          k = k+1
          mat(i,j) = y(k)
          mat(j,i) = mat(i,j)
       end do
    end do
  end subroutine set_sym_mat

  subroutine set_asym_mat(y, mat, dim)
    real(dl), intent(in) :: y(:)
    real(dl), intent(out) :: mat(:,:)
    integer, intent(in) :: dim

    integer :: i,j,k

    mat(:,:) = 0._dl
    k = 0
    do j=1,dim-1
       do i=j+1,dim
          k = k+1
          mat(i,j) = y(k)
          mat(j,i) = -mat(i,j)
       end do
    end do
  end subroutine set_asym_mat

  ! Set an array from the lower triangular part of a matrix.
  ! subdiag specifies how many (sub)diagonals should be skipped.
  subroutine set_lower_ar(y, mat, dim, subdiag)
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: mat(:,:)
    integer, intent(in) :: dim, subdiag

    integer :: i,j,k

    k = 0
    do j=1,dim-subdiag
       do i=j+subdiag,dim
          k = k+1
          y(k) = mat(i,j)
       end do
    end do
  end subroutine set_lower_ar

  ! Subroutine to symmetrize a matrix
  subroutine symmetrize_mat(mat, temp_mat, dim)
    real(dl) :: mat(:,:), temp_mat(:,:)
    integer, intent(in) :: dim

    integer :: i,j

    do j=1,dim
       do i=1,dim
          temp_mat(i,j) = mat(j,i)
       end do
    end do
    temp_mat = mat + temp_mat
  end subroutine symmetrize_mat
  
end module util_mod
