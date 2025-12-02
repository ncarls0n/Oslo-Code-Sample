program test_args
  use textlib
  implicit none
  integer :: maketable_flag, seedFFT
  maketable_flag = i4arg(1, 0)
  seedFFT  = i4arg(2, 13579)
  write(*,*) 'maketable: ', maketable_flag
  write(*,*) 'seedFFT: ', seedFFT
end program test_args
