module io

  use iso_c_binding

  implicit none

  interface

     ! Creates a Fortran wrapper for the C function myfwrite in src/modules/External/myio.C
     subroutine myfwrite(filename,buffer,size,offset) &
         bind(c,name="myfwrite") 

       use iso_c_binding

       character(C_CHAR), intent(in) :: filename
       type(C_PTR), value            :: buffer
       integer(C_LONG), value        :: size, offset

     end subroutine myfwrite

   end interface

end module io
