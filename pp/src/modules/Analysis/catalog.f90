module peakdata

  implicit none

  integer nhalo

  ! Halo arrays
  real, allocatable :: pos(:,:),rvir(:),rth(:),mass(:)

contains

  !=======================================================================

  subroutine loadpeaks(filename,form)
    
    implicit none
    
    integer i,j,idum,m

    real dum

    character *512 filename, form

    if(form='asci') then

       open(unit=1,file=filename)
       nhalo=0
10     continue
         read(1,*,end=20) dum
         nhalo = nhalo + 1
       goto 10
20     continue
       close(1)

       allocate(pos(3,nhalo))
       allocate(mass(nhalo))

       open(unit=1,file=filename)
       read(1,*) ( ( (pos(j,i),j=1,3), mass(i) ), i=1,nhalo)
       close(1)

    endif
       
    return

  end subroutine loadpeaks

end module peakdata

