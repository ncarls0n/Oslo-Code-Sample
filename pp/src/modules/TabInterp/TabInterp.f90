module TabInterp

  !=======================================================================
  !
  ! Tables for Interpolation
  !
  !                                                AUTHOR: MARCELO ALVAREZ
  !                                                LAST EDIT:     06.14.16
  !
  !=======================================================================

  ! This is a generic r3 --> r1 mapping 
  abstract interface
     function r3_to_r1(x,y,z)
       real :: r3_to_r1
       real, intent(in) :: x,y,z
     end function r3_to_r1
  end interface

  ! Initialize a pointer for aliasing the name of the function
  ! to use to make the table
!  procedure (r3_to_r1), pointer, intent(in) :: TabInterpFunc => null()

  ! Array for interpolation
  real, allocatable :: TabInterpArray(:,:,:)

  ! Table file 
  character*512 :: TabInterpFile

  ! Array dimensions
  integer :: TabInterpNx, TabInterpNy, TabInterpNz

  ! Array bounds and spacings
  real    :: TabInterpX1,TabInterpX2,TabInterpY1,TabInterpY2,TabInterpZ1,TabInterpZ2
  real    :: TabInterpDx,TabInterpDy,TabInterpDz

  ! Switches for I/O
  logical :: TabInterpWrite, TabInterpRead

  ! Out of bounds value
  real    :: TabInterpOut

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine TabInterpInit()

    implicit none
    
    allocate(TabInterpArray(TabInterpNx, TabInterpNy, TabInterpNz))
    TabInterpArray = 0

  end subroutine TabInterpInit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine TabInterpMake(TabInterpFunc)

    use mpivars
    
    implicit none
    integer i,j,k,dk
    real x,y,z

    real TabInterpFunc
    external TabInterpFunc
    
    TabInterpDx = (TabInterpX2-TabInterpX1)  / (TabInterpNx-1)
    TabInterpDy = (TabInterpY2-TabInterpY1)  / (TabInterpNy-1)
    TabInterpDz = (TabInterpZ2-TabInterpZ1)  / (TabInterpNz-1)

    do k=1,TabInterpNz
       write(*,*) k,'/',TabInterpNz
       z = TabInterpZ1 + (k-1) * TabInterpDz
       do j=1,TabInterpNy
          y = TabInterpY1 + (j-1) * TabInterpDy
          do i=1,TabInterpNx
             x = TabInterpX1 + (i-1) * TabInterpDx
             TabInterpArray(i,j,k) = TabInterpFunc(x,y,z)
          enddo
       enddo
    enddo

    if(TabInterpWrite) then
       open(unit=1,file=TabInterpFile,access='stream')
       write(1) TabInterpNx, TabInterpNy, TabInterpNz,TabInterpX1,TabInterpX2,&
            TabInterpY1,TabInterpY2,TabInterpZ1,TabInterpZ2,TabInterpArray
       close(1)
    endif

    return

  end subroutine TabInterpMake

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine TabInterpLoad()

    use mpivars

    implicit none
    
    logical :: FileExists

    inquire(file=TabInterpFile,exist=FileExists)

    if(.not.FileExists) then
       if(myid==0) then
          write(*,*) '  ********************************************'
          write(*,*) '  Error: file ',trim(TabInterpFile),' does not exist'
          write(*,*) 
          write(*,*) '  Try creating it first with '
          write(*,*) '      mpirun -np 1 ./hpkvd 1 '
          write(*,*) '  Exiting ...'
          write(*,*) '  ********************************************'
          write(*,*)
       endif
       call mpi_finalize(ierr)
       stop
    endif

    open(unit=1,file=TabInterpFile,access='stream')
    read(1) TabInterpNx, TabInterpNy, TabInterpNz,TabInterpX1,TabInterpX2,&
         TabInterpY1,TabInterpY2,TabInterpZ1,TabInterpZ2,TabInterpArray
    close(1)

    TabInterpDx = (TabInterpX2-TabInterpX1)  / (TabInterpNx-1)
    TabInterpDy = (TabInterpY2-TabInterpY1)  / (TabInterpNy-1)
    TabInterpDz = (TabInterpZ2-TabInterpZ1)  / (TabInterpNz-1)

    return

  end subroutine TabInterpLoad


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real function TabInterpInterpolate(x,y,z)

    implicit none
    integer i1,j1,k1,i2,j2,k2
    real fx,fy,fz
    real x,y,z
    
    TabInterpInterpolate = TabInterpOut ! TabInterpOut=-1, assigned in src/hpkvd/hpkvd.f90
    if(x>TabInterpX2 .or. x<TabInterpX1 .or. &
       y>TabInterpY2 .or. y<TabInterpY1 .or. &
       z>TabInterpZ2 .or. z<TabInterpZ1) return ! returns -1 if x,y,z outside bounds of interplation table

    ! Determine points (i1,j1,k1) & (i2,j2,k2) such that x lies between
    ! i1,i2 & y between j1,j2 & z between k1,k2 in TabInterpArray(i,j,k)
    TabInterpInterpolate = 0
    i1 = int ( (x - TabInterpX1) / TabInterpDx ) + 1
    i2 = i1 + 1
    j1 = int ( (y - TabInterpY1) / TabInterpDy ) + 1
    j2 = j1 + 1
    k1 = int ( (z - TabInterpZ1) / TabInterpDz ) + 1
    k2 = k1 + 1

    ! Interpolates value for TabInterpArray(x,y,z)
    fx = (x - (i1-1) * TabInterpDx - TabInterpX1) / TabInterpDx
    fy = (y - (j1-1) * TabInterpDy - TabInterpY1) / TabInterpDy
    fz = (z - (k1-1) * TabInterpDz - TabInterpZ1) / TabInterpDz
    TabInterpInterpolate = TabInterpInterpolate + & 
         TabInterpArray(i1,j1,k1) * (1-fx) * (1-fy) * (1-fz) + & 
         TabInterpArray(i2,j1,k1) * (  fx) * (1-fy) * (1-fz) + & 
         TabInterpArray(i1,j2,k1) * (1-fx) * (  fy) * (1-fz) + & 
         TabInterpArray(i2,j2,k1) * (  fx) * (  fy) * (1-fz) + & 
         TabInterpArray(i1,j1,k2) * (1-fx) * (1-fy) * (  fz) + & 
         TabInterpArray(i2,j1,k2) * (  fx) * (1-fy) * (  fz) + & 
         TabInterpArray(i1,j2,k2) * (1-fx) * (  fy) * (  fz) + & 
         TabInterpArray(i2,j2,k2) * (  fx) * (  fy) * (  fz)
         
    return ! returns interpolated value

  end function TabInterpInterpolate

end module TabInterp


