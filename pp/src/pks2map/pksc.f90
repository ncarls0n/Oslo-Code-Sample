module pksc  

  use cosmology
  use mpivars 

  implicit none

  real delta_pksc
  parameter(delta_pksc=200)

  ! Halo arrays
  real, allocatable :: posxyz(:,:),velxyz(:,:),rth(:),mass(:),vrad(:)
  real rhi
  integer nhalo, nhalol 
  real RTHmaxtot
  integer offset_num_floats

contains

  subroutine loadpksc(filename)
  !=======================================================================
  !LOAD IN HALOS
  !=======================================================================

    use profiles
    
    implicit none
    
    integer i,j,idum,m,nhalovars
    real dum
    
    character *512 filename

    if(testcase) then ! default testcase=.false.
       nhalo = 1
       allocate(posxyz(3,nhalo))
       allocate(rth(nhalo))
       allocate(mass(nhalo))
       posxyz(1,1) = 0
       posxyz(2,1) = 0
       posxyz(3,1) = rofz(0.05)
       return
    endif

    ! Open catalogue file

    ! Edit from April 2024 - trying to make it work on gcc and following modern standard
    ! (form='binary' was depreciated)
    !open(unit=1, file=filename, form='binary')
    open(unit=1, file=filename, form='unformatted', access='stream', status='replace')

    read(1) nhalo, RTHmaxtot ! Read header
    close(1)

    ! Determine the number of variables per halo in the catalogue: read
    ! catalogue size in bytes, divide by 4 to get total number of 32-bit
    ! numbers, subtract 3, the length of the header, and divide by the
    ! number of halos in the catalogue
    inquire(file=filename,size=nhalovars)
    nhalovars=(nhalovars/4-3)/nhalo

    if (myid==0) write(*,*) "Nhalo total = ", nhalo
    !Read in only the appropriate nhalo/ntasks, not whole catalogue
    nhalol = int((nhalo-1)/ntasks)+1          ! even number per processor
    nhalo  = min(nhalol, nhalo - nhalol*myid) ! last processor may not have 
                                              ! same number of halos

    if (myid==0) write(*,*) "Nhalo per processor = ", nhalo
    allocate(posxyz(3,nhalo))
    allocate(velxyz(3,nhalo))
    allocate(rth(nhalo))
    allocate(mass(nhalo))
    allocate(vrad(nhalo))

    offset_num_floats = nhalovars*nhalo*myid

    ! Edit from April 2024 - trying to make it work on gcc and following modern standard
    ! (form='binary' was depreciated)
    !open(unit=1, file=filename, form='binary')
    open(unit=1, file=filename, form='unformatted', access='stream', status='replace')

    read(1) idum, idum, idum ! read header into dummy variables
    read(1) (idum,i=1,offset_num_floats) !read into dummy all halos previous
    read(1) ( (posxyz(j,i),j=1,3), &
              (velxyz(j,i),j=1,3), rth(i),&
              (dum,j=8,nhalovars), i=1,nhalo)
    close(1)

    vrad = (posxyz(1,:)*velxyz(1,:) + posxyz(2,:)*velxyz(2,:)+ posxyz(3,:)*velxyz(3,:))
    vrad = vrad/(posxyz(1,:)**2 + posxyz(2,:)**2 + posxyz(3,:)**2)**(1./2)
    vrad = vrad / 299792.458 ! halo radial velocity in units of v/c 
    
    return

  end subroutine loadpksc


  subroutine scramble_halos
  !=======================================================================
  !Scramble halos in theta and phi
  !=======================================================================
    use random

    implicit none

    ! For random number generation   
    integer,dimension( IRandNumSize ) :: seed
    integer i

    real(kind=8) :: randnum
    real xh,yh,zh,rh,chih
    real muran, phiran

    seed = 13579 * (myid+1)
    do i=1,nhalo
        xh   = posxyz(1,i)
        yh   = posxyz(2,i)
        zh   = posxyz(3,i)
        chih = sqrt(xh**2+yh**2+zh**2)

        call ranf(seed,randnum)
        muran  = 2 * randnum - 1
        call ranf(seed,randnum)
        phiran = 2 * 3.14159 * randnum

        rh = sqrt(1-muran**2) * chih

        xh = rh    * cos(phiran)
        yh = rh    * sin(phiran)
        zh = muran * chih

        posxyz(1,i) = xh
        posxyz(2,i) = yh
        posxyz(3,i) = zh
     enddo

  end subroutine scramble_halos

  subroutine center_halos
  !=======================================================================
  !Center halos on the most massive object in the map (x,0,0)
  !=======================================================================

    implicit none

    ! For random number generation   
    integer i, imaxval(1)

    real xh,yh,zh,rh,chih
    real xmmax, xmmaxl, ymmax, ymmaxl, zmmax, zmmaxl, rmmax
    real RTHmax, RTHmaxl, rmaxtheta, rmaxphi

    !find largest remaining halo across all processors
    RTHmaxl = maxval(rth(1:nhalo))
    imaxval = maxloc(rth(1:nhalo))
    !write(*,*) "DEBUG maxloc, maxval RTH", imaxval, rth(imaxval(1))
    call mpi_allreduce(RTHmaxl,RTHmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    !Distribute position of largest halo to all processors
    xmmaxl = -1e5
    ymmaxl = -1e5
    zmmaxl = -1e5
    if(RTHmaxl == RTHmax) then
       xmmaxl = posxyz(1,imaxval(1))
       ymmaxl = posxyz(2,imaxval(1))
       zmmaxl = posxyz(3,imaxval(1))

       rmmax = sqrt(xmmaxl**2+ymmaxl**2+zmmaxl**2)
       write(*,19) RTHmax,rmmax, xmmaxl, ymmaxl, zmmaxl
19     format(/,'Before rotation largest halo is at',/,&
            'RTH, distance:    ',2(1pe9.3,1x),/,&
            'x, y, z:          ',3(1e10.3,1x),/)
       
    endif
    call mpi_allreduce(xmmaxl,xmmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    call mpi_allreduce(ymmaxl,ymmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    call mpi_allreduce(zmmaxl,zmmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)


    !rotate so largest halo is in first quadrant to make angles easier   
    if(xmmax <0.0)   posxyz(1,:) = -posxyz(1,:)
    if(ymmax <0.0)   posxyz(2,:) = -posxyz(2,:)
    if(zmmax <0.0)   posxyz(3,:) = -posxyz(3,:)
    xmmax = abs(xmmax)
    ymmax = abs(ymmax)
    zmmax = abs(zmmax)

    ! angle to rotate by in z-x plane
    rmaxtheta = asin(zmmax/sqrt(xmmax**2+zmmax**2))  
    xmmaxl    = xmmax/abs(xmmax)*sqrt(xmmax**2+zmmax**2)
    ! angle to rotate by in z-y plane
    rmaxphi   = asin(ymmax/sqrt(ymmax**2+xmmaxl**2))  

    ! Rotate by angles defined by largest halo                          
    do i=1,nhalo
       xh =   posxyz(1,i)*cos(rmaxtheta) + posxyz(3,i)*sin(rmaxtheta)
       zh =  -posxyz(1,i)*sin(rmaxtheta) + posxyz(3,i)*cos(rmaxtheta)

       yh =  -xh*sin(rmaxphi) + posxyz(2,i)*cos(rmaxphi)
       xh =   xh*cos(rmaxphi) + posxyz(2,i)*sin(rmaxphi)

       posxyz(1,i) = zh 
       posxyz(2,i) = yh
       posxyz(3,i) = xh

    enddo

    ! Double check rotation worked and find new distance to largest halo
    if(RTHmaxl == RTHmax) then
       xmmax = posxyz(3,imaxval(1))
       ymmax = posxyz(2,imaxval(1))
       zmmax = posxyz(1,imaxval(1))
       rmmax = sqrt(xmmax**2+ymmax**2+zmmax**2)
       write(*,20) RTHmax,rmmax, xmmax, ymmax, zmmax
20     format(/,'After rotation largest halo is at',/,&
            'RTH, distance:    ',2(1pe9.3,1x),/,&
            'x, y, z:          ',3(1e10.3,1x),/)
    endif

  end subroutine center_halos


end module pksc

