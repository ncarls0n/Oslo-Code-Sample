module exclusion

  !-----------------------------------------------------------------------
  ! Routines for Lagrangian exclusion
  !-----------------------------------------------------------------------
  use mpivars
  integer, allocatable :: hoc(:,:,:)

contains 

  subroutine pairwise(pairwise_code,x,y,z,r,vx1,vy1,vz1,vx2,vy2,vz2,harray,nharray,nin,nout,xbox,ybox,zbox,dL_box)

!    implicit none
    integer nin, nout, nharray, nc, i, m, ix,iy,iz,imx,jmx,kmx,ll(nin),id(nin)
    real    x(nin),y(nin),z(nin),vx1(nin),vy1(nin),vz1(nin),vx2(nin),vy2(nin),vz2(nin)
    real    r(nin),harray(nharray,nin),buffer(nin),f(nin),dV(nin)
    logical survived(nin), merged
    real    x1,x2,rmax,xmax
    real    xbox,ybox,zbox,dL_box,xboxmin,yboxmin,zboxmin
    
    parameter (nbr_mx=100000)
    real xlcl(nbr_mx),ylcl(nbr_mx),zlcl(nbr_mx),rlcl(nbr_mx)
    integer idlcl(nbr_mx)

    ! Edit from April 2024
    CHARACTER(len=7) :: reduction

    character *512 pairwise_code

    pi  = 3.14159
    ftp = 4./3.*pi

    nc = 256

    xboxmin = xbox-dL_box/2
    yboxmin = ybox-dL_box/2
    zboxmin = zbox-dL_box/2

    if(.not.allocated(hoc)) allocate(hoc(0:nc-1,0:nc-1,0:nc-1))

    do m=1,nin
       id(m) = m
    enddo

    buffer = r
    call sort(nin,-r,id)
    r = buffer

    survived = .true.

    ! BUILD LINKED LISTS

    ! CONVERT TO CELL SIZE LENGTH UNITS
    x = (x - xboxmin)/dL_box * nc
    y = (y - yboxmin)/dL_box * nc
    z = (z - zboxmin)/dL_box * nc
    r = (r     )/dL_box * nc 

    dV   = 0.
    hoc  = 0
    ll   = 0
    rmax = -1e10
    do mi=1,nin
       idi = id(mi)
       ix  = int(x(idi))
       iy  = int(y(idi))
       iz  = int(z(idi))
       ix  = min(ix,nc-1)
       iy  = min(iy,nc-1)
       iz  = min(iz,nc-1)

       ll(mi)        = hoc(ix,iy,iz)
       hoc(ix,iy,iz) = mi
    enddo


    ! LOOP OVER PEAKS, BEGINNING WITH MOST MASSIVE
    do mi=1,nin

!       if(myid==0 .and. mod(mi,10000)==0) write(*,*) mi,nin,(100.*mi)/nin,'%'

       idi = id(mi)
       fi  = 0.0

       ! DON'T DO ANYTHING IF THIS PEAK HAS ALREADY BEEN ANNIHILATED
       if(.not.survived(idi)) cycle

       xi  = x(idi)
       yi  = y(idi)
       zi  = z(idi)
       ri  = r(idi)

       i1 = max(0,   int(xi)-2*int(ri))
       i2 = min(nc-1,int(xi)+2*int(ri))

       j1 = max(0,   int(yi)-2*int(ri))
       j2 = min(nc-1,int(yi)+2*int(ri))

       k1 = max(0,   int(zi)-2*int(ri))
       k2 = min(nc-1,int(zi)+2*int(ri))

       nlcl = 0

       xlcl  = 0
       ylcl  = 0
       zlcl  = 0
       rlcl  = 0
       idlcl = 0

       ! FIRST MAKE LIST OF NEIGHBORS
       do i=i1,i2
          do j=j1,j2
             do k=k1,k2

                mj  = hoc(i,j,k)
                do while(mj>0)

                   idj = id(mj)
                   xj  =  x(idj)
                   yj  =  y(idj)
                   zj  =  z(idj)
                   rj  =  r(idj)
                   d = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
       
                   if (d < (ri + rj) .and. survived(idj) .and. ri >= rj .and. &
                        idj /= idi) then ! peaks i and j overlap
                      nlcl = nlcl + 1
                      if(nlcl > nbr_mx) then
                         write(*,*) 'nlcl > nbr_max', nlcl, nbr_mx
                         call mpi_finalize(ierr)
                         stop
                      endif
                      xlcl(nlcl)  = xj
                      ylcl(nlcl)  = yj
                      zlcl(nlcl)  = zj
                      rlcl(nlcl)  = rj
                      idlcl(nlcl) = idj

                   endif
                   mj = ll(mj)

                enddo

             enddo
          enddo
       enddo

       ! SORT NEIGHBORS BY MASS

       call sort(nlcl,-rlcl,idlcl)
       
       ! NOW GO DOWN THE LIST AND DO THE SPECIFIED PAIRWISE LOGIC
       do mj=1,nlcl

          idj = idlcl(mj)
          xj  = x(idj)
          yj  = y(idj)
          zj  = z(idj)
          rj  = r(idj)
          vi  = 0.
          vj  = 0.

          if(pairwise_code=='exclusion') call pairwise_exclusion(xi,yi,zi,ri,xj,yj,zj,rj,fi,merged)
          if(pairwise_code=='reduction') call pairwise_reduction(xi,yi,zi,ri,xj,yj,zj,rj,vi,vj,fi,merged)
          
          if(merged) survived(idj) = .false.

          reduction='shared'

          if(reduction=='shared')  then
             dV(idi) = dV(idi) + vi
             dV(idj) = dV(idj) + vj
          endif

          if(reduction=='smaller') then
             dV(idj) = dV(idj) + vi + vj
          endif

          x(idj)  = xj
          y(idj)  = yj
          z(idj)  = zj
          r(idj)  = rj


       enddo

       ! ASSIGN NEW VALUES TO PEAK i (displacements not re-calculated yet)

       x(idi) = xi
       y(idi) = yi
       z(idi) = zi

    enddo

    if(pairwise_code=='reduction') then 
       do mi=1,nin
          idi = id(mi)
          ! DON'T DO ANYTHING IF THIS PEAK HAS ALREADY BEEN ANNIHILATED
          if(.not.survived(idi)) cycle
          
          ! Subtract total overlapping volume 
          r(idi) = (r(idi)**3 - dV(idi)/ftp)**(1./3)

          if(r(idi)<0.) survived(idi) = .false.
       enddo
    endif

    ! FINAL CLEANUP
    nout = 0
    do i = 1,nin
       if(survived(i)) then
          nout        = nout + 1

          x(nout)     = x(i)
          y(nout)     = y(i)
          z(nout)     = z(i)
          r(nout)     = r(i)

          vx1(nout)   = vx1(i)
          vy1(nout)   = vy1(i)
          vz1(nout)   = vz1(i)
          vx2(nout)   = vx2(i)
          vy2(nout)   = vy2(i)
          vz2(nout)   = vz2(i)

          harray(:,nout) = harray(:,i) 
       endif
    enddo

    ! CHANGE BACK TO PHYSICAL COORDINATES
    x = xboxmin + dL_box / nc * x
    y = yboxmin + dL_box / nc * y
    z = zboxmin + dL_box / nc * z
    r =           dL_box / nc * r
       
    return

  end subroutine pairwise

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sort(nin,x,id)

    real x(nin),fid(nin)
    integer id(nin)

    fid = real(id)

    call sort2(nin,x,fid)

    id  = int(fid)

    return

  end subroutine sort

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pairwise_exclusion(xi,yi,zi,ri,xj,yj,zj,rj,f,merged)

    logical merged
    merged = .false.

    dij = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

    if(dij < ri    ) merged=.true. ! Kill halo if inside another 

    return

  end subroutine pairwise_exclusion

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pairwise_reduction(xi,yi,zi,ri,xj,yj,zj,rj,vi,vj,f,merged)

    logical merged
    merged = .false.

    pi = 3.14159
    ftp = 4./3.*pi

    dij = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
    eta = 1.0

    if(dij < ri+rj ) then
       call sphere_overlap(dij,ri,rj,vi,vj)

!       rj = (rj**3 - (v1 + v2)/ftp)**(1./3) ! Subtract overlapping volume from 
                                            ! the radius of the smaller of the 2 halos
!       dij_f = ri+rj
!       if(dij_f > dij) then 
          ! reduced peaks still overlap, as subtraction was in mass, not radius
          ! By 'conservation of mass' we must recenter if they still overlap
!          xj = xj + (xj-xi)/dij * (dij_f - dij)
!          yj = yj + (yj-yi)/dij * (dij_f - dij)
!          zj = zj + (zj-zi)/dij * (dij_f - dij)        
!       endif
    endif

!    if(dij**2 > ri**2 - rj**2) f = max(rj**3 / ri**3 * eta , f) 

    return

  end subroutine pairwise_reduction

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sphere_overlap(d,r1,r2,v1,v2)
    ! http://mathworld.wolfram.com/Sphere-SphereIntersection.html

    pi = 3.14159

    h1 = (r2 - r1 + d)*(r2 + r1 - d)/2./d
    h2 = (r1 - r2 + d)*(r1 + r2 - d)/2./d

    v1 = 1./3.*pi*h1**2*(3*r1-h1)
    v2 = 1./3.*pi*h2**2*(3*r2-h2)

    return 

  end subroutine sphere_overlap
         

end module exclusion
