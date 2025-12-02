program merge

  use mpivars

  implicit none 

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  call  merge_pkvd

  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_finalize(ierr)

end program merge

subroutine merge_pkvd

  use textlib
  use mpivars
  use arrays_params
  use exclusion

  ! DOMAIN DECOMPOSITION VARIABLES
  real    Rmin,  Rmax, Rmax_in
  real    Rminl, Rmaxl, Volhl, Volh, fField, fHalo, pi, ftpi, Mfield_p
  real, parameter :: epsilon=1e-6

  ! DATE AND TIME
  character(8)          :: dnt_date
  character(10)         :: dnt_time
  character(5)          :: dnt_zone
  integer, dimension(8) :: dnt_values
  integer initial_time, final_time, count_rate
  integer elapsed_hours, elapsed_minutes, elapsed_seconds
  double precision elapsed_time
  integer seedin
  character*128 seedinstr

  ! REPORT DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     write(*,12) dnt_values(3),dnt_values(2),&
          dnt_values(1),dnt_values(5),dnt_values(6),&
          dnt_values(7),dnt_zone
     call system_clock(initial_time,count_rate)
  endif

  pi   = 3.14159
  ftpi = 4./3.*pi
  ! READ IN PARAMETERS
  open(unit=1,file="merge_params.txt")

  read(1,*) boxsize,ntile,ievol,cenx,ceny,cenz,maximum_redshift,&
       Omx,OmB,Omvac,h,iZeld,&
       iLexc,iLmrg,iFexc,iFmrg,ioutshear,iwrap,iwant_field_part,mlatt,largerun
  read(1,'(a)') filein
  read(1,'(a)') fileout
  close(1)

  seedin = i4arg(1,13579)

  write (seedinstr, *) seedin
  seedinstr = adjustl(seedinstr)

  filein = trim(filein)//'.'//trim(seedinstr)
  fileout = trim(fileout)//'.'//trim(seedinstr)

  numtiles      = ntile**3
  dcore_box = boxsize / ntile

  call allocate_harray(ioutshear)


  ! SETUP COSMOLOGY
  call Dlinear_cosmology(Omx,OmB,Omvac,h,iamcurved,dcurv)

  ! RANDOM TILE LIST
  iseed=13580
  allocate(tilelist(numtiles))
  do i=1,numtiles
    tilelist(i)=i
  enddo
  do i=numtiles,2,-1
    xcur = ran(iseed)
    j = int(xcur * i) + 1
    if(j<1.or.j>i) then
      write(*,*) 'j out of bounds',x,j
      call mpi_finalize(ierr)
      stop
    endif
    m = tilelist(j)
    tilelist(j) = tilelist(i)
    tilelist(i) = m    
  enddo

  ! OPEN OUTPUT FILE 
  open(1,file=fileout,status='unknown',access='stream')
  ! -----------------------------------------------------------------------
  ! LOOP OVER TILES
  ! --------------------------------------------------------------------
  neach = ceiling((numtiles+0.0)/ntasks)
  noutl = 0
  Rmaxl = -1e10
  Rminl = 1e10
  Volhl = 0.

  do tile_index=myid+1,ntasks*neach,ntasks
  
  if(tile_index>numtiles) goto 201
  tile = tilelist(tile_index)

  ! index of tile
  ktile = (tile - 1 )                     / ntile**2 + 1
  jtile = (tile - ntile**2*(ktile-1) - 1) / ntile    + 1
  itile =  tile - ntile**2*(ktile-1) - (jtile-1)*ntile     

  ! center of tile in Mpc
  xtile = (itile - 0.5) * dcore_box - boxsize/2 + cenx
  ytile = (jtile - 0.5) * dcore_box - boxsize/2 + ceny
  ztile = (ktile - 0.5) * dcore_box - boxsize/2 + cenz

  ! distance in Mpc to nearest corner of tile
  xmintile = abs(xtile) - dcore_box/2 
  ymintile = abs(ytile) - dcore_box/2
  zmintile = abs(ztile) - dcore_box/2

  chitile = sqrt(xmintile**2+ymintile**2+zmintile**2)

  redshifttile = 1 / afn(chitile) - 1
  ! skip tile if beyond maximum redshift
  if(redshifttile > maximum_redshift.and.ievol==1) cycle

  !READ IN HALOS
  call read_halos

  !PERFORM EXCLUSION
  pairwise_code = 'exclusion'
  call pairwise(pairwise_code,x,y,z,r,vx1,vy1,vz1,vx2,vy2,vz2,harray,nharray,nin,nexclude,xtile,ytile,ztile,dL_box) 
  nmerge = nexclude

  pairwise_code = 'reduction'
  call pairwise(pairwise_code,x,y,z,r,vx1,vy1,vz1,vx2,vy2,vz2,harray,nharray,nexclude,nmerge,xtile,ytile,ztile,dL_box)

  xL = x
  yL = y
  zL = z
  do i=1,nmerge 
     ! Move halos to final Eulerian position
     x(i)   = x(i) + vx1(i) + vx2(i)
     y(i)   = y(i) + vy1(i) + vy2(i)
     z(i)   = z(i) + vz1(i) + vz2(i)

     if(ievol==1) then ! lightcone
        chipk = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
        apk   = afn(chipk)
        Zpk   = 1.0/apk
     else                    ! static box
        Zpk   = boxredshift+1
        apk   = 1./Zpk
     endif
     Dpk  = Dlinear(apk,chipk,HD_Ha,D_a)
     hubb = 100.0*hub_Dlin(apk)
     HDa   = hubb*HD_Ha*apk

     ! Change from Displacements to velocity
     vx1(i) = HDa * ( vx1(i) + 2*vx2(i) )
     vy1(i) = HDa * ( vy1(i) + 2*vy2(i) )
     vz1(i) = HDa * ( vz1(i) + 2*vz2(i) )
  enddo

  do i=1,nmerge  
     ! Keep only halos in original Lagrangian volume to save double counting
     fxpk = abs(2*(xL(i)-xtile)/dcore_box)
     if(fxpk<=1) then
        fypk = abs(2*(yL(i)-ytile)/dcore_box)
        if(fypk<=1) then
           fzpk = abs(2*(zL(i)-ztile)/dcore_box)
           if(fzpk<=1) then
              Volhl = Volhl + ftpi*r(i)**3

              noutl=noutl+1
              xo(noutl)       = x(i)
              yo(noutl)       = y(i)
              zo(noutl)       = z(i)
              vxo(noutl)      = vx1(i)
              vyo(noutl)      = vy1(i)
              vzo(noutl)      = vz1(i)
              ro(noutl)       = r(i)
              xLo(noutl)      = xL(i)
              yLo(noutl)      = yL(i)
              zLo(noutl)      = zL(i)
              harrayo(:,noutl) = harray(:,i) 
           endif
        endif
     endif
  enddo

  enddo
  ! -----------------------------------------------------------------------
  ! END LOOP OVER TILES 
  ! -----------------------------------------------------------------------

  201 continue

  call mpi_barrier(mpi_comm_world,ierr)

  allocate(Npk_eachtaskl(0:ntasks-1))
  allocate(Npk_eachtask (0:ntasks-1))
  allocate(Npk_begin    (0:ntasks-1))

  Npk_begin           = 0
  Npk_eachtaskl       = 0
  Npk_eachtaskl(myid) = noutl

  call mpi_allreduce(Npk_eachtaskl,Npk_eachtask,ntasks,mpi_integer,mpi_sum,&
                      mpi_comm_world,ierr)

  Npk_begin(0) = 0
  do i=1,ntasks-1
     Npk_begin(i)=Npk_eachtask(i-1)+Npk_begin(i-1)
  enddo

  Rminl = min(Rminl,minval(ro(:noutl)))
  Rmaxl = max(Rmaxl,maxval(ro(:noutl)))

  if(ioutshear==0) then
     outnum = 11 ! position + vel + radius + lag_pos + F                                                                
  elseif(ioutshear>=1) then
     outnum = 33 ! position + vel + radius + lag_pos + F + e + p + strain + d2F + zform + grad(1:3) + gradf(1:3) + Rf + FRf + d2FRf  + gradRf(1:3)                                                                                                   
  endif


  ! WRITE PEAKS 
  pos_offset = int(4,8)*(3+outnum*Npk_begin(myid))+1
  write(1,pos=pos_offset) (xo(i),yo(i),zo(i),vxo(i),vyo(i),vzo(i),&
       ro(i),xLo(i),yLo(i),zLo(i),harrayo(:,i),i=1,noutl)
  nout = Npk_begin(ntasks-1) + Npk_eachtask(ntasks-1) 


  deallocate(Npk_eachtaskl)
  deallocate(Npk_eachtask )
  deallocate(Npk_begin    )

  call mpi_reduce(Rminl,Rmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(Rmaxl,Rmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  call mpi_reduce(Volhl,Volh,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
  fHalo = Volh/boxsize**3

  if(iwant_field_part==1) then
     if(myid==0) then
        open(unit=4,file="field_particles.bin",access='stream')
        pos_offset = 1+8
        read(4,pos=pos_offset)fField
        
        Mfield_p =  mlatt * (1-fHalo)/fField

        write(4,pos=pos_offset) Mfield_p
        close(4)
     endif
  endif
  ! WRITE HEADER IF myid == 0       
  if(myid==0) write(1,pos=1) nout, Rmax, boxredshift

  ! CLOSE OUTPUT FILE
  close(1)

  ! PRINT FINAL STATS
  if(myid==0) then
     write(*,101) ntile,ntasks
     write(*,102) nout
     write(*,103) Rmin,Rmax
  endif
     
  ! DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     call system_clock(final_time,count_rate)

     elapsed_time    = (final_time - initial_time) / count_rate
     elapsed_hours   = int(elapsed_time / 3600)
     elapsed_minutes = int((elapsed_time - elapsed_hours * 3600)/60)
     elapsed_seconds = int(elapsed_time - elapsed_hours * 3600 &
                                        - elapsed_minutes * 60)

     write(*,401) dnt_values(3),dnt_values(2),dnt_values(1),&
       dnt_values(5),dnt_values(6),dnt_values(7),&
       dnt_zone,elapsed_hours,elapsed_minutes,elapsed_seconds

  endif

  return

101 format(3x,'ntile, ntasks:   ',2(i0,1x))
102 format(3x,'number of peaks: ',i0)
103 format(3x,'Rmin = ', f7.3,' Rmax = ',f7.3)
 12 format(/,3x,61('-'),/,3x,'Peak-patch MERGEPKVD running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)           
401 format(/,3x,61('-'),/,3x,&
         'Peak-patch MERGEPKVD finished running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         40x,'Total run time: ',i2.2,':',i2.2,':',i2.2,/&
         3x,61('-'),/)

end subroutine merge_pkvd


subroutine read_halos

  use mpivars
  use arrays_params

  ! Read in Nhalo from individual raw tiles of hpkvd if largerun==1
  ! Else read in halos from full raw catalogue
  
  real harraypk(nharray)
  integer err

  ! READ ONLY HALOS FROM NEIGHBOURING TILES IF LARGERUN=1
  nin = 0
  do kk=-largerun,largerun
     if( ((ktile+kk) < 1) .or. ((ktile+kk) > ntile) ) cycle
     do jj=-largerun,largerun
        if( ((jtile+jj) < 1) .or. ((jtile+jj) > ntile) ) cycle
        do ii=-largerun,largerun
           if( ((itile+ii) < 1) .or. ((itile+ii) > ntile) ) cycle

           tilej      = ((ktile-1)+kk)*ntile**2 + ((jtile-1)+jj)*ntile + (itile+ii)
           write(tilejout,'(I6.6)') tilej ! constructs 6 character entity name for box input

           if(largerun==0) fileini = filein
           if(largerun==1) fileini = trim(filein)//'_'//trim(tilejout)

           inquire(file=fileini, exist=file_exists)
           if(.not. file_exists) cycle ! if neighbour has no halos/raw file then skip

           ! OPEN RAW HALO CATALOGUE
           open(4,file=fileini,status='unknown',access='stream')
           read(4) Npk, Rmax_in, boxredshift
           Rbuff  = 2 * Rmax_in
           dL_box = dcore_box + 2 * Rbuff


           ! DOMAIN DECOMPOSITION FILTERING OF INPUT DATA
           do i=1,Npk

              if(nin+1 > nmax) then
                 if(myid==0) write(*,*) 'Tiles too large, exiting... nin, nmax = ',nin,nmax
                 call mpi_finalize(err)
                 stop
              endif

              read(4) xpk,ypk,zpk,vx1pk,vy1pk,vz1pk,rpk,vx2pk,vy2pk,vz2pk, (harraypk(pk),pk=1,nharray)

              if(rpk <= 0.) cycle
 
              fxpk = abs(2*(xpk-xtile)/dL_box)
              if(fxpk<=1) then
                 fypk = abs(2*(ypk-ytile)/dL_box)
                 if(fypk<=1) then
                    fzpk = abs(2*(zpk-ztile)/dL_box)
                    if(fzpk<=1) then                       
                       nin         = nin + 1
                       x(nin)      = xpk
                       y(nin)      = ypk
                       z(nin)      = zpk
                       vx1(nin)    = vx1pk
                       vy1(nin)    = vy1pk
                       vz1(nin)    = vz1pk              
                       r(nin)      = rpk 
                       vx2(nin)    = vx2pk
                       vy2(nin)    = vy2pk
                       vz2(nin)    = vz2pk              
                       harray(:,nin) = harraypk

                    endif
                 endif
              endif
              
              !PERIODIC WRAP HALOS
              if (iwrap==1) then
                 if(abs(xtile) > (boxsize - dL_box - cenx)/2) then  ! if xtile is on an edge of the large box
                    fxpk = abs(xpk) - ( boxsize/2 - 3./2*Rbuff - cenx) ! if xpk is within 2*Rbuff of the boxedge
                    if(fxpk>0) then
                       fypk = abs(2*(ypk-ytile)/dL_box)
                       if(fypk<=1) then                             ! if ypk is within y range of tile
                          fzpk = abs(2*(zpk-ztile)/dL_box)  
                          if(fzpk<=1) then                          ! if zpk is within z range of tile
!                             if(myid==0) write(*,*) xtile,ytile,ztile,xpk,ypk,zpk
                             nin        = nin + 1
                             x(nin)     = xpk - xpk/abs(xpk)*boxsize !gets the sign correct                            
                             y(nin)     = ypk
                             z(nin)     = zpk
                             vx1(nin)   = vx1pk
                             vy1(nin)   = vy1pk
                             vz1(nin)   = vz1pk
                             r(nin)     = rpk
                             vx2(nin)   = vx2pk
                             vy2(nin)   = vy2pk
                             vz2(nin)   = vz2pk
                             harray(:,nin) = harraypk
                          endif
                       endif
                    endif
                 endif
                 if(abs(ytile) > (boxsize - dL_box - ceny)/2) then  ! if ytile is on an edge of the large box
                    fypk = abs(ypk) - ( boxsize/2 - 3./2*Rbuff - ceny) ! if ypk is within 2*Rbuff of the boxedge
                    if(fypk>0) then
                       fxpk = abs(2*(xpk-xtile)/dL_box)
                       if(fxpk<=1) then                             ! if xpk is within x range of tile
                          fzpk = abs(2*(zpk-ztile)/dL_box)  
                          if(fzpk<=1) then                          ! if zpk is within z range of tile
                             nin        = nin + 1
                             x(nin)     = xpk 
                             y(nin)     = ypk - ypk/abs(ypk)*boxsize !gets the sign correct                            
                             z(nin)     = zpk
                             vx1(nin)   = vx1pk
                             vy1(nin)   = vy1pk
                             vz1(nin)   = vz1pk
                             r(nin)     = rpk
                             vx2(nin)   = vx2pk
                             vy2(nin)   = vy2pk
                             vz2(nin)   = vz2pk
                             harray(:,nin) = harraypk
                          endif
                       endif
                    endif
                 endif
                 if(abs(ztile) > (boxsize - dL_box - cenz)/2) then  ! if ztile is on an edge of the large box
                    fzpk = abs(zpk) - ( boxsize/2 - 3./2*Rbuff - cenz) ! if zpk is within 2*Rbuff of the boxedge
                    if(fzpk>0) then
                       fypk = abs(2*(ypk-ytile)/dL_box)
                       if(fypk<=1) then                             ! if ypk is within y range of tile
                          fxpk = abs(2*(xpk-xtile)/dL_box)  
                          if(fxpk<=1) then                          ! if xpk is within z range of tile               
                             nin        = nin + 1
                             x(nin)     = xpk
                             y(nin)     = ypk
                             z(nin)     = zpk - zpk/abs(zpk)*boxsize !gets the sign correct                            
                             vx1(nin)   = vx1pk
                             vy1(nin)   = vy1pk
                             vz1(nin)   = vz1pk
                             r(nin)     = rpk
                             vx2(nin)   = vx2pk
                             vy2(nin)   = vy2pk
                             vz2(nin)   = vz2pk
                             harray(:,nin) = harraypk
                          endif
                       endif
                    endif
                 endif

              endif

           enddo
           
        enddo
     enddo
  enddo

  close(4)
  if(myid==0) write(*,*) nin,'halos read...'

end subroutine read_halos
