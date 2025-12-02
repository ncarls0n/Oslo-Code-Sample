program padICs

  use textlib

  implicit none

  integer*4            nin,nout,i,j,k,kk,nbuff,pps,inflag
  integer*4            na1,na2,nb1,nb2,nc1,nc2,ma1,ma2,mc1,mc2
  character *512       filein, fileout, cnin, cnout, topf, midf, botf
  real, allocatable :: buffer(:,:,:)
  real                 ram

  call get_command_argument(1,filein)
  call get_command_argument(2,fileout)
  nin    = i4arg(3,4096)
  nbuff  = i4arg(4,32)
  ram    = r4arg(5,3.5) ! available RAM in GB
  inflag = i4arg(6,1)

  if(command_argument_count().lt.4) then
     write(*,99)
99      format('usage: padICs <filein> <fileout> <nin> <nbuff> <RAM (GB)> <inflag>')
     stop
  endif

  topf=trim(fileout)//'.top'
  midf=trim(fileout)//'.mid'
  botf=trim(fileout)//'.bot'

  open(1,file=  filein,access='stream')
  open(2,file=    topf,access='stream')
  open(3,file=    midf,access='stream')
  open(4,file=    botf,access='stream')
 
  nout = nin + 2 * nbuff

  na1 = 1
  na2 = nbuff
  nb1 = na2+1
  nb2 = nb1+nin-1
  nc1 = nb2+1
  nc2 = nout

  ma1 = nb2 - nbuff + 1
  ma2 = nb2

  mc1 = nbuff+1
  mc2 = 2*nbuff

  pps = int(1024.**3*ram/nout**2/4) - 1 ! pps = planes per slab
  pps = min(pps,nout)
  if(pps<1) then
     write(*,*) 'not enough memory, one plane is ',nout**2*4/1024.**3,'GB'
     write(*,*) 'limit: ',ram,'GB'
     stop
  endif

  allocate(buffer(nout,nout,pps))

  if(inflag==0) goto 11

  k=0 ! p = planes read since last slab
  do kk=1,nin
     write(*,*) kk,nin

     ! read a single plane and add buffer
     k = k + 1
     do j=nb1,nb2
        read(1) buffer(nb1:nb2,j,k)
        buffer(na1:na2,j,k)=buffer(ma1:ma2,j,k)
        buffer(nc1:nc2,j,k)=buffer(mc1:mc2,j,k)
     enddo
     buffer(:,na1:na2,k)=buffer(:,ma1:ma2,k)
     buffer(:,nc1:nc2,k)=buffer(:,mc1:mc2,k)
     
     if(kk> nin-nbuff) write(2) buffer(:,:,k) ! write this slab out to top
     if(kk<=    nbuff) write(4) buffer(:,:,k) ! write this slab out to bottom

     if(k==pps.or.kk==nin) then ! memory used up, write data to middle        
        write(3) buffer(:,:,1:k)
        k=0
     endif
        
  enddo

  close(1)
  close(2)
  close(3)

11 continue

  ! now read in top, middle and bottom and write to output file
  open(2,file=fileout,access='stream')

  ! TOP
  open(1,file=topf,access='stream')
  k=0 ! p = planes read since last slab
  do kk=1,nbuff

     ! read a single plane 
     k = k + 1    
     read(1) buffer(:,:,k)
     if(k==pps.or.kk==nbuff) then ! memory used up, write data to output file
        write(2) buffer(:,:,1:k)
        k=0
     endif

     write(*,*) kk,nbuff,'top'
        
  enddo
  close(1)

  ! MIDDLE
  open(1,file=midf,access='stream')
  k=0 ! p = planes read since last slab
  do kk=1,nin

     ! read a single plane 
     k = k + 1
     read(1) buffer(:,:,k)
     if(k==pps.or.kk==nin) then ! memory used up, write data to output file
        write(2) buffer(:,:,1:k)
        k=0
     endif

     write(*,*) kk,nin,'middle'
        
  enddo
  close(1)

  ! BOTTOM
  open(1,file=botf,access='stream')
  k=0 ! p = planes read since last slab
  do kk=1,nbuff

     ! read a single plane 
     k = k + 1
     read(1) buffer(:,:,k)
     if(k==pps.or.kk==nbuff) then ! memory used up, write data to output file
        write(2) buffer(:,:,1:k)
        k=0
     endif

     write(*,*) kk,nbuff,'bottom'
        
  enddo
  
  close(1)
  close(2)

  stop

end program padICs
