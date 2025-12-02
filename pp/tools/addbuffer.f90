program addbuffer

  implicit none

  integer*4 nin,nout,i,j,nbuff
  integer*4 na1,na2,nb1,nb2,nc1,nc2,ma1,ma2,mc1,mc2
  character *128 filein, fileout, cnin, cnout
  real, allocatable :: buffer(:,:,:)

  call getarg(1,filein)
  call getarg(2,fileout)
  call getarg(3,cnin)
  call getarg(4,cnout)

  if(iargc().ne.4) then
     write(*,99)
99      format('usage: addbuffer <filein> <fileout> <nin> <nout>')
     stop
  endif

  read(cnin,*) nin
  read(cnout,*) nout
 
  nbuff = (nout - nin) / 2

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

  allocate(buffer(nout,nout,nout))
  open(1,file=filein,access='stream')
  open(2,file=fileout,access='stream')
  do i=nb1,nb2
     do j=nb1,nb2
        read(1) buffer(nb1:nb2,j,i)
        buffer(na1:na2,j,i)=buffer(ma1:ma2,j,i)
        buffer(nc1:nc2,j,i)=buffer(mc1:mc2,j,i)
     enddo
     buffer(:,na1:na2,i)=buffer(:,ma1:ma2,i)
     buffer(:,nc1:nc2,i)=buffer(:,mc1:mc2,i)
  enddo
  buffer(:,:,na1:na2)=buffer(:,:,ma1:ma2)
  buffer(:,:,nc1:nc2)=buffer(:,:,mc1:mc2)

  write(2) buffer

  close(1)
  close(2)

  stop

end program addbuffer
