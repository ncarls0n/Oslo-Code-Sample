module arrays_params

  parameter (nmax=10000000)

  real x(nmax),y(nmax),z(nmax),xL(nmax),yL(nmax),zL(nmax)
  real vx1(nmax),vy1(nmax),vz1(nmax),vx2(nmax),vy2(nmax),vz2(nmax) ! 1LPT and 2LPT displacements           
  real r(nmax)
  real, allocatable :: harray(:,:), harrayo(:,:)

  real xo(nmax),yo(nmax),zo(nmax),xLo(nmax),yLo(nmax),zLo(nmax)
  real vxo(nmax),vyo(nmax),vzo(nmax)
  real ro(nmax)

  ! INPUT PARAMS
  real    Omx,OmB,Omvac,h,mlatt,boxredshift
  integer iZeld,iLexc,iLmrg,iFexc,iFmrg,iwrap,iwant_field_part,largerun
  integer iseed,iamcurved,dcurv,ioutshear

  ! June 2024 - New variables that need to be declared in order for gcc not to complain
  integer ievol
  real maximum_redshift

  ! DOMAIN DECOMPOSITION VARIABLES                                                        
  integer, allocatable :: tilelist(:)

  real    chipk, dcore_box, dL_box, cenx, ceny, cenz, xtile, ytile, ztile
  integer Npk, Non2, dummy, nin, nexclude, nmerge, noutl, nout, nharray
  real    boxsize, Rbuff, eff_read
  integer i,j,k,m,itile,jtile,ktile,ii,jj,kk,ntile,numtiles,tile,tilej,neach,tile_index,outnum

  integer, allocatable :: Npk_eachtask(:), Npk_eachtaskl(:), Npk_begin(:)

  ! I/O
  character *512 filein, fileini, fileout, pairwise_code, tilejout
  integer(kind=8) pos_offset
  logical  file_exists

contains 
  subroutine allocate_harray(ioutshear)

    if(ioutshear==0) then
       nharray = 1
    elseif(ioutshear>=1) then
       nharray = 23
    endif

    allocate(harray(nharray, nmax))
    allocate(harrayo(nharray, nmax))

  end subroutine allocate_harray

end module arrays_params
