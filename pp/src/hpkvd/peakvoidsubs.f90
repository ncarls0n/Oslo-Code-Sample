real function evolve_ellipse_function(Frho_in,e_in,poe_in) 
  
  use HomogeneousEllipsoid
  use input_parameters
  use params

  real, intent(in) :: Frho_in, e_in, poe_in
  real :: fcvir1p,Dvir1p,zvir1p
  integer :: iflagel,iuseinterp

  Frho = 10**(Frho_in)
  e_v   = e_in 
  p_v   = poe_in * e_in 

  zinit = zinit_fac_ell * Frho ! 20*Frho, zinit_fac_ell asssigned 20 in hpkvd.f90

  call evolve_ellipse_full(1,zvir1p,Dvir1p,fcvir1p)

  evolve_ellipse_function = zvir1p

  return

end function evolve_ellipse_function

real function fsc_of_z(z)

  use TabInterp

  implicit none

  real, intent(in) :: z
  real fv,dfv,zvir

  fv=TabInterpX2
  dfv=1e-4
  zvir=100

  ! A scale search, walking in toward z=0 to find redshift for 
  ! which peak collapses based on interpolation from collapse tables
  do while(zvir>1+z)
     ! Determine redshift at which given peak with strain given by trace(strain)=fv,
     ! ellipticity=0., prolateness=0.
     zvir=TabInterpInterpolate(fv,0.,0.) ! interpolates the value of TabInterpArray(fv,0.,0.), from src/modules/TabInterp/TabInterp.f90
     fv = fv - dfv ! next step for trace smaller, i.e. more collapsed
  enddo

  ! Now we have estimate of maximum fv=trace(strain) for which z>0 (i.e. a correlate to
  ! the average density of the largest peak to by the "present time").
  fsc_of_z = 10**fv ! return 10^fv

end function fsc_of_z
  
subroutine get_pks(z,Non,xbx,ybx,zbx,alattv,ired,Rsmooth)
  use input_parameters
  use mpivars
  use arrays

  USE intreal_types
  use params
  
  use HomogeneousEllipsoid
  use TabInterp

  implicit none

  real,    intent(in)    :: z, xbx, ybx, zbx, alattv(3), Rsmooth
  integer, intent(in)    :: ired
  integer, intent(inout) :: Non

  real :: cen1,cen2,cen3
  real :: xc,yc,zc,rc,ac,ff,rdsc,fac,fv,zvir,dfv
  real :: Dvir1p,fcvir1p,zvir1p
  integer :: i,j,k,ii,jj,kk,iii,jjj,kkk

  real :: afn, fsc_of_z

! Non is the old value coming in, new one comes out

  cen1=0.5*(n1+1) ! (cen1,cen2,cen3) defines the coordinate of the centre of 
  cen2=0.5*(n2+1) ! the parallelization tile, not the centre of the entire
  cen3=0.5*(n3+1) ! simulation volume (unless you're running it as one tile)

  fv = fsc_of_z(z)

  do k=nbuff+1,n3-nbuff ! cycle over lattice cells not in buffer
     zc=zbx+alattv(3)*(k-cen3) ! z coordinate of tile excluding buffer
     do j=nbuff+1,n2-nbuff ! cycle over lattice cells not in buffer
        yc=ybx+alattv(2)*(j-cen2) ! y coordinate of tile excluding buffer
        do i=nbuff+1,n1-nbuff ! cycle over lattice cells not in buffer
           xc=xbx+alattv(1)*(i-cen1) ! x coordinate of tile excluding buffer

           if(mask(i,j,k)==1) cycle ! skips over any lattice cells that are beyond the horizon in light-cone runs
           ff=delta(i,j,k) ! ff = overdensity of current lattice cell

!           if(ievol==1) then
!              rc = sqrt(xc**2+yc**2+zc**2)
!              ac = afn(rc)
!              rdsc = 1 / ac - 1              
!              if(rdsc > maximum_redshift) cycle
!              fv = fsc_of_z(rdsc)
!           endif

           if(ff.lt.fv) cycle ! skips over lattice cells that haven't reached sufficient overdensity to be virialized

           do kk=-1,1     ! 
             kkk=k+kk     ! Nearest neighbours search to
             do jj=-1,1   ! determine whether delta(i,j,k)
               jjj=j+jj   ! is a local overdensity maximum.
               do ii=-1,1 ! 
                 iii=i+ii ! 

                 if(ff<delta(iii,jjj,kkk)) goto 99 ! skips over any lattice cells that aren't delta peaks

               enddo
             enddo
           enddo
    
           Non                = Non+1 ! count of peaks found +1
           xpk(Non,ired)      = xc ! 
           ypk(Non,ired)      = yc ! coordinates of newly located peak
           zpk(Non,ired)      = zc ! 
           if (ioutshear>0) FcollvRf(Non,ired) = ff
           vzpk(Non,ired)     = Rsmooth 
           vypk(Non,ired)     = etay(i,j,k)*fac
           vxpk(Non,ired)     = etax(i,j,k)*fac
           lagrange(Non,ired) = i+(j-1)*n1+(k-1)*n1*n2 ! converts an element (i,j,k) in the n1 x n2 x n3 lattice to an element of a single array of length n1*n2*n3
           mask(i,j,k)        = 1

           if(ioutshear>=1) F_d2Rf(Non,ired)   = lapd(i,j,k)

           if (Non.ge.Npkmaxl) return ! ends subroutine if Non above max allowed value Npkmaxl=200000
                                      ! which is assigned in src/hpkvd/arrays.f90
99         continue
        enddo
     enddo
  enddo

 return

end subroutine get_pks

subroutine get_homel(npart,ipp,alatt,ir2min, &
    Sbar,RTHL,Srb,ZZon,Fnupk,Fevpk,Fpvpk,Fd2pk,strain_mat,Sbar2,Rfclvi,gradpkrf,gradpk,gradpkf,zvir_half)
  USE intreal_types
  use memory_management
  ! use bound
  use mpivars
  use HomogeneousEllipsoid

  !c  Radial integration of F profile to find Fbar=f_v crossing point
  !c  closest to R^2=ir2min*alatt fiducial. Return RTHL and Srb.

  !c  Use homogeneous ellipsoid model to truncate integration

  !c  Get back shear eigenvalues F_nu,F_ev,F_pv where
  !c              lam1 = F_nu/3( 1 - 3*F_ev + F_pv )
  !c              lam2 = F_nu/3( 1 - 2*F_pv )
  !c              lam3 = F_nu/3( 1 + 3*F_ev + F_pv )

  !c  So trace     F_nu = lam1 + lam2 + lam3 = F
  !c  ellipticity  F_ev = ( lam3 - lam1 )/ 2*F_nu
  !c  prolateness  F_pv = ( lam1 + lam3 - 2*lam2 )/ 2*F_nu

  !c  These are calculated over a top-hat filtering
  !c  IF zvir=-1 IT DID NOT COLLAPSE ALONG AXIS 1 (last axis)

  use input_parameters   ! from src/modules/GlobalVariables, declares input parameters
  use arrays             ! from src/hpkvd, subroutines: allocate_halos, allocate_boxes
  use timing_diagnostics ! from src/modules/External, subroutines: timer_begin, timer_end

  use TabInterp          ! from src/modules/TabInterp, tables for interpolation

  implicit none
  integer(i4b) iuseinterp
  real(sp) fcvir1p,Dvir1p

  integer(i4b) nn(3),npart,ipp,naverage
  real(sp) Fnupk,Fevpk,Fpvpk,Fd2pk
  real(sp) zvir_full,zvir_half
  real(sp) etavec(3)
  real(sp) eta2vec(3)
  integer ipk(3),iv(3)
  integer, allocatable :: nshell(:)
  real,    allocatable :: Sshell(:,:),Gshell(:,:),Gshellf(:,:),SRshell(:,:,:),Fbar(:),rad(:),S2shell(:,:)
  real(sp) Sbar(3),Ebar(3,3),strain_mat(3,3)
  real(sp) Ebar_mean_full(3,3),Ebar_mean_half(3,3)
  real(sp) gradpk(3),gradpkrf(3),gradpkf(3),gradpk_p(3),gradpkf_p(3)

  integer  cell(3),nmean_full, nmean_half
  real(sp) flocal
  real(sp) Sbar2(3)
  real(sp) s2local(3),slocal(3)
  real(sp) strlocal(6)
  real time1,time2,time3,time4,time5

  real rcur,jj

  integer(i4b) iblack
  integer(i4b) ir2min,n1xn2,ir2upp,nsh,m,m0,m1,ir2p,ifcrit,&
               iflag,iflagel
  integer(i4b) L,K,jp,ir2,j0,icon,i2c,nSbar
  integer(i4b) mupp,muppnew,mupp_p,mlow,mlownew,mstart,mp,mrf,mrfi
  real(sp) fsc_tabfn_of_ZZ
  real(sp) RTHL,Srb,ZZon,alatt,fcrit,fcritx,zvir1,zvir1p
  real(sp) Frhoh,Frhpk,hlatt,hlatt2,hlatt_1,hlatt_2,diff,con,aww
  real(sp) pi,fourpi,anor,aRnor,wnor,wRnor,one_third
  real(sp) rupp,Fshell,dFbar,dFbarp,rad3,rad3p,rlow,u0,u02,u12
  real(sp) wt,dZvir,RTHL3,RTHL5,drad3,Fbarx,rad5,rad5p
  real(sp) Frhoi,rmrf
  real(sp) Frho_mean_full, Frho_mean_half

  real(sp) Lam(3)

  !c Functions
  real(sp) hRinteg

  !G.S filterscale peak found at
  real(sp) Rfclvi

  real(sp) Frhoc, Evec(3,3), e_vc, p_vc, zvir1pc, zvir1e, zdynax, fcvir1, poe

  real afnofD,fsc_of_z

  ! G.S. 10/03/2017 change sign of Ebar for corrected displacements (x=q+psi instead of x=q-psi)
  ! This change was to compensate for the fixing of the fourier convention in IC generation,
  ! Where now have the proper e^(ikx) for inverse FFT, NOT e^(-ikx) as was being used
  ! so
  !  Ebar(L,K) = Ebar(L,K)+&
  !                  0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
  ! becomes
  !        
  !  Ebar(L,K) = Ebar(L,K)-&         
  !                  0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))   
  !
  ! everywhere in the code. See BM1 eq 2.16, where e_{b,ij} = -1/2 *( ds_{bi}/dr_j + ds_{bj}/dr_i)

  allocate(nshell(npart+1))
  allocate(Sshell(3,npart+1))
  allocate(S2shell(3,npart+1))
  allocate(SRshell(3,3,npart+1))
  allocate(Fbar(npart+1))
  allocate(rad(npart+1))
  allocate(Gshell(3,npart+1))
  allocate(Gshellf(3,npart+1))
  
  if(report_memory) then ! Reports maximum memory usage
     do m0=1,npart
        nshell(m0)=m0*ipp      ! 
        Fbar(m0)=m0*ipp**2     ! ipp (jp in hpkvd) is the lattice site of the current
        rad(m0)=m0*ipp**3      ! patch determined by lagrange(jpp,ired)
        Sshell(1,m0)=m0*ipp**4 ! 
        Sshell(2,m0)=m0*ipp**5 ! 
        Sshell(3,m0)=m0*ipp**6 ! 
        do m=1,3
           SRshell(1,m,m0)=m0*ipp**7+m*npart**2
           SRshell(2,m,m0)=m0*ipp**8+m*npart**3
           SRshell(3,m,m0)=m0*ipp**9+m*npart**4
        enddo
     enddo
     call ReportMemory(C_CHAR_"High Water Mark"//C_NULL_CHAR)
     nshell=0
     Sshell=0
     SRshell=0
     Fbar=0
     rad=0
  endif

  RTHL=0.0 ! Top-hat filter Radius (in Lagrangian space) at which peak virializes
  Srb=1.0 ! 
  hlatt=1.0 ! just 1 as far as I can tell...

  ! fcrit=1.686, derived in Gunn & Gott doi:10.1086/151605, set in <run_directory>/tables/filter.dat
  fcrit = fsc_of_z(ZZon-1) ! critical overdensity for a peak at this redshift
  ! fcrit=fsc_tabfn_of_ZZ(ZZon)

  !c In-box correction factor D_on
  pi        = 4.0*atan(1.0)
  fourpi    = 4.0*pi
  one_third = 1.0/3.0

  ! Note that alatt is the Lagrangian real-space sidelength of a lattice cell (in Mpc/h)
  anor   = 1.0/hlatt**3   ! = 1./1.**3        = 1.
  aRnor  = 3.0*anor/alatt ! = 3.*1./alatt     = 3./alatt
  wnor   = anor/fourpi    ! = 1./(4*pi)       = 7.95774683E-02         = 1/4pi
  wRnor  = 3.0*wnor/alatt ! = 3./(4*pi*alatt) = (7.95774683E-02)/alatt = 1/(4pi alatt)

  nn(1) = n1 ! 
  nn(2) = n2 ! usually n1=n2=n3=nmesh
  nn(3) = n3 ! 
  n1xn2 = nn(1)*nn(2) ! nmesh^2
  call get_ijk(n1xn2,nn(1),ipp,ipk(1),ipk(2),ipk(3)) ! defined below
  ! ipk(1),ipk(2),ipk(3) = i,j,k define the position of the current peak on the lattice

  !G.S 22/09/2015
  !if (masked by peak larger than max search radius for filterscale) then (quit get_homel)
  if (rmax2rs > 0.0) then ! note rmax2rs default value in parameter file is 0.0
     if ( mask(ipk(1),ipk(2),ipk(3)) > int((Rfclvi/alatt)*rmax2rs) ) then
        RTHL=-1.0
        Srb=0.0
        deallocate(nshell,Sshell,S2shell,Gshell,Gshellf,SRshell,Fbar,rad)
        return
     endif
  endif
  !if masked by peak larger than the minimum search radius but less than
  !max search radius set minimum search radius to masked value
!  if ( int((0.8*mask(ipk(1),ipk(2),ipk(3)))**2) > ir2min) then
!     ir2min = int((0.8*mask(ipk(1),ipk(2),ipk(3)))**2)
!  endif

  hlatt2  = hlatt*hlatt     ! = 1.*1.           = 1.
  hlatt_1 = 1.0/hlatt       ! = 1./1.           = 1.
  hlatt_2 = hlatt_1*hlatt_1 ! = (1./1.)*(1./1.) = 1.

  rupp   = sqrt(float(ir2min))+2.0*hlatt ! = 1.75*Rfclv(ic)+2.0
  ir2upp = int(rupp*rupp)+1 ! rounds up the square of rupp to the nearest integer

  Fshell=0.0
  nsh=0
  !c  First point = first shell
  m  = 1     ! 
  m0 = 1     ! Shell m=1 is the centre, so rad(1)=0.
  rad(1)=0.0 ! 
  call checkiv(ipk,nn,0) ! checks that 1 <= ipk(j) <= nn(j) for j={1,2,3}, defined below
  Fbar(1)   = delta(ipk(1),ipk(2),ipk(3)) ! overdensity of initial field at peak site
  nshell(1) = 1
  ir2p      = irs2(2) ! irs2(2) is square of distance between adjacent lattice cells (set in icloud)
  dFbarp    = Fbar(1) ! Previous value of Fbar, Fbarp, set to current value, Fbar(1)
  etavec(1) = etax(ipk(1),ipk(2),ipk(3)) ! etax,etay,etaz together define the 1LPT displacement vector field
  etavec(2) = etay(ipk(1),ipk(2),ipk(3)) ! for the entire simulaiton volume, so these three lines define
  etavec(3) = etaz(ipk(1),ipk(2),ipk(3)) ! etavec(1:3) as the 1LPT displacement vector at the current peak.
  if(ilpt==2) then
     eta2vec(1) = eta2x(ipk(1),ipk(2),ipk(3)) ! 
     eta2vec(2) = eta2y(ipk(1),ipk(2),ipk(3)) ! 2LPT displacement vector field at current peak
     eta2vec(3) = eta2z(ipk(1),ipk(2),ipk(3)) ! 
  else
     eta2vec(1) = 0.0 ! 
     eta2vec(2) = 0.0 ! if only 1LPT is to be used, sets 2LPT peak displacement vector to 0
     eta2vec(3) = 0.0 ! 
  endif
  do L=1,3
     Sshell(L,1)  = etavec(L)  ! Sshell(1:3,1) = 1LPT peak displacement
     S2shell(L,1) = eta2vec(L) ! S2shell(1:3,1) = 2LPT peak displacement
     Sshell(L,2)  = 0.0 ! s_1LPT^L(q) summed over q in the shell 2
     S2shell(L,2) = 0.0 ! s_2LPT^L(q) summed over q in the shell 2
     Gshell(L,1)  = 0.0 ! zero for shell 1 (see next)
     Gshell(L,2)  = 0.0 ! delta(q)*(q^L-q_pk^L) summed over q in the shell 2
     Gshellf(L,1) = 0.0 ! zero for shell 1 (see next)
     Gshellf(L,2) = 0.0
     do K=1,3
        SRshell(L,K,1) = 0.0 ! zero for shell 1 (see next)
        SRshell(L,K,2) = 0.0 ! s_1LPT^L(q) * (q^K-q_pk^K) summed over all q^K in the shell 2
     enddo
  enddo
  
  !c  Fill array out to ir2min (need to go out to ir2max for Ebar)
  ifcrit=1 
  do jp=2,npart ! cycle over lattice cells within ~ 7/4 the max filter radius (measured in lattice units)
     ir2=irs2(jp) ! radial distance squared of cell jp from cell 1
     if(ir2.ne.ir2p) then ! if (cell not at same radial distance as previous) then, i.e. creates radial shells
        m=m+1                                ! starting from m=1, counts upward
        rad(m)    = sqrt(float(ir2p))        ! floating point distance to cell jp from cell 1
        nshell(m) = nsh                      ! # lattice sites in shell, set below in if(iblack.eq.0)
        dFbar=0.0                            ! dFbar=0. for shell m=1
        if(nsh.gt.0) dFbar=Fshell/float(nsh) ! dFbar for shells m>1, Fshell = total overdensity of shell
        rad3p   = rad(m-1)**3                ! previus value of radial distance cubed
        rad3    = rad(m)**3                  ! new value of radial distance cubed
        Fbar(m) = (rad3p*Fbar(m-1)+&         ! extrapolates overdensity of shell via truncated series expansion
             0.5*(dFbarp+dFbar)*(rad3-rad3p))/rad3 ! dFbarp initially set to Fbar(1), dFbar initially set to 0
        dFbarp  = dFbar                      ! sets dFbarp to previous value of dFbar

        if(ifcrit.eq.1) then ! ifcrit operates as boolian initially 1 and set to 0 if previous squared radius ir2p greater than a minimum value ir2min passed to the subroutine
           m0=m
           if(ir2p.gt.ir2min) then ! if radius squared in lattice units as integer > min( int((1.75*Rfclv(ic)/alatt)**2), int((40.0/alatt-1)**2)) or very roughly ~ int((7/4*min(Rfclvi,22 h^-1 Mpc)/alatt)^2)
              rupp   = rad(m0)+2.0*hlatt ! defines new max radius rupp from ~1.75 R_filter of current peak to radius of shell m0 + a little bit
              ir2upp = int(rupp*rupp)+1  ! rupp^2 as integer
              ifcrit = 0                 ! so that subsequent passes will skip the statement if(ifcrit.eq.1)
              mupp=m                     ! defines mupp as m of smallest shell with radius larger than rupp?
           endif
        else
           if(rad(m).lt.rupp) mupp=m ! if (current radius <~ 1.75 R_filter) sets mupp to m
        endif

        ir2p=ir2   ! resets previous radial coordinate ir2p to the current ir2
        Fshell=0.0
        nsh=0      
        do L=1,3
           Sshell(L,m+1)  = 0.0
           S2shell(L,m+1) = 0.0
           Gshell(L,m)    = Gshell(L,m)/rad(m)  ! initially 0
           Gshellf(L,m)   = Gshellf(L,m)/rad(m) ! initially 0
           Gshell(L,m+1)  = 0.0
           Gshellf(L,m+1) = 0.0
           do K=1,3
              SRshell(L,K,m)   = SRshell(L,K,m)/rad(m) ! initially 0
              SRshell(L,K,m+1) = 0.0                   ! initially 0
           enddo
        enddo

     endif
     iblack=0
     do L=1,3
        iv(L)=ipk(L)+ixsvec(jp,L) ! defines the position (in lattice units) of the jpth cell (where 0<jp<=npart) within a distance nhunt (in lattice units) of the peak in question ipk(1:3)
        if((iv(L).lt.1).or.(iv(L).gt.nn(L))) iblack=1 ! if iv(L) is outside the simulation tile
     enddo
     if(iblack.eq.0) then ! if iv(L) is in the simulaiton tile
        nsh=nsh+1
        call checkiv(iv,nn,1) ! checks that 1 <= iv(j) <= nn(j) for j={1,2,3}, defined below
        Fshell    = Fshell+delta(iv(1),iv(2),iv(3)) ! summing up total overdensity of shell
        etavec(1) = etax(iv(1),iv(2),iv(3)) ! 
        etavec(2) = etay(iv(1),iv(2),iv(3)) ! Value of the 1LPT displacement field at position iv(1:3)
        etavec(3) = etaz(iv(1),iv(2),iv(3)) !                                      (in lattice units)
        if(ilpt==2) then
           eta2vec(1) = eta2x(iv(1),iv(2),iv(3)) ! 
           eta2vec(2) = eta2y(iv(1),iv(2),iv(3)) ! Value of 2LPT displacement field at position iv
           eta2vec(3) = eta2z(iv(1),iv(2),iv(3)) !                             (in lattice units)
        endif
        do L=1,3
           Sshell(L,m+1)  = Sshell(L,m+1)  + etavec(L) ! vector sum of 1LPT displaceemnt field contributions from all lattice sites in the shell, s_1LPT^L(q) summed over q in the shell
           S2shell(L,m+1) = S2shell(L,m+1) + eta2vec(L) ! vector sum of 2LPT displacement field contributions from all lattice sites in the shell, s_2LPT^L(q) summed over q in the shell
           Gshell(L,m+1)  = Gshell(L,m+1)  + delta(iv(1),iv(2),iv(3))*ixsvec(jp,L) ! Gshell(L,m) defines the Lth component of the vector sum (over points on shell m) of overdensity times vector from peak centre to point on shell, i.e. delta(q)*(q^L-q_pk^L) summed over q in the shell m
           Gshellf(L,m+1) = Gshellf(L,m+1) + delta(iv(3),iv(1),iv(2))*ixsvec(jp,L) ! similar to Gshell, but for some reason looking at delta(jki) rather than delta(ijk)
           do K=1,3
              SRshell(L,K,m+1) = SRshell(L,K,m+1)+& ! Matrix object s_1LPT^L(q) * (q^K-q_pk^K)
                   etavec(L)*ixsvec(jp,K)           ! summed over all q^K in the shell
           enddo
        enddo
     endif
     !c  Now bounce out if at ir2upp
     if (ir2.gt.ir2upp) goto 252
  enddo ! Finish cycling over shells

252 continue



  ! GET GRADIENT AT Rf
  ! finding gradient at shell closest to Rf by considering shells of radius
  ! Rf-2<r<Rf+2 (in lattice units). Works only if ir2min > (rf/alatt)^2
  rmrf = 10.0
  mrfi = 1 ! shell number
  ir2p=irs2(2) ! = lattice spacing
  do jp=2,npart ! cycle over all sites within nhunt of q_pk 
     ir2=irs2(jp) ! radial distance squared to jp (in lattice units)
     if(abs(ir2-(Rfclvi/alatt)**2) < rmrf) then ! if ir2 - (Rf^2 in lattice units) < rmrf
        mrf = mrfi                                  ! finds shell number associated with radius rmrf
        rmrf = min(rmrf,abs(ir2-(Rfclvi/alatt)**2)) ! finds closest shell radius (squared) ir2 to Rf^2
     endif
     if(ir2.ne.ir2p) then ! as long as Rf > lattice spacing
        mrfi=mrfi+1 ! shell number
        ir2p=ir2    ! ir2_previous
     endif
  enddo
  ! Now we have:
  ! rmrf = closest shell radius squared to Rf^2
  ! mrf  = the shell number corresponding to mrf

  !c  Find mlow
  rlow=max(rad(mrf)-2.0*hlatt,0.0) ! max of (radius in lattice units of shell w/ radius nearest Rf)-2, or zero
  mlow=mrf                         ! shell number closest in radius to Rf
  if(mrf.gt.1) then ! if not the central cell
     do m1=mrf-1,1,-1
        if(rad(m1).gt.rlow) mlow=m1 ! sets mlow to shell of radius just larger than rlow
     enddo
  endif
  !c  Find mupp
  rupp=rad(mrf)+2.0*hlatt ! (radius in lattice units of shell w/ radius nearest Rf)+2
  mupp=mrf                ! shell number closest in radius to Rf
  do m1=mrf+1,mupp
     if(rad(m1).lt.rupp) mupp=m1 ! sets mupp to shell of radius just smaller than rupp
  enddo

  !c  Sum over particles on shells mlow to mupp
  gradpkrf = 0
  do m1=mlow,mupp
     u0  = hlatt_1*(rad(mrf)-rad(m1)) ! = rad(mrf)-rad(m1) since hlatt_1 is by default equal to 1.0
     u02 = u0*u0
     if(u02.lt.4.0) then
        u12 = hlatt_2*(rad(mrf)**2+rad(m1)**2) ! = rad(mrf)^2 + rad(m1)^2
        wt  = hRinteg(u0,u12)/(u12-u02)**2 ! function hRinteg defined below, integrates spherical kernel from |u0| to 2
        do L=1,3
           gradpkrf(L) = gradpkrf(L)+wt*Gshell(L,mrf) 
        enddo
     endif
  enddo
  gradpkrf = wRnor*gradpk/rad(mrf) ! Normalisation: gradpk/rad(mrf)/(4pi*alatt)

  !c  Check at ir2min (m0) for Fbar=fcrit
  if(Fbar(m0).ge.fcrit) then
     !c  Go outward to Fbar<fcrit
     ifcrit=1
     j0=jp+1
     do jp=j0,npart
        ir2=irs2(jp)
        if(ir2.ne.ir2p) then
           m=m+1
           rad(m)=sqrt(float(ir2p))
           nshell(m)=nsh
           dFbar=0.0
           if(nsh.gt.0) dFbar=Fshell/float(nsh)
           rad3p   = rad(m-1)**3
           rad3    = rad(m)**3
           Fbar(m) = (rad3p*Fbar(m-1)+&
                0.5*(dFbarp+dFbar)*(rad3-rad3p))/rad3
           dFbarp  = dFbar

           !c  check for Fbar<fcrit : make sure go out to r+2h
           if(ifcrit.eq.1) then
              m0=m
              rupp=rad(m0)+2.0*hlatt
              if(Fbar(m).lt.fcrit) ifcrit=0
           else
              if(rad(m).ge.rupp) then
                 mupp=m-1
                 goto 254
              endif
           endif

           ir2p=ir2
           Fshell=0.0
           nsh=0
           do L=1,3
              Sshell(L,m+1)  = 0.0
              S2shell(L,m+1) = 0.0
              Gshell(L,m)    = Gshell(L,m)/rad(m) ! Gshell is as above but divided by rad(m)
              Gshellf(L,m)   = Gshellf(L,m)/rad(m)
              Gshell(L,m+1)  = 0.0
              Gshellf(L,m+1) = 0.0
              do K=1,3
                 SRshell(L,K,m)   = SRshell(L,K,m)/rad(m)
                 SRshell(L,K,m+1) = 0.0
              enddo
           enddo
        endif
        iblack=0
        do L=1,3
           iv(L)=ipk(L)+ixsvec(jp,L)
           if((iv(L).lt.1).or.(iv(L).gt.nn(L))) iblack=1
        enddo
        if(iblack.eq.0) then
           nsh=nsh+1
           call checkiv(iv,nn,2) ! checks that 1 <= iv(j) <= nn(j) for j={1,2,3}, defined below
           Fshell    = Fshell+delta(iv(1),iv(2),iv(3))
           etavec(1) = etax(iv(1),iv(2),iv(3))
           etavec(2) = etay(iv(1),iv(2),iv(3))
           etavec(3) = etaz(iv(1),iv(2),iv(3))
           if(ilpt==2) then
              eta2vec(1) = eta2x(iv(1),iv(2),iv(3))
              eta2vec(2) = eta2y(iv(1),iv(2),iv(3))
              eta2vec(3) = eta2z(iv(1),iv(2),iv(3))
           endif
           do L=1,3
              Sshell(L,m+1)  = Sshell(L,m+1) + etavec(L)
              S2shell(L,m+1) = S2shell(L,m+1) + eta2vec(L)
              Gshell(L,m+1)  = Gshell(L,m)   + delta(iv(1),iv(2),iv(3))*ixsvec(jp,L)
              Gshellf(L,m+1) = Gshellf(L,m)  + delta(iv(3),iv(1),iv(2))*ixsvec(jp,L)
              do K=1,3
                 SRshell(L,K,m+1) = SRshell(L,K,m+1)+&
                      ixsvec(jp,K)*etavec(L)
              enddo
           enddo
        endif
     enddo

     !c  To reach here, end of array was encountered without fcrit
     
     if(ifcrit.eq.1) then
        mupp=m     
     else
        deallocate(nshell,Sshell,S2shell,Gshell,Gshellf,SRshell,Fbar,rad)
        return
     endif
254  continue

  else
     !c  Back down through shells to find fcrit

     mstart=max(m0-1,1)
     do mp=mstart,1,-1
        if(Fbar(mp).ge.fcrit) goto 256
        m0=mp
     enddo

     !c  To reach here, first shell was reached without downcross

     RTHL=-1.0
     Srb=0.0
     deallocate(nshell,Sshell,S2shell,Gshell,Gshellf,SRshell,Fbar,rad)
     return

256  continue

     !c  Find mupp
     rupp=rad(m0)+2.0*hlatt
     muppnew=m0
     do m1=m0+1,mupp
        if(rad(m1).lt.rupp) muppnew=m1
     enddo
     mupp=muppnew

  endif
  !c  Now have Fbar(m0-1) >= fcrit > Fbar(m0)
  !c  Check inward for homeoellipse virialization
  !c  Find mlow
  rlow=max(rad(m0)-2.0*hlatt,0.0)
  mlow=m0
  if(m0.gt.1) then
     do m1=m0-1,1,-1
        if(rad(m1).gt.rlow) mlow=m1
     enddo
  endif
  mupp_p=mupp

  !c  Zvir at m0
  !c  Need to sum over particles on shells mlow to mupp
  Ebar    = 0
  gradpk  = 0
  gradpkf = 0
  do m1=mlow,mupp
     u0  = hlatt_1*(rad(m0)-rad(m1))
     u02 = u0*u0
     if(u02.lt.4.0) then
        u12 = hlatt_2*(rad(m0)**2+rad(m1)**2)
        wt  = hRinteg(u0,u12)/(u12-u02)**2
        do L=1,3           
           gradpk(L)  = gradpk(L)+wt*Gshell(L,m1) 
           gradpkf(L) = gradpkf(L)+wt*Gshellf(L,m1) 
           do K=1,3
              Ebar(L,K) = Ebar(L,K)-&
                   0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
           enddo
        enddo
     endif
  enddo

  gradpk = wRnor*gradpk/rad(m0)
  gradpkf = wRnor*gradpkf/rad(m0)
  do L=1,3
     do K=1,3
        Ebar(L,K) = wRnor*Ebar(L,K)/rad(m0)
     enddo
  enddo
  !c  Solve for eigenvalues
  call get_evals(3,Ebar,Lam,iflag,Evec) ! computes eigenvalues Lam(1:3) of Ebar(3:3), defined below
  !c  These should be in INCREASING order. 
  Frho = Lam(1) + Lam(2) + Lam(3)
  if(Frho.gt.0.0) then
     e_v = 0.5*(Lam(3) - Lam(1))/Frho
     p_v = 0.5*(Lam(3) + Lam(1) - 2.0*Lam(2))/Frho
  else
     e_v = 0.0
     p_v = 0.0
  endif
  !c  Renormalize to Fbar
  Frhoh= Frho           ! 
  Frho = Fbar(m0)       ! Fbar(m0-1) >= fcrit > Fbar(m0), so now Frho<fcrit
  strain_mat = Ebar     ! 
  gradpk_p   = gradpk   ! 
  gradpkf_p   = gradpkf ! 
  zvir1p=-1
  if(iflag.eq.0.and.Frho>0) then ! iflag==0 if eigenvalues found by get_evals, overdensity should be >0
     e_vc = e_v
     p_vc = p_v
     poe = p_vc/e_vc
     if(e_vc < 1e-5) poe=0
     zvir1p = TabInterpInterpolate(log10(Frhoc),e_vc,poe) ! interpolates redshift of virialization given eigenvalues, from src/modules/TabInterp/TabInterp.f90, Frhoc not defeined yet, so this just returns zvir1p=-1 for every peak
     !     call gethom_ellipse(fcvir1p,Dvir1p,zvir1p,iflagel,iuseinterp)
  endif

  !c  Check to see that things haven't virialized yet
  if(zvir1p.ge.ZZon) then ! zvir1p should be -1 (because Frhoc=0. above), ZZon is 1 + collapse redshift for filter at which peak was found, so we expect ZZon == 1. for non-lightcone runs
     zvir1=-1.0
     goto 300
  endif
  zvir1=zvir1p ! = -1.
  Frhpk=Frhoh  ! Don't understand how this is different from Frhoh defined in the next loop (see line 766)
  Fnupk=Frho
  Fevpk=e_v
  Fpvpk=p_v

  if(m0.eq.1) goto 299 ! if the shell closest in radius to Rf is the zeroth shell
  !c  Step down through shells
  mstart=m0-1
  do mp=mstart,1,-1
     rupp=rad(mp)+2.0*hlatt
     rlow=max(rad(mp)-2.0*hlatt,0.0)
     if(rad(mupp).gt.rupp) then
        !c  New mupp
        muppnew=mp
        do m1=mp+1,mupp
           if(rad(m1).lt.rupp) muppnew=m1
        enddo
        mupp=muppnew
     endif

     if(mlow.gt.1) then
        if(rad(mlow-1).gt.rlow) then
           !c  New mlow
           mlownew=mlow
           do m1=mlow-1,1,-1
              if(rad(m1).gt.rlow) mlownew=m1
           enddo
           mlow=mlownew
        endif
     endif
     !c  Check ellipsoid
     Ebar    = 0
     gradpk  = 0
     gradpkf = 0
     if(mp.gt.1) then
        do m1=mlow,mupp
           u0  = hlatt_1*(rad(mp)-rad(m1))
           u02 = u0*u0
           if(u02.lt.4.0) then
              u12 = hlatt_2*(rad(mp)**2+rad(m1)**2)
              wt  = hRinteg(u0,u12)/(u12-u02)**2
              do L=1,3
                 gradpk(L)  = gradpk(L)+wt*Gshell(L,m1)
                 gradpkf(L) = gradpkf(L)+wt*Gshellf(L,m1)
                 do K=1,3
                    Ebar(L,K) = Ebar(L,K)-&
                         0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
                 enddo
              enddo
           endif
        enddo

        Ebar    = wRnor*Ebar/rad(mp)
        gradpk  = wRnor*gradpk/rad(mp)
        gradpkf = wRnor*gradpkf/rad(mp)
     else
        do m1=mlow,mupp
           u0=hlatt_1*rad(m1)
           con=100.0*u0*u0
           icon=int(con)
           diff=con-icon
           i2c=icon+1
           !c  linear interp
           aww=akk(i2c) + diff*(akk(i2c+1) - akk(i2c))
           wt=0.5*aww/rad(m1)
           do L=1,3
              gradpk(L)  = gradpk(L)+wt*Gshell(L,m1)
              gradpkf(L) = gradpkf(L)+wt*Gshellf(L,m1)
              do K=1,3
                 Ebar(L,K) = Ebar(L,K)-&
                      0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
              enddo
           enddo
        enddo

        Ebar    = aRnor*Ebar
        gradpk  = aRnor*gradpk
        gradpkf = aRnor*gradpkf

     endif

     call get_evals(3,Ebar,Lam,iflag,Evec) ! computes eigenvalues Lam(1:3) of Ebar(1:3,1:3), defined below
     Frho = Lam(1) + Lam(2) + Lam(3)
     if(Frho.gt.0.0) then
        e_v = 0.5*(Lam(3) - Lam(1))/Frho
        p_v = 0.5*(Lam(3) + Lam(1) - 2.0*Lam(2))/Frho
     else
        e_v = 0.0 ! doesn't bother with ellipticity or prolateness
        p_v = 0.0 ! for underdense regions
     endif

     Frhoh= Frho     ! Don't understand how this is different from Frhoh defined in previous loop (see line 675)
     Frho = Fbar(mp)

     zvir1p=-1
     if(iflag.eq.0.and.Frho>0) then
        Frhoc = Frho
        e_vc = e_v
        p_vc = p_v
        poe = p_vc/e_vc
        if(e_vc < 1e-5) poe=0
        zvir1p = TabInterpInterpolate(log10(Frhoc),e_vc,poe)
     endif
     m0=mp+1
     if(zvir1p.ge.ZZon) goto 300 ! if peak will collapse go to 300
     mupp_p=mupp
     zvir1=zvir1p
     Frhpk=Frhoh
     Fnupk=Frho
     Fevpk=e_v
     Fpvpk=p_v
     gradpk_p  = gradpk
     gradpkf_p = gradpkf
     strain_mat = Ebar
  enddo

299 continue ! if the shell closest in radius to Rf is the zeroth shell at m=1
             ! i.e. if the peak does not collapse by end redshift
  RTHL=-1.0
  Srb=0.0
  deallocate(nshell,Sshell,S2shell,Gshell,Gshellf,SRshell,Fbar,rad)

  return

300 continue ! jumps to here if peak is found to collapse by end redshift
  
  ! Now have localized zvir1p >= ZZon > zvir1, where zvir1p is the redshift
  ! of the last shell to collapse (shell # m0-1), ZZon is the redshift of
  ! the end of the simmulation (0.0 for non-lightcone runs) and zvir1 is
  ! the redshift of the first shell not to collapse (shell m0)
  dZvir=zvir1p-zvir1
  if(zvir1.gt.0.0.and.dZvir.ne.0.0) then
     ! RTHL=rad(m0-1)+(rad(m0)-rad(m0-1))*(zvir1p-ZZon)/dZvir
     RTHL3=rad(m0-1)**3+(rad(m0)**3-rad(m0-1)**3)*(zvir1p-ZZon)/dZvir ! interpolating RTHL3
     RTHL=RTHL3**(1./3.)
  else
     RTHL=rad(m0-1)
  endif
  if(RTHL.le.0.0) then
     RTHL=-1.0
     Srb=0.0
     deallocate(nshell,Sshell,S2shell,Gshell,Gshellf,SRshell,Fbar,rad)
     return
  endif
  RTHL3=RTHL**3 ! this seems unnecessary since RTHL is already defined as the cube root of RTHL3...
  RTHL5=RTHL3*RTHL*RTHL

  rad3p=rad(m0-1)**3
  rad3=rad(m0)**3
  dFbar=Fbar(m0)-Fbar(m0-1)
  drad3=rad3-rad3p
  if(zvir1.gt.0.0.and.drad3.ne.0.0) then ! Interpolates halo characteristics from values on last shell to collapse and the values found at the peak centre (maybe? ...)
     Fbarx      = Fbar(m0-1)+(RTHL3-rad3p)*dFbar/drad3              ! 
     Frhpk      = Frhoh+(RTHL3-rad3p)*(Frhpk-Frhoh)/drad3           ! not sure how Frhpk and Frhoh are different
     Fnupk      = Frho+(RTHL3-rad3p)*(Fnupk-Frho)/drad3             ! 
     Fevpk      = e_v+(RTHL3-rad3p)*(Fevpk-e_v)/drad3               ! 
     Fpvpk      = p_v+(RTHL3-rad3p)*(Fpvpk-p_v)/drad3               ! 
     strain_mat = Ebar+(RTHL3-rad3p)*(strain_mat-Ebar)/drad3        ! 
     gradpk     = gradpk + (RTHL3-rad3p)*(gradpk_p-gradpk)/drad3    ! 
     gradpkf    = gradpkf + (RTHL3-rad3p)*(gradpkf_p-gradpkf)/drad3 ! 
     !write(*,*) Frhpk, Frhoh
  else
     Fbarx      = Fbar(m0-1)
     Fnupk      = Frho
     Fevpk      = e_v
     Fpvpk      = p_v
     strain_mat = Ebar
     gradpk     = gradpk_p
     gradpkf    = gradpkf_p
  endif
  call normalize_strain(Fbarx,strain_mat)
  call get_evals(3,strain_mat,Lam,iflag,Evec)
  Frhoc = Lam(1) + Lam(2) + Lam(3)
  e_v = 0.5*(Lam(3) - Lam(1))/Frhoc
  p_v = 0.5*(Lam(3) + Lam(1) - 2.0*Lam(2))/Frhoc

  strain_mat = Ebar * Fbarx/Frhoh

  Sbar = 0
  Sbar2 = 0
  nSbar=0
  do m1=1,m0-1
     if(nshell(m1).gt.0) then
        Sbar  = Sbar + Sshell(:,m1)
        Sbar2 = Sbar2 + S2shell(:,m1)
        nSbar = nSbar+nshell(m1)
     endif
  enddo
  Sbar  = Sbar / float(nSbar)
  Sbar2 = Sbar2 / float(nSbar)

  !c  Integrate Fbar to get energy factor     
  Srb=0.0
  rad5p=0.0
  if(m0.gt.2) then
     do mp=2,m0-1
        rad5=rad(mp)**5
        Srb=Srb+0.5*(Fbar(mp-1)+Fbar(mp))*(rad5-rad5p)
        rad5p=rad5
     enddo
  endif
  if(zvir1.gt.0.0.and.dZvir.ne.0.0) &
       Srb=Srb+0.5*(Fbar(m0-1)+Fbarx)*(RTHL5-rad5p)

  Srb=Srb/(Fbarx*RTHL5)

  ! ----------------------------------------------------------
  ! Now make local measurements by looping over all cells such 
  ! that r<RTHL
  ! 

  ! Mask particles
  ! Sets mask(i,j,k)=1 if i,j,k is not in the buffer, leaves as 0 if it is in the buffer
  do jp=1,npart
     if(sqrt(float(irs2(jp)))>RTHL) cycle
     do L=1,3
        iv(L)=ipk(L)+ixsvec(jp,L)
        if(iv(L)<nbuff+1.or.iv(L)>nn(L)-nbuff) goto 87
     enddo
     mask(iv(1),iv(2),iv(3)) = 1 !max(mask(iv(1),iv(2),iv(3)),int(RTHL))
87   continue
  enddo


  if(ioutshear>=1) then
     Fd2pk=0
!     Fnupk=0
     naverage=0
     do jp=2,npart
        if(sqrt(float(irs2(jp)))>RTHL) cycle
        do L=1,3
           iv(L)=ipk(L)+ixsvec(jp,L)
           if(iv(L)<=1.or.iv(L)>=nn(L)) goto 88 ! break out of if statment if not in main volume or buffer volume
        enddo
     naverage = naverage + 1
!     Fnupk = Fnupk + delta(iv(1),iv(2),iv(3))
     Fd2pk = Fd2pk + lapd(iv(1),iv(2),iv(3))
     enddo
     if(naverage.ne.0) then
        Fd2pk = Fd2pk / naverage
!        Fnupk = Fnupk / naverage
     endif
  endif
88 continue

  ! calculate "formation redshift" = when halo of size Rhalo/2 collapsed
  ! goto 101 !uncomment if you want formation redshift
  zvir_half = -1
  if(RTHL/hlatt < 3) goto 101

  ebar_mean_full = 0
  ebar_mean_half = 0
  
  do jj=1,10
     rcur = RTHL * ( (jj-1) / 9. * 0.5 + 0.5 )**(1./3.)
     !c  Find mupp
     rupp=rcur+2.0*hlatt
     do m1=2,npart-1
        if(rad(m1).gt.rupp) then
           mupp=m1
           goto 33
        elseif(rad(m1)==0) then
           mupp=m1-1
           goto 33
        endif
     enddo
33   continue
     
     rlow=max(rcur-2.0*hlatt,0.0)
     do m1=1,npart-1
        if(rad(m1).gt.rlow) goto 34
     enddo
34   continue
     mlow=m1
     
     do m1=1,npart-1
        if(rad(m1).gt.rcur) goto 35
     enddo
35   continue
     m0=m1

     !c  Zvir at m0
     !c  Need to sum over particles on shells mlow to mupp
     do m1=mlow,mupp
        u0  = hlatt_1*(rad(m0)-rad(m1))
        u02 = u0*u0
        if(u02.lt.4.0) then
           u12 = hlatt_2*(rad(m0)**2+rad(m1)**2)
           wt  = hRinteg(u0,u12)/(u12-u02)**2
           do L=1,3           
              do K=1,3
                 ebar_mean_full(L,K) = ebar_mean_full(L,K)-&
                      0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
              enddo
           enddo
        endif
     enddo
     
     do L=1,3
        do K=1,3
           ebar_mean_full(L,K) = wRnor*ebar_mean_full(L,K)/rad(m0)
        enddo
     enddo
     !c  Solve for eigenvalues
     call get_evals(3,ebar_mean_full,Lam,iflag,Evec)
     !c  These should be in INCREASING order. 
     Frho = Lam(1) + Lam(2) + Lam(3)
     if(Frho.gt.0.0) then
        e_v = 0.5*(Lam(3) - Lam(1))/Frho
        p_v = 0.5*(Lam(3) + Lam(1) - 2.0*Lam(2))/Frho
     else
        e_v = 0.0
        p_v = 0.0
     endif
     !c  Renormalize to Fbar
     Frhoc = Fbar(m0)
     e_vc = e_v
     p_vc = p_v
     poe = p_vc/e_vc
     if(e_vc < 1e-5) poe=0
     zvir_half = max(TabInterpInterpolate(log10(Frhoc),e_vc,poe)-1,zvir_half)
  enddo

!  write(*,*) 'zform',RTHL*alatt,zvir_half

101 continue

  deallocate(nshell,Sshell,S2shell,Gshell,Gshellf,SRshell,Fbar,rad)
  return

end subroutine get_homel

function hRinteg(u0,u12)
  USE intreal_types

  !c Integral of SPH kernel 4pi*W(x)dx*x(x^2-u12) from x0 to 2

  real(sp) hRinteg,u0,u12,x0,x02,x04,h1,h2,one7

  hRinteg=0.0
  x0=abs(u0)
  if(x0.ge.2.0) return
  hRinteg=1.4*u12-31.0/35.0
  if(x0.eq.0.0) return
  x02=x0*x0
  x04=x02*x02
  one7=1.0/7.0

  !c value at upper boundary x=2
  hRinteg=1.6*u12-6.4*one7
  if(x0.ge.1.0) then
     h1=x04*(2.0-x0*(2.4-x0*(1.0-x0*one7)))
     h2=x02*(4.0-x0*(4.0-x0*(1.5-x0*0.2)))
     hRinteg=hRinteg-h2*u12+h1
     return
  else
     hRinteg=hRinteg-0.2*(u12-one7)
     if(x0.gt.0.0) then
        h1=x04*(1.0-x02*(1.0-3.0*x0*one7))
        h2=x02*(2.0-x02*(1.5-0.6*x0))
        hRinteg=hRinteg-h2*u12+h1
     endif
  endif

  return
end function hRinteg


function hinteg(u0)
  USE intreal_types

  !c Integral of SPH kernel 4pi*W(x)dx*x from x0 to 2

  real(sp) hinteg,u0,x0,x02,x04,h1

  hinteg=0.0
  x0=abs(u0)
  if(x0.ge.2.0) return
  hinteg=1.4
  if(x0.eq.0.0) return
  x02=x0*x0
  x04=x02*x02

  !c value at upper boundary x=2
  hinteg=1.6
  if(x0.ge.1.0) then
     h1=x02*(4.0-x0*(4.0-x0*(1.5-x0*0.2)))
     hinteg=hinteg-h1
     return
  else
     hinteg=hinteg-0.2
     if(x0.gt.0.0) then
        h1=x02*(2.0-x02*(1.5-0.6*x0))
        hinteg=hinteg-h1
     endif
  endif

  return
end function hinteg

subroutine get_ijk(n1xn2,n1,j,j1,j2,j3) ! this is the inverse operation to
  USE intreal_types                     ! assignment of array lagrange()
  !c Corrected version stm 3 Jun 1993   ! in subroutine get_pks()
  jj=j-1
  !c this used to be just j
  j3=1+jj/n1xn2
  je=jj-n1xn2*(j3-1)
  j2=1+je/n1
  j1=je-n1*(j2-1)+1
  !c this used to not have the +1 (since it was j not j-1)
  return
end subroutine get_ijk

subroutine lagrint(xa,val,est,n,mdim,x)
  USE intreal_types
  dimension xa(n),est(mdim),val(n,mdim)
  !C Lagrange 4 pt interpolation or 2 point 

  if(n.eq.4) then
     a1=(x-xa(2))/(xa(1)-xa(2))*(x-xa(3))/(xa(1)-xa(3))*(x-xa(4))/(xa(1)-xa(4))
     a2=(x-xa(1))/(xa(2)-xa(1))*(x-xa(3))/(xa(2)-xa(3))*(x-xa(4))/(xa(2)-xa(4))
     a3=(x-xa(1))/(xa(3)-xa(1))*(x-xa(2))/(xa(3)-xa(2))*(x-xa(4))/(xa(3)-xa(4))
     a4=(x-xa(1))/(xa(4)-xa(1))*(x-xa(2))/(xa(4)-xa(2))*(x-xa(3))/(xa(4)-xa(3))
     do m=1,mdim
        est(m)=a1*val(1,m)+a2*val(2,m)+a3*val(3,m)+a4*val(4,m)
     enddo
  else 
     a1=(x-xa(2))/(xa(1)-xa(2))
     a2 = 1-a1
     do m=1,mdim
        est(m) =a1*val(1,m)+a2*val(2,m)
     enddo
  endif
  return
end subroutine lagrint

subroutine checkiv(iv,nn,i)

  integer iv(3),nn(3)

  do L=1,3
     if(iv(L)<1.or.iv(L)>nn(L)) then
        write(*,*) 'call = ',i
        write(*,*) 'iv out of bounds:',iv(L),L,nn(L)
        call mpi_finalize(ierr)
        stop
     endif
  enddo
  
end subroutine checkiv

real function nth_order_derivative(f,order)

  implicit none

  integer order
  real f(-4:4),derivative

  if(order==8) then
     derivative = &
          1./280. * ( f(-4) - f(4) ) + &
          4./105. * (-f(-3) + f(3) ) + &
          1./5.   * ( f(-2) - f(2) ) + &
          4./5.   * (-f(-1) + f(1) ) 
  else
     derivative = &
          1./2. * (-f(-1) + f(1) ) 
  endif

  nth_order_derivative = derivative

  return 
  
end function nth_order_derivative

subroutine icloud(npart,nhunt)
  USE intreal_types
  use arrays  
  use mpivars

  integer, allocatable :: indx(:)

  ! THIS SUBROUTINE GENERATES Lattice POSITIONS 
  ! irs2 is the (r/alatt)^2 ixs,iys,izs are the positions in alatt units     

  n11=n1+1 ! 
  n21=n2+1 ! usually = nmesh+1
  n31=n3+1 ! 

  n1_21=n1/2+1 ! 
  n2_21=n2/2+1 ! usually 1/2 nmesh rounded up to nearest integer
  n3_21=n3/2+1 ! 

  irad2=nhunt**2 ! Here we set largest halo radius squared in lattice units
  if(nhunt>nhalomax) then
    if(myid==0) write(*,*) 'nhunt > nhalomax, exiting',nhunt,nhalomax
    call mpi_finalize(ierr)
    stop
  endif
  m=0
  do  i = 1,n31          !
     izz=(i-n3_21)       !
     do  j = 1,n21       ! cycles over all lattice sites in a tile
        iyy=(j-n2_21)    !
        do  k = 1,n11    !
           ixx=(k-n1_21) !
           irr=ixx*ixx+iyy*iyy+izz*izz ! distance (in lattice units) from the centre squared
           if(irr.le.irad2) then
              m = m + 1 ! if distance less than nhunt, m+=1
           endif
        enddo
     enddo
  enddo
  npart=m ! number of lattice sites within a radius nhunt of the centre of the simulaiton tile

  allocate(indx(npart+20))
  allocate(ixsvec(npart+20,3))
  allocate(irs2(npart+20))

  m=0
  do  i = 1,n31          !
     izz=(i-n3_21)       !
     do  j = 1,n21       ! cycles over all lattice sites in a tile
        iyy=(j-n2_21)    !
        do  k = 1,n11    !
           ixx=(k-n1_21) !
           irr=ixx*ixx+iyy*iyy+izz*izz
           if(irr.le.irad2) then
              m = m + 1                         ! if distance greater than nhunt, m+=1
              irs2(m) = ixx*ixx+iyy*iyy+izz*izz ! make list of radial distance squared
              indx(m)  = m                      ! make list of indices m
           endif 
        enddo    ! 
     enddo       ! Saves a map from number m in the list indx(m) to radial coordinate squared irs2(m)
  enddo          ! 

  call sort2(npart,irs2,indx) ! sorts irs2 from least to greatest and rearranges elements in indx so they still match irs2
  do i=1,npart
     irs2(indx(i))=i ! irs2 = (/ 1, 2, 3, ..., npart-2, npart-1, npart /)
  enddo              ! 

  mold=0
  do  i = 1,n31
    izz=(i-n3_21)
    do  j = 1,n21
       iyy=(j-n2_21)
       do  k = 1,n11
           ixx=(k-n1_21)
           irr=ixx*ixx+iyy*iyy+izz*izz
           if(irr.le.irad2) then
              mold = mold + 1
              m = irs2(mold)

              ixsvec(m,1) = ixx ! Creates an array ixsvec(m,:) such that indices
              ixsvec(m,2) = iyy ! m are ordered in the same way as indx(m)
              ixsvec(m,3) = izz ! 
           endif
        enddo
     enddo
  enddo

  do m=1,npart                    ! Restores irs2 to the radial
    irs2(m) = sum(ixsvec(m,:)**2) ! distance to the centre of a
  enddo                           ! tile squared

  deallocate(indx) ! nukes indx

  return
end subroutine icloud


subroutine atab4
  USE intreal_types
  use arrays
  parameter (two_thirds=2.0/3.0)
  pi=4.0*atan(1.0)

  !C  This subroutine provides the kernel akk and the derivative of the kernel,
  !C  dkk. The latter is divided by r. The grid density assignment kernel is akg.
  !C  The kernel akk is a generalization of schoenbergs m4 kernel with continuous 
  !C  derivatives up to the second. 

  ca=3./(2.*pi)
  cb=1./(4.*pi)

  do i=1,1100
     u=(i-1.)*0.01
     ru=sqrt(u)
     akk(i)=0.
     dkk(i)=0.
     if(ru.lt.1.)then
        akk(i)=ca*(two_thirds-u+0.5*u*ru)
        dkk(i)=ca*(-2.0+1.5*ru)
     else if(ru.ge.1.0.and.ru.le.2.0)then
        akk(i)=cb*(2.0-ru)**3
        dkk(i)=-3.0*cb*(2.0-ru)**2/ru
     end if
  enddo
  return
end subroutine atab4


subroutine get_evals(Nu,Qin,rdiag,iflag,Qout)
  USE intreal_types

  PARAMETER (nkrp=3)
  real(sp) Qin(nkrp,nkrp),rdiag(nkrp),Qout(nkrp,nkrp)
  real(dp) evec(nkrp,nkrp)
  real(dp) diag(nkrp),offdiag(nkrp)
  save nuphys
  data nuphys/3/

  iflag=0

  do ib=1,Nu
     do ia=1,Nu
        evec(ia,ib)=dble(Qin(ia,ib))
     enddo
  enddo

  call TRED2(evec,Nu,nuphys,diag,offdiag)
  call TQLInm(diag,offdiag,Nu,nuphys,evec,iflag)

  !C Kth COLUMN OF evec IS THE NORMALIZED EVECTOR OF EVALUE diag(K)
  !C   THAT IS   evec(ia,K),ia=1,3 
  !c IFLAG=0 normally, IFLAG=-1 for bad solution

  if (iflag.eq.0) then
     do ib=1,Nu
        rdiag(ib)=sngl(diag(ib))
        do ia=1,Nu
           Qout(ia,ib)=sngl(evec(ia,ib))
        enddo
     enddo
     call eval_sort(Nu,nuphys,rdiag,Qout)
  else
     !c TQLI failed
     do ib=1,Nu
        rdiag(ib)=-1
     enddo
  endif

  return
end subroutine get_evals


SUBROUTINE TRED2(A,N,NP,D,E)
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  REAL(DP) A(NP,NP),D(NP),E(NP)
  IF(N.GT.1)THEN
     DO I=N,2,-1  
        L=I-1
        H=0.0d0
        SCALE=0.0d0
        IF(L.GT.1)THEN
           DO K=1,L
              SCALE=SCALE+ABS(A(I,K))
           enddo
           IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
           ELSE
              DO K=1,L
                 A(I,K)=A(I,K)/SCALE
                 H=H+A(I,K)**2
              enddo
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.0d0
              DO J=1,L
                 A(J,I)=A(I,J)/H
                 G=0.0d0
                 DO K=1,J
                    G=G+A(J,K)*A(I,K)
                 enddo
                 IF(L.GT.J)THEN
                    DO K=J+1,L
                       G=G+A(K,J)*A(I,K)
                    enddo
                 ENDIF
                 E(J)=G/H
                 F=F+E(J)*A(I,J)
              enddo

              HH=F/(H+H)
              DO J=1,L
                 F=A(I,J)
                 G=E(J)-HH*F
                 E(J)=G
                 DO K=1,J
                    A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                 enddo
              enddo
           ENDIF
        ELSE
           E(I)=A(I,L)
        ENDIF
        D(I)=H
     enddo
  ENDIF
  D(1)=0.0d0
  E(1)=0.0d0
  DO I=1,N
     L=I-1
     IF(D(I).NE.0.)THEN
        DO J=1,L
           G=0.0d0
           DO K=1,L
              G=G+A(I,K)*A(K,J)
           enddo
           DO K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
           enddo
        enddo
     ENDIF
     D(I)=A(I,I)
     A(I,I)=1.0d0
     IF(L.GE.1)THEN
        DO J=1,L
           A(I,J)=0.0d0
           A(J,I)=0.0d0
        enddo
     ENDIF
  enddo
  RETURN
END SUBROUTINE TRED2

SUBROUTINE TQLInm(D,E,N,NP,Z,iflag)
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  REAL(DP) D(NP),E(NP),Z(NP,NP)
  real(dp), parameter :: EPS=epsilon(x)
  iflag=0
  IF (N.GT.1) THEN
     DO I=2,N
        E(I-1)=E(I)
     enddo
     E(N)=0.0d0
     DO L=1,N
        ITER=0
1       DO M=L,N-1
           DD=ABS(D(M))+ABS(D(M+1))
           ! SEE http://www.nr.com/forum/showthread.php?p=5251
           ! IF (ABS(E(M))+DD.EQ.DD) GO TO 2 ! OLD WAY
           if(ABS(E(M))<=EPS*DD) GO TO 2     ! NEW WAY 
        enddo
        M=N
2       IF(M.NE.L)THEN
           !c            IF(ITER.EQ.30)PAUSE 'too many iterations'
           IF(ITER.EQ.60) THEN
              write(*,*) 'TQLI : too many iterations = ',ITER
              iflag=-1
              return
           ENDIF
           ITER=ITER+1
           G=(D(L+1)-D(L))/(2.0d0*E(L))
           R=SQRT(G**2+1.0d0)
           G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
           S=1.0d0
           C=1.0d0
           P=0.0d0
           DO I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                 C=G/F
                 R=SQRT(C**2+1.0d0)
                 E(I+1)=F*R
                 S=1.0d0/R
                 C=C*S
              ELSE
                 S=F/G
                 R=SQRT(S**2+1.0d0)
                 E(I+1)=G*R
                 C=1.0d0/R  
                 S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.0d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO K=1,N
                 F=Z(K,I+1)
                 Z(K,I+1)=S*Z(K,I)+C*F
                 Z(K,I)=C*Z(K,I)-S*F
              enddo
           enddo
           D(L)=D(L)-P
           E(L)=G
           E(M)=0.0d0
           GO TO 1
        ENDIF
     enddo
  ENDIF
  RETURN
END SUBROUTINE TQLInm

subroutine eval_sort(Nu,nuphys,diag,evec)
  USE intreal_types
  !C SORTS THE EVALUES SO BIGGEST IN 3, SMALLEST IN 1
  !C SORTS THE EVECTORS AS WELL 
  parameter (n=3)
  real(sp) diag(nuphys),evec(nuphys,nuphys)
  real(sp) wksp(n)
  integer(i4b) iwksp(n)

  call rindexx(Nu,diag,iwksp)

  !C THE NEXT STEP CHANGES TO DECREASING ORDER FROM INCREASING ORDER
  !c      Nu_2=Nu/2
  !c      Nu1=Nu+1
  !c      do 10 j=1,Nu_2
  !c         idum=iwksp(j)
  !c         iwksp(j)=iwksp(Nu1-j)
  !c         iwksp(Nu1-j)=idum
  !c 10   continue

  do j=1,Nu
     wksp(j)=diag(j)
  enddo
  do j=1,Nu
     diag(j)=wksp(iwksp(j))
  enddo
  do j=1,Nu
     wksp(j)=evec(1,j)
  enddo
  do j=1,Nu
     evec(1,j)=wksp(iwksp(j))
  enddo
  do j=1,Nu
     wksp(j)=evec(2,j)
  enddo
  do j=1,Nu
     evec(2,j)=wksp(iwksp(j))
  enddo
  do j=1,Nu
     wksp(j)=evec(3,j)
  enddo
  do j=1,Nu
     evec(3,j)=wksp(iwksp(j))
  enddo
  return
end subroutine eval_sort

subroutine normalize_strain(F,strain)

  implicit none
  real             :: strain(3,3)
  real, intent(in) :: F

  strain = strain * F / (strain(1,1)+strain(2,2)+strain(3,3))

  return

end subroutine normalize_strain

SUBROUTINE sort2(n,arr,brr)         ! This subroutine sorts the first n elements
      INTEGER n,M,NSTACK            ! of array arr from least to greatest and
      INTEGER(kind=4) arr(n),brr(n) ! rearranges array brr so as to maintain the
      PARAMETER (M=7,NSTACK=50)     ! relationship arr(i):brr(i).
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then           ! if n-1 < 7
        do 12 j=l+1,ir            !     do j=2,n
          a=arr(j)                !         a=arr(j)
          b=brr(j)                !         b=brr(j)
          do 11 i=j-1,l,-1        !         do i=j-1,1,-1
            if(arr(i).le.a)goto 2 !             if
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          temp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
END SUBROUTINE SORT2

SUBROUTINE RINDEXX(N,ARRIN,INDX)
  USE intreal_types
  REAL(SP) ARRIN(*),Q
  INTEGER(I4B) INDX(*)
  DO J=1,N
     INDX(J)=J
  enddo
  L=N/2+1
  IR=N
10 CONTINUE
  IF(L.GT.1)THEN
     L=L-1
     INDXT=INDX(L)
     Q=ARRIN(INDXT)
  ELSE
     INDXT=INDX(IR)
     Q=ARRIN(INDXT)
     INDX(IR)=INDX(1)
     IR=IR-1
     IF(IR.EQ.1)THEN
        INDX(1)=INDXT
        RETURN
     ENDIF
  ENDIF
  I=L
  J=L+L
20 IF(J.LE.IR)THEN
     IF(J.LT.IR)THEN
        IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
     ENDIF
     IF(Q.LT.ARRIN(INDX(J)))THEN
        INDX(I)=INDX(J)
        I=J
        J=J+J
     ELSE
        J=IR+1
     ENDIF
     GO TO 20
  ENDIF
  INDX(I)=INDXT
  GO TO 10
END SUBROUTINE RINDEXX
