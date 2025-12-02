program pks2map

  ! This code reads in a peak patch catalog and outputs a map.
  ! For Pasted Profiles the contribution from each halo to a given pixel is 
  ! stored in a 3D table as a function of halo mass, redshift, and l.o.s. angle

  use cosmology
  use line_profile
  use pksc
  use haloproject
  use profiles

  use mpivars
  use fitsvars
  use healpixvars
  use flatskyvars

  use fitstools
  use healpix_types
  use head_fits

  use textlib
  
  use random

  implicit none

  ! Filenames
  character *512 filein, fileout, fileout_bn,fileout_fs,fileout_hp,&
       fileout_dp,tablefile,freqout
  character *4 outcode
  integer ilast, ilastnu

  ! Halo properties
  real xh,yh,zh,mh,rh,redshifth,chih,vh,mfh
  real xmmax,xmmaxl, ymmax,ymmaxl, zmmax,zmmaxl, rmmax

  ! dPdy
  integer, parameter :: ndpdy = 10000
  real,    parameter :: dpdy_max = 1e-4,dpdy_dy=dpdy_max/ndpdy
  real dpdy(ndpdy)

  ! Other variables
  integer i,j,k,jmin,jmax,kmin,kmax,m

  real mmin,chihview,zview,cut_low,cut_high,zmin,zmax,mmax
  real mypi
  real RTHmaxl,RTHmax,rmaxtheta,rmaxphi

  integer nhalotot, nhalomap 
  integer model, centre, scramble, kappa, PSZcut

  ! profile
  character *3 profilecode

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  ! Usage check
  if(command_argument_count()<3) then
     if(myid==0) write(*,11) 
11   format('Usage: pks2map <filein> <fileout> <tablefile> <profilecode> [<zmin> <nside> <scramble> <center> <npix> <fov> <zmax> <chihview> <model> <PSZcut>]  profile choices: tsz, ksz, kap, tau, tco, thi')
     call mpi_finalize(ierr)
     stop
  endif

  ! Get commandline
  call      get_command_argument(1,filein)
  call      get_command_argument(2,fileout)
  call      get_command_argument(3,tablefile)
  call      get_command_argument(4,profilecode)
  zmin     = r4arg(5,0.0)
  hpnside  = i4arg(6,4096)
  scramble = i4arg(7,0)
  centre   = i4arg(8,0) 
  npix     = i4arg(9,16)
  fov      = r4arg(10,10.0)
  zmax     = r4arg(11,1.25)
  chihview = r4arg(12,0.0)
  model    = i4arg(13,1)
  PSZcut   = i4arg(14,0)

  if(profilecode == 'tsz') profile = itsz ! tSZ map
  if(profilecode == 'ksz') profile = iksz ! kSZ map
  if(profilecode == 'kap') profile = ikap ! kappa lensing convergence map
  if(profilecode == 'tau') profile = itau ! tau optical depth to polarization
  if(profilecode == 'tco') profile = itco ! T_CO carbon monoxide line temperature
  if(profilecode == 'thi') profile = ithi ! T_HI 21-cm hydrogen line temperature

  mmin = 2.5e10

  if (profile == itco) then
     !set CO line emission parameters
     nu_rest = 115.27    ! rest frame frequency in GHz (CO(1-0)=115.27GHz)
     dnu     = 0.0078125 ! bandwidth in GHz
     nuobs_i = 34.0      ! frequency range in GHz
     nuobs_f = 30.0
  elseif (profile == ithi) then
     !set HI line emission parameters
     nu_rest = 1.420405751786  ! rest frame frequency in GHz 
     dnu     = 9.626953125e-05 ! 0.0004929 bandwidth in GHz to match 100 CO slices
     nuobs_i = 0.41896         ! frequency range in GHz to match CO
     nuobs_f = 0.36967
  elseif (profile == ikap) then
     kappa = 1
  endif

  ! number of maps to create
  if ((profile ==  itco) .or. (profile ==  ithi)) then
     num_nuobs = int((nuobs_i - nuobs_f)/dnu) ! = 512 for hardcoded values of nu's above
  else
     num_nuobs = 1
  endif

  testcase = .false.
  mypi = 4.*atan(1.0)

  ! Set field of view and resolution
  fov     = fov / 360 * 2 * mypi ! radians
  dpix    = fov     / npix
  Ompix   = (dpix)**2  ! pixel size to convert to brightness temp

  if(myid==0) write(*,12) npix,npix,dpix,fov/2/mypi*360.
12 format(/,/,&
          'Resolution:         ',i0,' x ',i0,' pixels (',1pe8.2,' radians)',/,&
          'Field of view:      ',1pe8.2,' degrees')


  ! Healpix map resolution
  hpnpix = nside2npix(hpnside) ! 12 * hpnside

  if(myid==0) write(*,13) hpnside,hpnpix
13 format('Healpix Nside:      ',i0,/,&
          'Number of pixels:   ',i0)

  ! Allocate maps  
  if(myid==0) allocate(fsmap(npix*npix))
  allocate(fsmapl(npix*npix))

  if(myid==0) allocate(hpmap(0:hpnpix-1,1:1)) 
  allocate(hpmapl(0:hpnpix-1,1:1))
  allocate(hplist(0:hpnpix-1))

  ! Set background cosmology 
  omegam = 0.3138 ! Matter energy density fraction
  omegab = 0.0493 ! Baryonic matter energy density fraction
  omegal = 0.6862 ! Cosmological constant energy density fraction
  h      = 0.6735 ! Dimensionless Hubble constant H_0/(100 km/s/Mpc)
  sigma8 = 0.8111 ! Standard deviation of delta at z=0 smoothed at R=8 Mpc/h
  ns     = 0.9649 ! scalar spectral index n_s
  w      = -1     ! equation of state w=p/rho
  fb     = omegab/omegam ! baryon fraction

  rhocrit_0 = 2.775e11*h**2    ! critical energy density at z=0 3H_0^2/8piG [M_sol/Mpc^3]
  rho_0     = rhocrit_0*omegam ! average matter energy density at z=0

  ! Create distance-redshift tables
  call set_rofztable
  call set_zofrtable

  cut_low  = rofz(zmin)
  cut_high = rofz(zmax)

  if(myid==0) write(*,14) cut_low/1e3,cut_high/1e3
14 format('cut_low, cut_high:  ',f4.2,1x,f4.2,' Gpc')

  ! Load halos
  call loadpksc(filein)

  ! Cut all the halos not within range
  m=0
  do i=1,nhalo 
     xh  = posxyz(1,i)
     yh  = posxyz(2,i) 
     zh  = posxyz(3,i)           
     xh = xh - chihview

     chih = sqrt(xh**2+yh**2+zh**2)
     if(chih>cut_high.or.chih<cut_low) cycle

     m = m+1
     posxyz(3,m) = posxyz(3,i)
     posxyz(2,m) = posxyz(2,i)
     posxyz(1,m) = posxyz(1,i)
     rth(m)      = rth(i)
  enddo
  nhalo = m

  !Center and Scramble halo switch
  if(scramble==1) call scramble_halos ! scramble halos in theta and phi
  if(centre==1  ) call center_halos   ! recentre halos around the most massive object in the catalogue
  
  mass = 4./3.*mypi * rho_0 * rth**3  ! approximate mass within Lagrangian sphere of halo

  ! Load profile table
  if (profile==itco)then
     call SFRtabload(tablefile)   ! load SFR Behroozi table
  else
     call loadmaptable(tablefile) ! load BBPS table created from make_maptable
     if (profile==ikap) then
        ! table has units of rhobar. Need to subtract 1 to get delta 
        ! add (gas profile) + (DM PROFILE) together in table level to make kappa_tot
        table(3,:,:,:) = (table(2,:,:,:) + table(3,:,:,:))-1
     endif
     if (profile==itau) profile=iksz
  endif

  ! theta0 and phi0 are the center of fov
  ! here we set the center to lie along z-axis
  theta0 = 0
  phi0   = 0
  theta1 = mypi/2 - theta0
  
  ! Loop over halos and project each into map
  do j=1,num_nuobs ! loop over 512 frequency bands for CO maps (otherwise num_nuobs=1)
  fsmapl=0.
  hpmapl=0.
  meantaul=0.
  maxtheta=0.
  m=0
  nuobsj = nuobs_i - (j-1./2)*dnu
  if( (profile == itco .or. profile == ithi) .and. myid==0) write(*,*) "nuobs = ", nuobsj

  do i=1,nhalo

     if(mod(i,10000)==0) write(*,*) i,nhalo,real(i)/real(nhalo)

     xh  = posxyz(1,i) ! 
     yh  = posxyz(2,i) ! position of halo i from catalogue
     zh  = posxyz(3,i) ! 

     vh  = vrad(i) ! halo radial velocity in units of v/c (determined in subroutine loadpksc from pksc.f90)

     chih = sqrt(xh**2+yh**2+zh**2)
     redshifth=zofr(chih)
     xh = xh - chihview
     chih = sqrt(xh**2+yh**2+zh**2)

     !Mass Conversions and Halo Selection
     if ((profile==itco) .or. (profile ==  ithi)) then
        nuobsi = nu_rest/( (1+redshifth) *(1+vh))   ! observed frequency of halo
        if ( abs(nuobsj-nuobsi)>dnu/2 ) cycle  ! only use halos in frequency slice
        mh = mass(i) ! use m200
     else ! estimate mass based on BBPS model and earlier estimate of halo volume * \bar{rho}
        mh = sqrt(deltacrit(redshifth)/bbps_delta(model)) * mass(i) ! Convert to m500 assuming SIS
     endif

     !Cut halos below polynomial fit to Planck 2015 selection function
     !From figure 26 of Planck 2015 XXVII (50% survey completeness contour)
     if (PSZcut == 1) mmin = sqrt(500./200)*(redshifth**3 - 2.73*redshifth**2 + 2.52*redshifth + 0.11)*6.78e14 
    
     !Mass cut
     if(mass(i)<mmin) cycle
     m = m+1
     rh  = (3*mh/4/mypi/delta_pksc/rho_0)**(1./3.) !comoving virial radius

     ! Multiplication Factors for different profiles
     mfh = 1.0 
     if (profile==itsz .or. profile==ikap) vh = 1.0
     ! Integrated Lensing Kernel for Kappa  
     if (kappa==1) mfh = rh*W_Kappa(redshifth,chih,14.2e3) !hardwired to lensing of CMB
     ! tau0persigma * W_tau^tilde * Mvir^(1/3) * Sigmatilde_gas. negative for -vh/c
     if (profile==iksz) mfh = -tau0persigma*(omegam*h**2)**(2./3)*(1+redshifth)**2*(mh/bbps_delta(model))**(1./3)
     
     ! Project halo onto the HEALpix skymap using subroutines from src/pks2cmb/haloproject.f90
     if(chih>rh) then ! exclude local halo if one exists
        call projecthalo_healpix(xh,yh,zh,rh,chih,mh,redshifth,vh,mfh)
        call projecthalo_flatsky(xh,yh,zh,rh,chih,mh,redshifth,vh,mfh)
     endif

  enddo

  nhalol        = m
  maxthetal     = maxtheta
  maxthetamassl = maxthetamass
  maxthetachihl = maxthetachih
  call mpi_reduce(nhalol,nhalomap,1,mpi_int,mpi_sum,0,mpi_comm_world,ierr) ! sum number of halos from each parallel task
  call mpi_reduce(meantaul,meantau,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr) ! sum mean tau from each parallel task

  ! Gather maps from processors
  call mpi_reduce(fsmapl,fsmap,npix**2,mpi_real,mpi_sum,0,mpi_comm_world,ierr) ! sum npix^2 elements of flat sky map from each parallel task
  call mpi_reduce(hpmapl,hpmap,hpnpix, mpi_real,mpi_sum,0,mpi_comm_world,ierr) ! sum hpnpix elements of HEALPix map from each parallel task

  ! Output maps from processor 0
  if(myid==0) then

  write(*,15) nhalomap
15 format('Number of halos in map:    ',i0)

  if(profile.ne.itsz.and.profile.ne.itco.and.profile.ne.ithi) write(*,16) meantau/hpnpix
16 format('Mean tau:    ',1pe9.3)

  fsmom1 = sum(fsmap**1) /npix / npix
  fsmom2 = sum(fsmap**2) /npix / npix
  fsmom3 = sum(fsmap**3) /npix / npix

  hpmom1 = sum(hpmap**1) / hpnpix
  hpmom2 = sum(hpmap**2) / hpnpix
  hpmom3 = sum(hpmap**3) / hpnpix
  
  write(*,17) fsmom1,fsmom2,fsmom3,hpmom1,hpmom2,hpmom3
17 format('Flatsky moments:    ',3(1pe9.3,1x),/,&
          'Healpix moments:    ',3(1pe9.3,1x))

  write(*,18) sqrt(sum((fsmap-fsmom1)**2)/npix**2),fsmom1,&
              sqrt(sum((hpmap-hpmom1)**2)/hpnpix) ,hpmom1
18 format('Flatsky RMS:                ',1pe9.3,/,&
          'Flatsky Mean:               ',1pe8.2,/,&
          'Healpix RMS:                ',1pe9.3,/,&
          'Healpix Mean:               ',1pe8.2,/)
     

  ! Get dp/dy
  dpdy=0
  do i=0,hpnpix-1
     m = int(hpmap(i,1)/dpdy_dy)+1
     if(m>0.and.m<=ndpdy) dpdy(m)=dpdy(m)+1
  enddo
  dpdy = dpdy / real(hpnpix)

  ! Check mean
  fsmom1=0.0
  do i=1,ndpdy
     fsmom1=fsmom1+(i-0.5)*dpdy_dy*dpdy(i)
  enddo

  do i=2,ndpdy
     dpdy(i)=dpdy(i-1)+dpdy(i)
  enddo  

  ! Create filenames 
  ilast   = indlnb(fileout) !index last non-blank character in fileout (from src/modules/External/textlib.f90)
  write (freqout, "(I0.3)") j
  if ((profile==itco) .or. (profile ==  ithi)) then
     fileout_fs=fileout(1:ilast)//'_fs.map'//trim(freqout)
     fileout_hp=fileout(1:ilast)//'_hp.fits'//trim(freqout)
     fileout_dp=fileout(1:ilast)//'_dp.bin'//trim(freqout)
  else
     fileout_fs=fileout(1:ilast)//'_fs.map'
     fileout_hp=fileout(1:ilast)//'_hp.fits'
     fileout_dp=fileout(1:ilast)//'_dp.bin'
  endif

  ! P(<y) file
  open(unit=1,file=fileout_dp,form='binary')
  write(1) ndpdy,dpdy_max
  write(1) dpdy
  close(1)

  ! Flat sky binary file
  open(unit=1,file=fileout_fs,form='binary')
  write(1) npix,npix,fov,fov
  write(1) fsmap
  close(1)
     
  ! Healpix FITS file
  hpheader(:)=''
  call add_card(hpheader,'NSIDE',hpnside,'the nside of the map')
  call add_card(hpheader,'ORDERING','RING','the nside of the map')
  call output_map(hpmap,hpheader,fileout_hp)

  endif

  enddo !end loop over freq for profile==itco and profile==ithi
  
  if(myid==0) write(*,*)

  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_finalize(ierr)

  stop

81 format('000',i1)
82 format( '00',i2)
83 format(  '0',i3)
84 format(      i4)

end program pks2map
