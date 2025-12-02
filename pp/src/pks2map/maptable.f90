module maptable

  use cosmology
  use bbps_profile
  use integrate_profiles
  use profiles

  implicit none

  ! Table bounds and interval
  real   chit,   chimint,   chimaxt
  real    mht,    mhmint,    mhmaxt
  real     rt,     rmint,     rmaxt

  ! Table sizes and intervals
  integer nrt, nchit, nmht, nprofiles
  real    drt, dchit, dmht

  parameter (nprofiles=3)

  ! Table
  real, allocatable :: table(:,:,:,:)

contains

  subroutine makemaptable(fileout,model)

    implicit none 

    character *512 fileout
    character *512 filein
    integer model

    ! Compton y-parameter
    real y

    ! Thomson scattering optical depth
    real tau

    ! Integrated gas
    real Sigmatilde_gas
    ! Integrated DM
    real Sigmatilde_dm

    ! Other variables
    integer i,j,k
    real z,theta,Delta,zfacy,zfact,zfactdm,m2rfac,mht3
    real r200c

    ! Debugging
    real tau_mean, tau_mean_expected, thetap, dtheta

    ! Create distance-redshift tables
    call set_rofztable
    call set_zofrtable
  
    bbps_model = model
    Delta = bbps_delta(bbps_model)
    
    m2rfac = (3./4./pi/Delta)**(1./3.) 

    y0      = fb * Delta * h**2 * ypermass  ! in 1/Msun
    y0      = y0/2 ! to undo for the multiplicative factor of 2 in the generic integrate_profile(...),
                   ! As P0_Am has already accounted for this (BBPS eq. 9). 
                   ! This is not true for generic normalizations

    ! Setup table bounds
    
    ! Distance
    nchit   = 200
    chimint = 1.842
    chimaxt = 1.393e4
    
    dchit = (log(chimaxt)-log(chimint))/(nchit-1)
    
    ! Mass
    nmht   = 200
    mhmint = 1e10
    mhmaxt = 1e16
    
    dmht = (log(mhmaxt)-log(mhmint))/(nmht-1)

    ! Radius
    nrt   = 200 
    rmint = 1e-2
    rmaxt = 4
    
    drt = (log(rmaxt)-log(rmint))/(nrt-1)

    allocate(table(nprofiles,nrt,nmht,nchit))
    
    ! Calculate tabulated values and output
    
    ! Edit from April 2024 - trying to make it work on gcc and following modern standard
    ! (form='binary' was depreciated)
    !open(unit=1, file=filein, form='binary')
    open(unit=1, file=filein, form='unformatted', access='stream', status='replace')

    write(1) nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt
    
    do i=1,nchit
       chit    = exp( log(chimint) + (i-1) * dchit )
       z       = zofr(chit)
       zfacy   = (omegam*(1+z)**3+omegal)

       do j=1,nmht
          mht  = exp( log(mhmint) + (j-1) * dmht )
          mht3 = mht**(1./3.)

          do k=1,nrt
             rt = exp( log(rmint) + (k-1) * drt )
             theta = atan(rt/chit) ! NJC 6apr23 - should that be asin(rt/chit) ? Because the angle subtended by sphere of radius R
             ! centred at x is arcsin(R/x) not arctan(R/x) though they converge for x >> R

             y   =   y0 * mht  * zfacy * integrate_profile(theta,z,mht,1)

             Sigmatilde_gas = integrate_profile(theta,z,mht,2)
             Sigmatilde_dm  = integrate_profile(theta,z,mht,3) 
             
             table(1,k,j,i) = y 
             table(2,k,j,i) = fb*Sigmatilde_gas
             table(3,k,j,i) = Sigmatilde_dm

          enddo
       enddo
    enddo

    write(1) table
    
    close(1)
    
  end subroutine makemaptable
  
  !=======================================================================

  subroutine loadmaptable(filein)

    implicit none

    character *512 filein

    ! Edit from April 2024 - trying to make it work on gcc and following modern standard
    ! (form='binary' was depreciated)
    !open(unit=1, file=filein, form='binary')
    open(unit=1, file=filein, form='unformatted', access='stream', status='replace')

    read(1) nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt

    allocate(table(nprofiles,nrt,nmht,nchit))

    read(1) table

    close(1)

    dchit = (log(chimaxt) - log(chimint)) / (nchit - 1)
    drt   = (  log(rmaxt) -   log(rmint)) / (nrt - 1)
    dmht  = ( log(mhmaxt) -  log(mhmint)) / (nmht - 1)

    return

  end subroutine loadmaptable

  !=======================================================================

  real function interpolate_table(theta,z,mh)

    implicit none

    real theta,z,chi,r,mh,pad
    real fr,fc,fm
    integer ir,ichi,imh

    chi = rofz(z)
    r   = chi * tan(theta)

    interpolate_table = 0
    if(log(r)>log(rmaxt)-drt .or. log(chi)>log(chimaxt)-dchit) return

    if(   r < rmint   )   r =   rmint + 1e-5
    if( chi < chimint ) chi = chimint + 1e-5
    if(  mh < mhmint  )  mh =  mhmint + 1e-5

    if(  mh > mhmaxt  )  mh =  mhmaxt - 1e-5

    ir   = int( (   log(r) -   log(rmint) ) / drt   ) + 1
    ichi = int( ( log(chi) - log(chimint) ) / dchit ) + 1
    imh  = int( (  log(mh) -  log(mhmint) ) / dmht  ) + 1

    fr = log(r)   - (   log(rmint) + (  ir - 1) *   drt )
    fc = log(chi) - ( log(chimint) + (ichi - 1) * dchit )
    fm = log(mh)  - (  log(mhmint) + ( imh - 1) *  dmht )
    
    if(ir>=nrt.or.ichi>=nchit.or.imh>=nmht) then
       write(*,*) 'error: ir,ichi,imh',ir,&
            r,rmint,rmaxt,drt
       stop
    endif
    interpolate_table = &
         table(profile,ir  ,imh  ,ichi  ) * (1-fr) * (1-fc) * (1-fm) + &
         table(profile,ir+1,imh  ,ichi  ) * (  fr) * (1-fc) * (1-fm) + &
         table(profile,ir  ,imh+1,ichi  ) * (1-fr) * (  fc) * (1-fm) + &
         table(profile,ir+1,imh+1,ichi  ) * (  fr) * (  fc) * (1-fm) + &
         table(profile,ir  ,imh  ,ichi+1) * (1-fr) * (1-fc) * (  fm) + &
         table(profile,ir+1,imh  ,ichi+1) * (  fr) * (1-fc) * (  fm) + &
         table(profile,ir  ,imh+1,ichi+1) * (1-fr) * (  fc) * (  fm) + &
         table(profile,ir+1,imh+1,ichi+1) * (  fr) * (  fc) * (  fm) 
    
    interpolate_table = interpolate_table 
            
  end function interpolate_table

end module maptable
