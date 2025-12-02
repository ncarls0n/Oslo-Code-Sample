module haloproject

  use cosmology

  use mpivars
  use flatskyvars
  use healpixvars
  use maptable

  use bbps_profile
  use line_profile
  use profiles

  use pix_tools
  use healpix_types
  use head_fits

  implicit none

  real, parameter :: rmax = 4.0
  double precision theta_halo,phi_halo,dtheta_halo
  real angle

  ! Halo angular properties
  real maxtheta,maxthetamass,maxthetachih  
  real maxthetal,maxthetamassl,maxthetachihl  

  ! Map moments
  double precision fsmom1,fsmom2,fsmom3
  double precision hpmom1,hpmom2,hpmom3
  real curval

  real interpval, sfri
  ! Weird stuff
  real localpi

contains

!=======================================================================

  subroutine projecthalo_healpix(x,y,z,r,chi,m,redshift,v,mf)

    implicit none

    real x,y,z,r,chi,m,redshift,v,mf
    integer i,j,profile

    localpi = 4*atan(1.0)

    dtheta_halo = rmax*asin(r/chi) ! note rmax is parameter set to 4.0
    if (profile==itco .or. profile==ithi) dtheta_halo = 0.    

    ! Get unit vector    
    hpvec1(1) = x / chi
    hpvec1(2) = y / chi
    hpvec1(3) = z / chi

    ! For Pasted Luminosities
    if ((profile==itco) .or. (profile==ithi)) then
       call vec2pix_ring(hpnside,hpvec1,j)
       if (profile==itco) then
          sfri   = 10**SFR(log10(m),redshift)
          curval = L_CO(sfri,redshift)
       elseif (profile==ithi) then
          curval = MHtoMHI(m,redshift)
          curval = L_HI(curval)
       endif
       curval = I_line(curval,redshift,chi,dnu)
       curval = T_line(curval,nuobsj)

       meantaul    = meantaul + curval
       hpmapl(j,1) = hpmapl(j,1) + curval
    else
    ! For Pasted Profiles call
    call query_disc(hpnside,hpvec1,dtheta_halo,hplist,nlist) ! subroutine in src/modules/External/pixel_routines.f90, returns list of HEALPix pixels hplist that the halo overlaps and the length of that list nlist
    do i=0,nlist-1 ! Loop over all pixels in HEALPix map

       j=hplist(i)
       call pix2vec_ring(hpnside,j,hpvec2) ! sets hpvec2 to the vector corresponding the the jth pixel in the HEALPix map
       call angdist(hpvec1,hpvec2,hpangle) ! sets hpangle as the angle between the jth pixel and the angular position of the halo centre at (x,y,z)/chi
       
       angle = hpangle

       if(j>hpnpix-1.or.j<0) write(*,*) 'error',i,j,hpnpix,nlist
       
       if(angle<dtheta_halo) then ! if the pixel j is within the halo, add the intensity
          interpval   = interpolate_table(angle,redshift,m)
          curval      = mf * v * interpval 
          meantaul    = meantaul + interpval
          hpmapl(j,1) = hpmapl(j,1) + curval
       endif

    enddo
    endif
    return

  end subroutine projecthalo_healpix

!=======================================================================

  subroutine projecthalo_flatsky(x,y,z,r,chi,m,redshift,v,mf)

    implicit none

    real x,y,z,r,chi,m,redshift,v,mf

    integer i,j,k,jmin,jmax,kmin,kmax


    localpi = 4*atan(1.0)

    dtheta_halo = rmax*asin(r/chi)
    if (profile==itco) dtheta_halo = 0.

    ! Get unit vector
    x = x / chi
    y = y / chi
    z = z / chi

    if(z<0.0) return

    theta_halo = asin(x) 
    phi_halo   = asin(y)

    thetamin = theta_halo - 1.1 * dtheta_halo
    thetamax = theta_halo + 1.1 * dtheta_halo
    
    phimin   = phi_halo - 1.1 * dtheta_halo
    phimax   = phi_halo + 1.1 * dtheta_halo
    
    if(thetamax < -fov/2 .or. thetamin > fov/2) return
    if(  phimax < -fov/2 .or.   phimin > fov/2) return

    ! For Pasted Luminosities
    if ((profile==itco) .or. (profile==ithi)) then
       if(theta_halo < -fov/2 .or. theta_halo > fov/2) return
       if(  phi_halo < -fov/2 .or.   phi_halo > fov/2) return

       j = int((theta_halo + fov/2)/dpix) + 1
       k = int((phi_halo + fov/2)/dpix) + 1

       if (profile==itco) then
          sfri   = 10**SFR(log10(m),redshift)
          curval = L_CO(sfri,redshift)

       elseif (profile==ithi) then
          curval = MHtoMHI(m,redshift)
          curval = L_HI(curval)
       endif
       curval = I_line(curval,redshift,chi,dnu)
       curval = T_line(curval,nuobsj)
       fsmapl(j+npix*(k-1)) = fsmapl(j+npix*(k-1)) + curval
    else
    
    ! For Pasted Profiles
    if(thetamin <= -fov/2) thetamin = -fov/2 + 1e-10
    if(thetamax >=  fov/2) thetamax =  fov/2 - 1e-10
    
    if(  phimin <= -fov/2) phimin = -fov/2 + 1e-10
    if(  phimax >=  fov/2) phimax =  fov/2 - 1e-10
    
    jmin = max(1,int((thetamin + fov/2)/dpix)+1)
    jmax = min(int((thetamax + fov/2)/dpix)+1,npix)
    
    kmin = max(1,int((phimin + fov/2)/dpix)+1)
    kmax = min(int((phimax + fov/2)/dpix)+1,npix)
    
    if(dtheta_halo > maxtheta) then
       maxtheta = dtheta_halo
       maxthetachih = chi
       maxthetamass = m
    endif
    
    do k=kmin,kmax
       phip   = -fov/2 + (k-0.5)*dpix
       do j=jmin,jmax           
          thetap   = -fov/2 + (j-0.5)*dpix
          
          dphi   = phip   - phi_halo
          dtheta = thetap - theta_halo
          
          angle = sqrt(dphi**2+dtheta**2)
          
          if(angle<dtheta_halo) then
             curval =  mf * v *interpolate_table(angle,redshift,m)              
             fsmapl(j+npix*(k-1)) = fsmapl(j+npix*(k-1)) + curval
if(curval/=0.) write(*,*) 'badger 1 - haloproject.f90, mf,v=', mf, v
          endif
          
       enddo
    enddo
    endif
    return 
    
  end subroutine projecthalo_flatsky

end module haloproject
