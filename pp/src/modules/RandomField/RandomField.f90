module RandomField

  !=======================================================================
  !
  ! GENERATES UNIGRID INITIAL CONDITIONS FOR COSMOLOGICAL SIMULATIONS
  !
  !                                    AUTHOR:             MARCELO ALVAREZ
  !                                 LAST EDIT:  9 May 2024 by Nate Carlson
  !
  !=======================================================================

contains

  subroutine RandomField_seed(seed_in)
    
    ! INCLUDE NECESSARY MODULES
    use gaussian_field

    ! IMPLICIT STATEMENT    
    implicit none

    ! INPUT PARAMETERS
    integer seed_in

    seed      = seed_in

    return

  end subroutine RandomField_seed

  subroutine RandomField_init(filepktab_in,ireadfield_in,outcode_in,incode_in,fielddir_in,&
       boxsize_in,nmesh_in,nbuff_in,ntile_in,ainit_in,NonGauss_in,fNL_in,A_nG_in,B_nG_in,R_nG_in,&
       H_e_in,ng_seed_in,fileTktab_in,h_in)

    !=====================================================================
    ! DECLARATIONS BEGIN  
    
    ! INCLUDE NECESSARY MODULES
    use grid
    use tiles
    use mpivars
    use gaussian_field
    use pktable
    use globalvars
    use fftw_interface

    ! IMPLICIT STATEMENT    
    implicit none

    double precision, parameter :: pi=3.1415926535897932384626433832795_dp

    ! INPUT PARAMETERS
    character(len=512) filepktab_in, fileTktab_in, outcode_in, incode_in, fielddir_in
    real               boxsize_in, ainit_in, fNL_in, A_nG_in, B_nG_in, R_nG_in, H_e_in, h_in
    integer            nmesh_in, nbuff_in, ntile_in, NonGauss_in, ireadfield_in,ng_seed_in

    ! OUTPUT PARAMETERS
    integer total_local_sizes_out
    integer fieldid
    
    ! DECLARATIONS END
    !=====================================================================

    !=====================================================================
    ! EXCUTABLE BEGIN

    filepktab = filepktab_in
    fileTktab = fileTktab_in
    outcode   = outcode_in
    incode    = incode_in
    fielddir  = fielddir_in
    boxsize   = boxsize_in
    nmesh     = nmesh_in
    nbuff     = nbuff_in
    ntile     = ntile_in
    ainit     = ainit_in
    NonGauss  = NonGauss_in
    fNL       = fNL_in
    A_nG      = A_nG_in
    B_nG      = B_nG_in
    R_nG      = R_nG_in
    H_e       = H_e_in
    ng_seed   = ng_seed_in
    h         = h_in
    nsub = nmesh - 2 * nbuff
    n    = nsub * ntile + 2 * nbuff    

    if(ireadfield_in<2) then !if reading in density and displacements (ireadfield==2) do not care about global ffts
       if(mod(n,ntasks).ne.0) then
          if(myid==0) then
             write(*,*) '  ERROR: Number n not evenly divisible by ntasks in FFTW'
             write(*,*) '         exiting ...'
          endif
          call mpi_finalize(ierr)
          stop
       endif
    endif
    z = 1./ainit - 1.

    n12  = n / 2   ! n/2 rounded down to nearest integer
    n21  = n12 + 1 ! this is the length of the last axis in a real-to-complex FFT
    n2p1 = 2 * (n12 + 1) ! this is the length of the last axis in an in-place real-to-complex FFT
    dk   = 2 * pi / boxsize ! Fourier-space line element
    kmax = n * dk           ! max wavenumber k in the box where k\in[0,kmax]
    d3k  = dk**3            ! Fourier-space volume element
    V_k  = d3k*n**3
    ! Note that boxsize is the sidelength in Mpc/h of the simulation volume excluding buffers (i.e. n-2*nmesh), so shouldn't V_k here be d3k*(n-2*nmesh)**3 because the sidelength of an n^3 box is boxsize+2*nbuff*alatt?

    ! INITIALIZE FFTW  
    nfft   = n 
    nfft_s = nmesh
    call fftw_initialize()

    ! READ POWER SPECTRUM, transfer function and chi FROM TABLE
    call read_pktable(A_nG,B_nG,R_nG) ! subroutine in RandomField/pktable.f90 
    pksav       = sqrt(pksav_in    *V_k)
    tsav        = sqrt(tsav_in     *V_k)
    pkchisav    = sqrt(pkchisav_in *V_k)
    pkchixsav   = sqrt(pkchixsav_in*V_k)

    ! Read the \Delta\phi \to \Delta\zeta transfer function, setting array T_dphi2dzeta from
    ! src/modules/RandomField/pktable.f90
    if((NonGauss>=8).and.(NonGauss<=10)) then
       call read_Tktable(B_nG,H_e)
    elseif(NonGauss>=11) then
       call read_Tktable(B_nG,H_e,.True.,.True.)
    endif

    return

  end subroutine RandomField_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine RandomField_make(fieldid,fieldp,fieldpc,fieldpgc)
    use grid
    use mpivars
    use gaussian_field
    use pktable
    use globalvars
    !for spherical profile
    use fftw_interface
    use mpivars
    use openmpvars
    use tiles
    use chi2zeta
    implicit none

    integer fieldid
    !for single spherical overdensity
    real alpha, kR
    real ro, rmax, ar, dr, ff, rx,ry,rz, fcrit
    real*8 lavg,avg,lsig,sig,lsig8,sig8

    real(C_FLOAT),            pointer :: fieldp(:,:,:)
    complex(C_FLOAT_COMPLEX), pointer :: fieldpc(:,:,:)
    complex(C_FLOAT_COMPLEX), pointer :: fieldpgc(:,:,:)

    real sigmachi

    ! For intermittent, spatially-localized non-gaussianity (NonGauss=5)
    integer ng_peak_rate ! expected number of non-guassianity peaks per (Mpc/h)^3
    integer Rpk_seed  ! seed for peak scales

    ! --------------------------------------------------------------------
    ! GENERATE WHITE NOISE, CONVOLVE WITH P(K), AND OUTPUT
    ! --------------------------------------------------------------------

    delta   => fieldp   ! pointer to real field
    deltac  => fieldpc  ! pointer to complex field
    deltagc => fieldpgc ! pointer to global complex field
    
    deltac = deltagc

    if(fieldid>0.or.fieldid==-1) then

       !Assign each delta(i,j,k) gaussian random deviate with sig=1, mu=0
       call generate_noise

       if(NonGauss == 0) then
          code='rho'
          call convolve_noise !returns deltac as complex DFT of delta convolved with power spectrum, delta at this stage is still just gaussian noise

       elseif(NonGauss > 0) then
          !1 - classic fNL with correlated seed
          !2 - classic fNL with uncorrelated seed
          !3 - intermittent NG (Gaussian Spike)
          !4 - intermittent NG (Chaotic Billiards)
          !5 - spatially localized intermittent NG with random peaks
          !6 - spatially localized intermittent NG from chi power spectrum

          if(NonGauss==1) then
             ! Finds zeta = zeta_G + f_NL (zeta_G^2 - <zeta_G^2>) 
             ! Where the gaussian zeta_G is used in the non-gaussian term
             code='zetag'
             call convolve_noise !returns 'delta' and 'deltac' as zeta_lin(x) and zeta_lin(k)
!             delta=delta*V_k
             call RandomField_Output(-15,delta) !outputs zeta_lin to <run_dir>/fields/zetag_...
!             delta=delta/V_k
             !now have zeta_lin field. 
             !zeta_NL = zeta_lin + f_NL * [zeta_lin**2 - <zeta_lin**2>]
             lavg = sum(delta(1:n,:,:)**2) / real(n)**3 !=<zeta_lin^2>
             call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr) !outputs avg as lavg summed over parallel nodes
             delta = delta + fNL * (delta**2 - avg) ! zeta_NL(x)

             ! Output zeta = zeta_g + fNL(zeta_g^2 - <zeta_g^2>)
             code='zetang'
!             delta=delta*V_k
             call RandomField_Output(-17,delta)
!             delta=delta/V_k

             code='zeta2delta'
             call convolve_noise !sets deltac to delta_NL(k)
          
          elseif(NonGauss==2) then
             ! Finds zeta = zeta_G + f_NL (chi^2 - <chi_G^2>)
             ! where chi is a gaussian random field generated using the same powerspectrum P(k) as zeta_G, but with a different initial random seed
             ! should check normalization and replace the delta = delta * 0.9925 with something acting on both Gaussian and non-Gaussian parts
             code='rho'
             call convolve_noise !returns deltac as complex DFT of delta convolved with power spectrum, delta at this stage is still just gaussian noise
             call iFT_field ! returns delta as iDFT of deltac (i.e. delta is now gaussian overdensity)
             !Output field and create uncorrelated intermittent NG field
             delta = delta * 0.9925
             call RandomField_Output(-1,delta)

             code='zetag'
             delta = 0.0
             seed  = ng_seed
             call generate_noise
             call convolve_noise
!             delta=delta*V_k
             call RandomField_Output(-15,delta)
!             delta=delta/V_k
             !zeta_NL = zeta_lin + f_NL * [zeta_lin**2 - <zeta_lin**2>]
             delta = delta**2
             lavg = sum(delta(1:n,:,:)) / real(n)**3 !=<delta^2>
             call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr) !outputs avg as lavg summed over parallel nodes
             !get zeta**2 contribution
             delta = fNL * (delta - avg)
             
             code='zeta2delta'
             call convolve_noise
             call iFT_field
             !add uncorrelated delta field in place
             call RandomField_Input(-16, delta) ! reads in gaussian part and adds it to delta
             call FT_field

          elseif(NonGauss==3.or.NonGauss==4) then
             code='rho'
             call convolve_noise
             call iFT_field
 !            delta = delta !* 0.99
             !now have delta_g field. 
             !Output field and create uncorrelated intermittent NG field
             call RandomField_Output(-1, delta)

             seed  = ng_seed ! using a new seed number gives us an uncorrelated field
             delta = 0.0
             code  = 'chi'
             call generate_noise
             call convolve_noise ! convolve_noise with code 'chi' convolves the white noise field with an independent chi power spectrum that is not related to the CMB power spectrum used above
             if(NonGauss==3) then
                sigmachi = 1.1e-6
                delta    = delta/sigmachi !change units of field to sigmas  
             elseif(NonGauss==4) then
                delta    = delta + 5.8e-8 ! this is some sort of fudge factor to account for normalization of P_chi, I think
             endif
             !now have chi field
             !run through chi to zeta mapping function
             call chi2zetaFNL ! creates an intermittent NG part of the field by mapping zeta to chi (i.e. see figure 1 in PRL 103, 071301 (2009) "Non-Gaussian Curvature Spikes from Chaotic Billiards in Inflation Preheating")
             !code="zetang"
             !call RandomField_Output(delta)

             code='zetang'
!             delta=delta*V_k
             call RandomField_Output(-17,delta) ! outputs intermittent zeta field
!             delta=delta/V_k
             code='zeta2delta'
             call convolve_noise ! gets delta_NL(k)
             call iFT_field      ! gets delta_NL(x)
!             delta=delta*fNL
!             delta = delta/10 
             call RandomField_Input(-16, delta)! reads in gaussian part and adds it to delta
             !delta_NL. change code incase fields needs to be written out
             !delta=(.8111/.8304)*delta
             call FT_field

          elseif(NonGauss==5) then
             ! Primordial intermittent non-Gaussianity (PING) from peak statistics
             ! For spatially localized, intermittent NG, we start with a
             ! field that is mostly gaussian with a few spatially localized
             ! non-gaussianities (indexed c) at positions x_c. The non-
             ! gaussian part of zeta is summed over these non-gaussian
             ! peaks c and composed of 3 factors (where zeta_c is
             ! zeta_NG(x=x_c):
             !   1)  Z_c ~ zeta_c V_c is like a luminosity for NG peak c,
             !       the product of the peak height in zeta and the peaks
             !       approximate volume in x
             !   2)  chi_{zeta_c} has units of 1/volume and characterizes
             !       the shape of zeta at x_c, its form will eventually
             !       come from Tom's lattice simmulations, but in the first
             !       instance we may simply choose a simple form to see
             !       what it looks like.
             !   3)  n_c(x) looks like a numberdensity, probably with
             !       n_c(x_c)=1 and n_c(x_n)=0 for any x_n that is not near
             !       any of the peaks x_c.
             ! This means that zeta ~ zeta_G except at the peaks, hence
             ! "spatially localized". The form of zeta is
             ! zeta(x) = zeta_G(x) + sum_c Z_c chi_{zeta_c}(x-x_c) n_c(x)

             ! Convolve white-noise Gaussian random field with matter power
             ! spectrum and zeta transfer function to get Gaussian
             ! component zeta_G(k) of zeta, and write it to a file
             code='zetag'
             call convolve_noise ! sets deltac equal to zeta_G(k)
             call iFT_field      ! sets delta equal to zeta_G(x)
!             delta=delta*V_k
             call RandomField_Output(-15, delta) ! outputs zeta_G(x)
!             delta=delta/V_k

             ! With zeta_G(x) written to a file, we zero the Fourier
             ! conjugate so that it can be populated with intermittent,
             ! spatially-localised non-Gaussian peaks, giving the non-
             ! Gaussian component zeta_NG(x) of zeta
             deltac = deltac*0.0
             call zetang_sli(deltac,ng_peak_rate,ng_seed,Rpk_seed,n,&
                boxsize)    ! sets deltac equal to zeta_NG(k)
             call iFT_field ! sets delta equal to zeta_NG(x)

             ! Reads zeta_G(x) from file "in place", setting the array
             ! variable delta to the sum zeta(x)=zeta_G(x)+zeta_NG(x)
!             delta=delta*V_k
             call RandomField_Input(-15, delta) ! sets delta to zeta(x)
             call Randomfield_Output(-17, delta) ! writes out zeta(x)
!             delta=delta/V_k
             call FT_field                      ! sets deltac to zeta(k)

          elseif(NonGauss==6) then
             ! A toy model of spatially-localized, intermittent NG
             ! from an uncorrelated inflationary chi field with a bump in
             ! it's primordial power spectrum

             ! Set arrays delta and deltac to the Gaussian part of zeta
             ! field, zeta_G(x) and its Fourier transform zeta_G(k)
             code='zetag'        ! Sets delta equal to zeta_G(x)
             call convolve_noise ! and deltac equal to zeta_G(k)

             ! Writes the Gaussian part of zeta, zeta_G(x), to a file in
             ! the fields directory so that we can reuse the array variable
             ! delta in determining the non-Gaussian part of zeta
!             delta=delta*V_k
             call RandomField_Output(-15, delta)
!             delta=delta/V_k

             ! Get the density field corresponding to the Gaussian
             ! component of zeta and output it
             code='zeta2delta'
             call convolve_noise ! sets deltac to rho_G(k)
             call iFT_field ! sets delta to rho_G(x)
             call RandomField_Output(-98, delta) ! writes out delta_G(x)

             ! Zeros the array variable delta, then uses new seed to
             ! populate delta with a Gaussian random deviate, producing a
             ! white noise field that is uncorrelated with that used in
             ! generating zeta_G
             delta = delta*0.0
             seed  = ng_seed
             call generate_noise ! sets delta to uncorrelated white noise

             ! Convolve the white-noise field in delta with the power
             ! spectrum from the transverse inflaton field chi(k), then
             ! transform this to the form of a zeta field to get the non-
             ! Gaussian component of zeta(x), which produces intermittent,
             ! spatially-localised non-Gaussian peaks. Finally, we multiply
             ! by fNL (which, for simplicity is the same variable used in
             ! previous runs, but doesn't actually represnt a true f_NL
             ! which can be constrained through bispectrum measurements).
             code='chi_extreme'
             call convolve_noise ! sets delta to real-space chi(x)
             call RandomField_Output(-97,delta) ! writes out chi(x)
             call chi_to_zeta
             delta=fNL*delta
             ! sets delta to zeta_NG(x)=f_nG [ chi_E^2(x) - <chi_E^2> ]

             ! Read zeta_G(x) from file written out above "in place" so
             ! that the array varaible delta is set to
             !     delta(:,:,:) --> zeta(x)=zeta_G(x)+fNL*zeta_NG(x)
             ! Finally, we convolve with the zeta transfer function to
             ! convert this to a Fourier-space density field
!             delta=delta*V_k
             call RandomField_Input(-15, delta) ! sets delta to zeta(x)
             call RandomField_Output(-17, delta) ! writes zeta(x) to file
!             delta=delta/V_k
             code='zeta2delta'
             call convolve_noise ! sets deltac to rho(k)
             deltac = (.8111/1.5251)*deltac

          elseif(NonGauss==7.or.NonGauss==9) then
             ! Primordial intermittent non-Gaussianity (PING) from a quartic instability in the 
             ! inflationary potential generating a pulse of non-Gaussian structure formation at a
             ! characteristic scale.
             ! Normalised so that sigma_8 = sigma_8 from Planck 2018

             ! Set array delta to the Gaussian part of zeta field zeta_G(x), write zeta_G(x) to file
             code='zetag'
             call convolve_noise                 ! set delta to zeta_G(x)
             call RandomField_Output(-15, delta) ! write zeta_G(x) to file

             ! Set array delta to correpsonding density field rho_G(x), write rho_G(x) to file
             code='zeta2delta'
             call convolve_noise                 ! set deltac to rho_G(k)
             call iFT_field                      ! set delta to rho_G(x)
             call RandomField_Output(-98, delta) ! write rho_G(x) to file 

             ! Smooth field on scale R=8 h^-1 Mpc to measure sigma_8,G for the rho_G(x)
             call FT_field ! set deltac to rho_G(k)
             do k=1,local_nz
                do j=1,n
                   do i=1,n12+1
                      kR = sqrt( float(i-1)**2. + float(j-1)**2. + float(k-1)**2. ) * dk * 8./h
                      if(kR/=0.) deltac(i,j,k) = deltac(i,j,k)*3./kR**3. * (sin(kR)-kR*cos(kR))
                   enddo
                enddo
             enddo
             call iFT_field ! set delta to smoothed field <rho_G,8>(x)

             ! Calculate sigma_8,G for Gaussian part of density <rho_G,8>(x)
             lavg = sum(delta(1:n,:,:)) / real(n)**3
             call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
             lsig = sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
             call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
             if(myid==0) write(*,*) 'sigma8 <delta_G,8> = ', sig

             ! Overwrite delta to white noise uncorrelated with delta_G(x)
             delta = delta*0.0
             seed  = ng_seed ! change this to read in from parameter file
             call generate_noise ! sets delta to uncorrelated white noise

             ! Set array delta to inflationary field chi(x) then perform
             ! non-linear transformation to non-Gaussian part of zeta(x)
             code='chi'
             call convolve_noise                ! set delta chi(x)
             call RandomField_Output(-97,delta) ! write chi(x) to file
             call chi_to_zeta ! NL transformation to set delta to
                               ! zeta_NG(x) = chi_E^2(x) - <chi_E^2>

             ! Read in-place to add zeta_G(x) from file to arrray delta
             !     delta(:,:,:) --> zeta(x) = zeta_G(x) + zeta_NG(x)
             ! convolve with the transfer function to set delta to rho(x)
             call RandomField_Input(-15, delta)  ! set delta to zeta(x)
             call RandomField_Output(-17, delta) ! write zeta(x) to file ! NOTE THAT THIS ZETA IS NOT NORMALIZED TO SIGMA8
             code='zeta2delta'
             call convolve_noise                 ! set deltac to rho(k)

             ! Smooth field on scale R=8 h^-1 Mpc
             do k=1,local_nz
                do j=1,n
                   do i=1,n12+1
                      kR = sqrt( float(i-1)**2. + float(j-1)**2. + &
                           float(k-1)**2. ) * dk * 8./h
                      if(kR/=0.) deltac(i,j,k) = deltac(i,j,k) *3./kR**3. &
                                                 * (sin(kR)-kR*cos(kR))
                   enddo
                enddo
             enddo
             call iFT_field ! sets delta to <rho,8>(x)
             
             ! Calculate sigma8 for full density field <rho,8>(x)
             lavg = sum(delta(1:n,:,:)) / real(n)**3
             call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                                mpi_comm_world,ierr)
             lsig8 = sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
             call mpi_allreduce(lsig8,sig8,1,mpi_double_precision,mpi_sum,&
                                mpi_comm_world,ierr)

             if(myid==0) write(*,*) 'sigma8 <delta_nG,8> = ', sig8
             if(myid==0) write(*,*) '<delta_G,8>/<delta_nG,8> = ', sig/sig8

             ! Normalize to sigma8, overwrite zeta(x) file with normalized
             ! array, and convolve to get normalized rho(k)
             call RandomField_Input(-17, delta) ! sets delta to zeta(x)
             delta=delta*sig/sig8               ! normalize zeta(x)
             call RandomField_Output(-17,delta) ! overwrite zeta(x) file
             code='zeta2delta'
             call convolve_noise                ! set deltac to rho(k)

          elseif(NonGauss==8.or.NonGauss==10.or.NonGauss==11) then
             ! Primordial intermittent non-Gaussianity (PING) from a quartic instability in the 
             ! inflationary potential generating a pulse of non-Gaussian structure formation at a
             ! characteristic scale.

             ! Set array delta to the Gaussian part of zeta field zeta_G(x), write zeta_G(x) to file
             code='zetag'
             call convolve_noise                 ! set delta to zeta_G(x)
             call RandomField_Output(-15, delta) ! write zeta_G(x) to file

             ! Set array delta to correpsonding density field rho_G(x), write rho_G(x) to file
             code='zeta2delta'
             call convolve_noise                 ! set deltac to rho_G(k)
             call iFT_field                      ! set delta to rho_G(x)
             call RandomField_Output(-98, delta) ! write rho_G(x) to file 

             ! Overwrite delta to white noise uncorrelated with delta_G(x)
             delta = delta*0.0
             seed  = ng_seed ! change this to read in from parameter file
             call generate_noise ! sets delta to uncorrelated white noise

             ! Set array delta to inflationary field chi(x) then perform
             ! non-linear transformation to non-Gaussian part of zeta(x)
             code='chi'
             call convolve_noise                ! set delta chi(x)
             call RandomField_Output(-97,delta) ! write chi(x) to file
             call chi_to_zeta ! NL transformation to set delta to
                              ! zeta_NG(x) = chi_E^2(x) - <chi_E^2>

             ! Read in-place to add zeta_G(x) from file to arrray delta
             !     delta(:,:,:) --> zeta(x) = zeta_G(x) + zeta_NG(x)
             ! convolve with the transfer function to set delta to rho(x)
             call RandomField_Input(-15, delta) ! set delta to zeta(x)
             call RandomField_Output(-17,delta) ! write zeta(x) to file
             code='zeta2delta'
             call convolve_noise                ! set deltac to rho(k)
          endif
       endif
    endif

    if(fieldid==-2)  code='xlpt1'
    if(fieldid==-3)  code='ylpt1'
    if(fieldid==-4)  code='zlpt1'
    if(fieldid==-5)  code='phi11'
    if(fieldid==-6)  code='phi22'
    if(fieldid==-7)  code='phi33'
    if(fieldid==-8)  code='phi21'
    if(fieldid==-9)  code='phi31'
    if(fieldid==-10) code='phi32'
    if(fieldid==-11) code='xlpt2'
    if(fieldid==-12) code='ylpt2'
    if(fieldid==-13) code='zlpt2'
    if(fieldid==-14) code='lapld'
    
    if(fieldid<-1) call convolve_noise

    !testing single spherical overdensity
    if(fieldid==-15) then
       !set parameters for profile
       ro = 35.
       rmax = 100. 
       alpha = 1.5
       fcrit = 1.686
       dr = boxsize/n
       ff = 4*3.14159/3 * (rmax/(dr*n))**3
       do i=1,n
          rx = (i-n/2) * dr
          do j=1,n
             ry = (j-n/2) * dr
             do k=1,local_nz!local_nz assigned in fftw_initialize() as mpi slab size, equal to n in runs with ntile=1
                rz = (k+local_z_start - n/2)*dr
                ar=sqrt(rx**2+ry**2+rz**2+dr)
                if(ar < rmax) then
                   delta(i,j,k) = (3-alpha)/3 * fcrit * (ar/ro)**(-alpha)
                else
                   delta(i,j,k) = ff/(ff-1) * fcrit * (rmax/ro)**(-alpha) 
                endif
             enddo
          enddo
       enddo

       lavg = sum(delta(1:n,:,:)) / real(n)**3
       call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
            mpi_comm_world,ierr)
       lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
       call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
            mpi_comm_world,ierr)
      
       if(myid==0) write(*,*) "avg, sigma of spherical profile = ", avg,sqrt(sig)
          
       call mpi_barrier(mpi_comm_world,ierr)

    endif

    return

  end subroutine RandomField_make
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Non-Gaussian component zeta_NG(k) of zeta field with spatially-
  ! localised, intermittent non-Gaussian peaks
  subroutine zetang_sli( zetak_ng, ng_peak_rate, ng_seed, Rpk_seed, n, s )
    use, intrinsic :: iso_c_binding
    use random
    use mpivars
    implicit none

    ! Input variables
    complex(c_float_complex), pointer :: zetak_ng(:,:,:) ! zeta_NG(k)
    integer ng_peak_rate ! expected number of non-guassianity peaks per Mpc/h
    integer ng_seed   ! seed for non-gaussian peak generation
    integer Rpk_seed  ! seed for peak scales
    integer(kind=8) n ! sidelength of simulation box in cells (w/ buffers)
    real s

    ! Random number generation
    integer ng_seeds(MaxRandNumStreams,IRandNumSize),&
            Rpk_seeds(MaxRandNumStreams,IRandNumSize)
    real(kind=8) ng_dev ! uniform random deviate

    ! Variables for each localised peak
    integer q
    real cutoff ! 
    real k_1,k_2,k_3  ! wavevector
    real(kind=8),allocatable,dimension(:) :: Z0_c ! height of ng peak c
    real,allocatable,dimension(:) :: R_c          ! scale of ng peak c
    real,allocatable,dimension(:) :: x_c,y_c,z_c  ! array of peak mean values for each ng peak c
    integer ng_peak_positions(n,n,n)
    integer N_ng_peaks

    ! We start with a Poisson process, populating the field with a uniform
    ! random deviate, and taking points in the field for which this deviate
    ! is greater than a cutoff to be the sites of the non-Gaussian peaks. 
    ! This cutoff is a function of expected rate of peaks per unit volume
    cutoff = 1-ng_peak_rate*s**3*n**-3

    ! Initialize random number seeds
    call rans(ntasks,ng_seed ,ng_seeds )
    call rans(ntasks,Rpk_seed,Rpk_seeds)

    ! Find locations of non-Gaussianity peaks
    do k_1=1,n
      do k_2=1,n
        do k_3=1,n
          call ranf(ng_seeds(1,:),ng_dev)
          if (ng_dev>=cutoff) then
            ng_peak_positions(k_1,k_2,k_3)=1
          else
            ng_peak_positions(k_1,k_2,k_3)=0
          endif
        enddo
      enddo
    enddo
    N_ng_peaks = sum(ng_peak_positions) ! number of non-gaussianity peaks

    ! Allocate arrays describing:
    allocate(Z0_c(N_ng_peaks)) ! relative prominence of non-G peaks
    allocate(x_c(N_ng_peaks))  ! 
    allocate(y_c(N_ng_peaks))  ! location of non-G peak centres
    allocate(z_c(N_ng_peaks))  ! 

    ! Get peak positions from ng_peak_positions(:,:,:) and use second
    ! deviate to determine Z_0 for each peak
    q=0 ! use integer q as counter
    do k_1=1,n
      do k_2=1,n
        do k_3=1,n
          if (ng_peak_positions(k_1,k_2,k_3)==1) then
            x_c(q)=k_1*s/n
            y_c(q)=k_2*s/n
            z_c(q)=k_3*s/n
            call ranf(Rpk_seeds(1,:),Z0_c(q))
            q=q+1
          endif
        enddo
      enddo
    enddo

!    ! Next, in parallel, we create the non-Gaussian part of zeta
!!$OMP PARALLEL DO &
!!$OMP DEFAULT(FIRSTPRIVATE) &
!!$OMP SCHEDULE(STATIC) &
!!$OMP SHARED(deltac) 
!    do k=1,local_nz!local_nz assigned in fftw_initialize() as mpi slab size, equal to     n in runs with ntile=1
!      kz=(k+local_z_start-1)*dk ! local_z_start assigned in fftw_initialize() in modu    le fftw_interface.f90,
!      if (k+local_z_start.gt.n12) kz=kz-kmax        ! index of slab start, equal to 0     for runs with ntile=1
!
!      do j=1,n
!        ky=(j-1)*dk
!        if (j.gt.n12) ky=ky-kmax
!
!        do i=1,n12+1
!          kx=(i-1)*dk
!          if (i.gt.n12) kx=kx-kmax
!
!          ak=sqrt(kx**2+ky**2+kz**2)
!
!        enddo
!      enddo
!    enddo

  end subroutine zetang_sli

!  ! Fourier transform of normalised 3D Gaussian
!  real function gaussk(k_1,k_2,k_3,x_0,y_0,z_0,R_0)
!    implicit none
!    real k_1,k_2,k_3 ! wavevector
!    real x_0,y_0,z_0 ! mean of real-space Gaussian
!    real R_0         ! standard deviation of real-space Gaussian
!    gaussk = gaussk_j(k_1,x_0,R_0) &
!           * gaussk_j(k_2,y_0,R_0) &
!           * gaussk_j(k_3,z_0,R_0)
!    return
!  end function gaussk
!   
!  ! Fourier transfrom of normalised 1D Gaussian
!  real function gaussk_j(k_j,x_0j,R_0)
!    implicit none
!    double precision, parameter :: pi=3.1415926535897932384626433832795_dp
!    real k_j  ! wavenumber
!    real x_0j ! mean of real-space Gaussian
!    real R_0  ! standard deviation of real-space Gaussian
!    gaussk_j = x_0j/(2*pi*R_0) * cmplx( cos(k_j*x_0j) , sin(k_j*x_0j) ) & 
!      * exp( -x_j**2*x_0j**2/2 )
!    return  
!  end function gaussk_j
!   
!  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real function T_interp(k_in,T_in,ii,jj,kk)
    use gaussian_field ! Needed to access dk, the fourier space voxel grid size
    implicit none

    ! Input values
    real, intent(in) :: k_in(:), T_in(:)
    integer             ii, jj, kk

    ! Variables
    integer len_k, i_floor, i_ceil
    real    i_interp, k_vox

    ! Length of input arrays
    len_k = size(k_in)

    ! Determine the wavenumber of the Fourier space voxel (i,j,k)
    k_vox = sqrt( float(ii-1)**2. + float(jj-1)**2. + float(kk-1)**2. )*dk

    ! Given spacing of k table, determine effective index in list k (as a floating point number) for
    ! Fourier space voxel (ii,jj,kk) with wavenumber k_vox
    if(k_vox==0.0) then
       i_interp = 1
    else
       i_interp = 1 + (len_k-1) * log10(k_vox/k_in(1)) / log10(k_in(len_k)/k_in(1))
    endif
    
    ! Find the integer index such that k(i_floor) <= k_vox < k(i_floor+1)
    i_floor = floor(i_interp)
    i_ceil  = i_floor+1

    ! Extrpolate as power law if k_vox is lower than any value in the table k
    if( i_floor < 1 ) then
       T_interp = T_in(1) * 10**(   log10( T_in(2) / T_in(1) ) &
                                  / log10( k_in(2) / k_in(1) ) &
                                  * log10( k_vox   / k_in(1) ) )

    ! Extrapolate as power law if k_vox is higher than any balue in the table k_tab
    elseif( i_floor+1 > len_k ) then
       T_interp = T_in(len_k) * 10**(   log10( T_in(len_k) / T_in(len_k-1) ) &
                                      / log10( k_in(len_k) / k_in(len_k-1) ) &
                                      * log10( k_vox       / k_in(len_k  ) ) )

    ! Interpolate if k_vox within range spanned by table k_tab (linear in logT v logk)
    else
       if( (T_in(i_floor)<=0.0) .or. (T_in(i_floor+1)<=0.0) .or. (k_vox==0.0) ) then
          T_interp = 0.0 ! Set to zero if negative values or k_vox = 0
       else
          T_interp = T_in(i_floor) * 10**(   log10( T_in(i_ceil) / T_in(i_floor) ) &
                                           / log10( k_in(i_ceil) / k_in(i_floor) ) &
                                           * log10( k_vox        / k_in(i_floor) ) )
       endif
    endif

    return  
  end function T_interp

  subroutine chi_to_zeta
    ! Used in spatially-localized, intermittent non-gaussianity, where
    ! zeta(x)=zeta_G(x) + f_nG * zeta_NG(x)
    ! and
    ! zeta_NG(x) = m ( A chi^2 + B chi + C ) + b
    ! where B~C~b~0
    use grid
    use mpivars
    use timing_diagnostics
    use fftw_interface
    use globalvars !uncommented 18 Nov 2022
    use pktable
    use gaussian_field ! Needed to access dk, the fourier space voxel grid size
    implicit none
    real(kind=8) lavg,avg,lsig,sig
    real slope, i_interp, Tk_interp, Wchi_interp
    integer len_Tk, i_low
    real, allocatable :: k_tab(:), Tk_tab(:), Tchi2phi_tab(:)

    ! -------------------------------------------------------------------------------------------- !
    ! The first step is to convert from \chi to \Delta\phi                                         !
    ! -------------------------------------------------------------------------------------------- !

    ! In non-Gaussianity models 7-10, \Delta\phi \propto \chi^2
    if((NonGauss>=7).and.(NonGauss<=10)) then

       ! The ratio \Delta\zeta/\Delta\phi is proportional to \lambda_\chi,
       ! Note that in non-Gaussian model 7, R_nG = \lambda_\chi
       slope=(3.84158538e+2)*R_nG**(-4.39886747e-1)!+(-2.40220085)
       ! Fit to a few samples of the transfer function, requires more rigorous
       ! treatment.  -- Nate Carlson, January 2023

       if(myid==0) call timer_begin
 
       ! Assuming delta(:,:,:) is the chi_E(x) field, determines the
       ! zeta_nG(x) field element-by-element and computes its average
       lavg=0.
       do k=1,local_nz ! Loop over the local mpiFFTW slab (n x n x local_nz
          do j=1,n     ! slice of the full n x n x n field)
             do i=1,n  !
                delta(i,j,k) = chi2Deltaphi(delta(i,j,k),A_nG,0.,0.)

                ! For the simplified model, approximate deltaphi2zeta with a
                ! functional
                if((NonGauss==7).or.(NonGauss==8)) then
                    delta(i,j,k) = Deltaphi2zeta(delta(i,j,k),slope,0.)
                    lavg         = lavg + delta(i,j,k)
                endif
             enddo
          enddo
       enddo

    ! In non-Gaussianity model 11, \Delta\phi is related to \chi by a transfer function
    elseif(NonGauss==11) then

       ! Square the smoothed chi field
       delta = delta**2

    endif

    ! Apply \Delta\phi \to \Delta\zeta transfer function for NonGaussianity models 9, 10 and 11
    if((NonGauss>=9).and.(NonGauss<=11)) then

       ! Set up wavenumber and \Delta\phi(H_e,x) to \Delta\zeta(H_f,x) transfer function arrays
       len_Tk = size( T_dphi2dzeta(:,1) )
       allocate( k_tab(len_Tk))
       allocate(Tk_tab(len_Tk))
       k_tab  = T_dphi2dzeta(:,1)
       Tk_tab = T_dphi2dzeta(:,2)

       ! set deltac to \delta\phi(k)
       call FT_field

       ! In model 11, we also have a transfer function from \chi^2(H_e,x) to \Delta\phi(H_e,x)
       if(NonGauss==11) then

          ! Set up (W*\chi)^2(H_e,x) to \Delta\phi(H_e,x) transfer function
          allocate(Tchi2phi_tab(len_Tk))
          Tchi2phi_tab = T_dphi2dzeta(:,3)

          ! Loop over all values of deltac, interpolte the (W*\chi)^2(H_e,x) to \Delta\phi(H_e,x)
          ! transfer function
          do k=1,local_nz
             do j=1,n
                do i=1,n12+1
                   Tk_interp     = T_interp(k_tab,Tchi2phi_tab,i,j,k)
                   deltac(i,j,k) = deltac(i,j,k) * Tk_interp
                enddo
             enddo
          enddo
       endif

       ! Loop over all values of deltac, interpolte the \Delta\phi to \Delta\zeta transfer function
       ! from the tablulated values Tk_tab and convolve deltac with transfer function
       do k=1,local_nz
          do j=1,n
             do i=1,n12+1
                Tk_interp     = T_interp(k_tab,Tk_tab,i,j,k)
                deltac(i,j,k) = deltac(i,j,k) * Tk_interp
             enddo
          enddo
       enddo
       call iFT_field ! set delta to \Delta\zeta(x), the non-Gaussian perturbation to \zeta

       ! Calculate the average value of \zeta(x)
       do k=1,local_nz ! Loop over the local mpiFFTW slab (n x n x local_nz
          do j=1,n     ! slice of the full n x n x n field)
             do i=1,n  !
                lavg = lavg + delta(i,j,k)
             enddo
          enddo
       enddo
       
    endif

    ! Divide by n^3 to get slab averages
    lavg=lavg/real(n)**3

    ! Sum slab averages across MPI tasks
    call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
    if(myid==0) write(*,*) "Average of $\\delta\\zeta$ is: ",avg

    ! Enforce zeta_nG(x) has 0 mean
    delta = delta - avg
    if(myid==0) then
       timing_diagnostics_code='chi to zeta to delta_ng'
       call timer_end
    endif

  end subroutine chi_to_zeta

  real function chi2Deltaphi(chi,a,b,c)
    implicit none
    real chi,a,b,c
    chi2Deltaphi=a*chi**2.+b*chi+c
    return
  end function chi2Deltaphi

  real function Deltaphi2zeta(Deltaphi,a,b)
    implicit none
    real Deltaphi,a,b
    Deltaphi2zeta=a*Deltaphi+b
    return
  end function Deltaphi2zeta

  subroutine RandomField_FT(fieldp,fieldpc,ifft)

    use grid
    use mpivars
    use gaussian_field
    use globalvars
    use fftw_interface
    use mpivars
    use openmpvars
    use tiles

    implicit none
    
    real(C_FLOAT),            pointer :: fieldp(:,:,:)
    complex(C_FLOAT_COMPLEX), pointer :: fieldpc(:,:,:)

    integer ifft

    delta   => fieldp
    deltac  => fieldpc

    if(ifft==1) call FT_field   ! subroutine in src/modules/RandomField/gaussian_field.f90
    if(ifft==-1) call iFT_field ! subroutine in src/modules/RandomField/gaussian_field.f90

    return

  end subroutine RandomField_FT

  subroutine RandomField_Output(fieldid, fieldp)

    use intreal_types
    use cosmology
    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use timing_diagnostics
    use globalvars

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    integer ifile
    character *512 outputfile

    real a(1)

    integer iorig
    
    real r4arg
    integer i4arg, indlnb, command_argument_count
    integer fieldid

    real(C_FLOAT), pointer :: fieldp(:,:,:)
    delta  => fieldp

    if(myid==0) call timer_begin

    if(fieldid==-1)  code='rho'
    if(fieldid==-2)  code='xlpt1'
    if(fieldid==-3)  code='ylpt1'
    if(fieldid==-4)  code='zlpt1'
    if(fieldid==-11) code='xlpt2'
    if(fieldid==-12) code='ylpt2'
    if(fieldid==-13) code='zlpt2'
    if(fieldid==-14) code='lapld'
    if(fieldid==-15) code='zetag'
    if(fieldid==-17) code='zetang'
    if(fieldid==-97) code='chi'
    if(fieldid==-98) code='rhog'
    if(fieldid==-99) code='rho_smooth'

    if(code=='rho') then
       outputfile=trim(fielddir)//'Fvec_'//trim(outcode)
    elseif(code=='rhog') then
       outputfile=trim(fielddir)//'rhog_'//trim(outcode)
    elseif(code=='rho_smooth') then
       outputfile=trim(fielddir)//'Fvec_smooth_'//trim(outcode)
    elseif((code=='rho_fNL').or.(code=='lapld')) then
       outputfile=trim(fielddir)//'Fvec_fNL_'//trim(outcode) ! you should probably change this name to like laplacian or somehting because Fvec_fNL sounds like it's a delta
    elseif(code=='zetang') then
       outputfile=trim(fielddir)//'zetang_'//trim(outcode)
    elseif(code=='zetag') then
       outputfile=trim(fielddir)//'zetag_'//trim(outcode)
    elseif(code=='chi') then
       outputfile=trim(fielddir)//'chi_'//trim(outcode)
    elseif(code=='xlpt1') then
       outputfile=trim(fielddir)//'etax1_'//trim(outcode)
    elseif(code=='ylpt1') then
       outputfile=trim(fielddir)//'etay1_'//trim(outcode)
    elseif(code=='zlpt1') then
       outputfile=trim(fielddir)//'etaz1_'//trim(outcode)
    elseif(code=='xlpt2') then
       outputfile=trim(fielddir)//'etax2_'//trim(outcode)
    elseif(code=='ylpt2') then
       outputfile=trim(fielddir)//'etay2_'//trim(outcode)
    elseif(code=='zlpt2') then
       outputfile=trim(fielddir)//'etaz2_'//trim(outcode)
    else
       write(*,*) 'error in code'
       write(*,*) code
       call mpi_finalize(ierr)
       stop
    endif

    offset_tile = n*n*local_z_start
    offset_bytes = offset_tile*int(4,8)+1

    open(unit=33,file=outputfile,access='stream')
    write(33,pos=offset_bytes) (((delta(i,j,k),i=1,n),j=1,n),k=1,local_nz)
    close(33)

    if(myid==0) then
       timing_diagnostics_code = 'writing data'
       call timer_end
    endif
    
    return
    
  end subroutine RandomField_Output

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine RandomField_input(fieldid, fieldp)

    use intreal_types
    use cosmology
    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use timing_diagnostics
    use globalvars

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    integer ifile
    character *512 inputfile

    real a(1), deltaini

    integer iorig
    
    real r4arg
    integer i4arg, indlnb, command_argument_count
    integer fieldid

    double precision lavg,avg,lsig,sig

    real(C_FLOAT), pointer :: fieldp(:,:,:)
    delta  => fieldp

    if(myid==0) call timer_begin

    if(fieldid==-1)  code='rho'
    if(fieldid==-2)  code='xlpt1'
    if(fieldid==-3)  code='ylpt1'
    if(fieldid==-4)  code='zlpt1'
    if(fieldid==-5)  code='phi11'
    if(fieldid==-6)  code='phi22'
    if(fieldid==-7)  code='phi33'
    if(fieldid==-8)  code='phi21'
    if(fieldid==-9)  code='phi31'
    if(fieldid==-10) code='phi32'
    if(fieldid==-11) code='xlpt2'
    if(fieldid==-12) code='ylpt2'
    if(fieldid==-13) code='zlpt2'
    if(fieldid==-14) code='lapld'
    if(fieldid==-15) code='zetag'
    if(fieldid==-16) code='rho_inplace'
    if(fieldid==-17) code='zetang'

    if((code=='rho').or.(code=='rho_inplace')) then
       inputfile=trim(fielddir)//'Fvec_'//trim(incode) 
    elseif((code=='rho_fNL').or.(code=='lapld')) then
       inputfile=trim(fielddir)//'Fvec_fNL_'//trim(incode)
    elseif(code=='zetang') then
       inputfile=trim(fielddir)//'zetang_'//trim(incode)
    elseif((code=='zetag').or.(code=='zetag_inplace')) then
       inputfile=trim(fielddir)//'zetag_'//trim(incode)
    elseif(code=='xlpt1') then
       inputfile=trim(fielddir)//'etax1_'//trim(incode)
    elseif(code=='ylpt1') then
       inputfile=trim(fielddir)//'etay1_'//trim(incode)
    elseif(code=='zlpt1') then
       inputfile=trim(fielddir)//'etaz1_'//trim(incode)
    elseif(code=='xlpt2') then
       inputfile=trim(fielddir)//'etax2_'//trim(incode)
    elseif(code=='ylpt2') then
       inputfile=trim(fielddir)//'etay2_'//trim(incode)
    elseif(code=='zlpt2') then
       inputfile=trim(fielddir)//'etaz2_'//trim(incode)
    else
       write(*,*) 'error in code'
       call mpi_finalize(ierr)
       stop
    endif

    !Find offsets
    offset_tile = n*n*local_z_start
    offset_bytes = offset_tile*int(4,8)+1

    open(unit=33,file=inputfile,access='stream')
    if((code=='zetag').or.(code=='rho_inplace')) then
       !Read in data in place
       offset_i = 0
       do k=1,local_nz
          do j=1,n
             do i=1,n
                read(33,pos=offset_bytes+offset_i) deltaini             
                delta(i,j,k) = delta(i,j,k)+deltaini
                offset_i     = offset_i+4 !4 bytes
             enddo
          enddo
       enddo

    else
       !Read in data
       read(33,pos=offset_bytes) (((delta(i,j,k),i=1,n),j=1,n),k=1,local_nz)

    endif
    close(33)
    
    !Calculate mean and std
    lavg = sum(delta(1:n,:,:)) / real(n)**3
    call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)
    lsig=sum( (delta(1:n,:,:)-avg)**2 ) / real(n)**3
    call mpi_allreduce(lsig,sig,1,mpi_double_precision,mpi_sum,&
                       mpi_comm_world,ierr)

    !Print mean and std
    if(myid==0.and.rf_report==1) then
       if(code.eq.'rho')         write(*,96) avg,sqrt(sig)
       if(code.eq.'xlpt1')       write(*,97) avg,sqrt(sig)
       if(code.eq.'ylpt1')       write(*,98) avg,sqrt(sig)
       if(code.eq.'zlpt1')       write(*,99) avg,sqrt(sig)
       if(code.eq.'zetag')       write(*,100) avg,sqrt(sig)
       if(code.eq.'zeta2delta')  write(*,101) avg,sqrt(sig)
       if(code.eq.'chi')         write(*,102) avg,sqrt(sig)
       if(code.eq.'lapld')       write(*,103) avg,sqrt(sig)
       if(code.eq.'phi11')       write(*,104)
       if(code.eq.'phi22')       write(*,105)
       if(code.eq.'phi33')       write(*,106)
       if(code.eq.'phi21')       write(*,107)
       if(code.eq.'phi31')       write(*,108)
       if(code.eq.'phi32')       write(*,109)
       if(code.eq.'xlpt2')       write(*,110) avg,sqrt(sig)
       if(code.eq.'ylpt2')       write(*,111) avg,sqrt(sig)
       if(code.eq.'zlpt2')       write(*,112) avg,sqrt(sig)
    endif

96   format(3x,'density:',       /5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
97   format(3x,'x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
98   format(3x,'y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
99   format(3x,'z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
100  format(3x,'zeta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
101  format(3x,'zeta2delta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
102  format(3x,'chi:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
103  format(3x,'density laplacian:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
104  format(3x,'phi,11:')
105  format(3x,'phi,22:')
106  format(3x,'phi,33:')
107  format(3x,'phi,21:')
108  format(3x,'phi,31:')
109  format(3x,'phi,32:')
110  format(3x,'2lpt x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
111  format(3x,'2lpt y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
112  format(3x,'2lpt z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)

    if(myid==0) then
       timing_diagnostics_code = 'reading data'
       call timer_end
    endif
    
    return
    
  end subroutine RandomField_input

  subroutine RandomField_input_pointer(fieldid, fieldp, input_field)

    use intreal_types
    use cosmology
    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use timing_diagnostics
    use globalvars

    implicit none

    integer, intent(in) :: fieldid
    real(C_FLOAT), pointer, intent(inout) :: fieldp(:,:,:)
    real(C_FLOAT), pointer, intent(in) :: input_field(:,:,:)

    double precision :: lavg, avg, lsig, sig

    if(myid==0) call timer_begin

    if(fieldid==-1)  code='rho'
    if(fieldid==-2)  code='xlpt1'
    if(fieldid==-3)  code='ylpt1'
    if(fieldid==-4)  code='zlpt1'
    if(fieldid==-5)  code='phi11'
    if(fieldid==-6)  code='phi22'
    if(fieldid==-7)  code='phi33'
    if(fieldid==-8)  code='phi21'
    if(fieldid==-9)  code='phi31'
    if(fieldid==-10) code='phi32'
    if(fieldid==-11) code='xlpt2'
    if(fieldid==-12) code='ylpt2'
    if(fieldid==-13) code='zlpt2'
    if(fieldid==-14) code='lapld'
    if(fieldid==-15) code='zetag'
    if(fieldid==-16) code='rho_inplace'
    if(fieldid==-17) code='zetang'

    ! Copy data from input_field to fieldp
    fieldp = input_field

    ! Calculate mean and std
    lavg = sum(fieldp) / real(n)**3
    call mpi_allreduce(lavg, avg, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    lsig = sum((fieldp - avg)**2) / real(n)**3
    call mpi_allreduce(lsig, sig, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

    ! Print mean and std
    if (myid==0 .and. rf_report==1) then
        select case (code)
            case ('rho')
                write(*,96) avg, sqrt(sig)
            case ('xlpt1')
                write(*,97) avg, sqrt(sig)
            case ('ylpt1')
                write(*,98) avg, sqrt(sig)
            case ('zlpt1')
                write(*,99) avg, sqrt(sig)
            case ('zetag')
                write(*,100) avg, sqrt(sig)
            case ('zeta2delta')
                write(*,101) avg, sqrt(sig)
            case ('chi')
                write(*,102) avg, sqrt(sig)
            case ('lapld')
                write(*,103) avg, sqrt(sig)
            case ('phi11')
                write(*,104)
            case ('phi22')
                write(*,105)
            case ('phi33')
                write(*,106)
            case ('phi21')
                write(*,107)
            case ('phi31')
                write(*,108)
            case ('phi32')
                write(*,109)
            case ('xlpt2')
                write(*,110) avg, sqrt(sig)
            case ('ylpt2')
                write(*,111) avg, sqrt(sig)
            case ('zlpt2')
                write(*,112) avg, sqrt(sig)
        end select
    endif

96  format(3x,'density:',       /5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
97  format(3x,'x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
98  format(3x,'y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
99  format(3x,'z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
100 format(3x,'zeta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
101 format(3x,'zeta2delta:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
102 format(3x,'chi:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
103 format(3x,'density laplacian:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
104 format(3x,'phi,11:')
105 format(3x,'phi,22:')
106 format(3x,'phi,33:')
107 format(3x,'phi,21:')
108 format(3x,'phi,31:')
109 format(3x,'phi,32:')
110 format(3x,'2lpt x-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
111 format(3x,'2lpt y-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)
112 format(3x,'2lpt z-displacement:',/5x' mean = ',1pe13.6,/5x,'sigma = ',1pe13.6)

    if(myid==0) then
        timing_diagnostics_code = 'reading data'
        call timer_end
    endif
    
  end subroutine RandomField_input_pointer
  
end module RandomField


