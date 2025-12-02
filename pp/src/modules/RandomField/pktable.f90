module pktable

  use globalvars
  implicit none
  character(len=512) filepktab, fileTktab, filechi2phitab, file_inst_transfer_tab
  real pktabmin,pktabmax,pktabminl,pktabmaxl,dpktab,dtab
  real, allocatable:: tsav(:), tsav_in(:)
  real, allocatable:: pksav(:), pksav_in(:)
  real, allocatable:: pkchisav(:), pkchisav_in(:)
  real, allocatable:: pkchixsav(:), pkchixsav_in(:)
  real, allocatable:: T_chi2phi(:), T_chi2phi_in(:)
  real, allocatable:: T_dphi2dzeta(:,:)
  integer npktab,itab,nTktab
  integer:: ios

contains

  subroutine read_pktable(A_in,B_in,R_in)
    ! This subroutine reads a file containing an interpollation table for a series of power spectra
    ! and transfer functions as a function of the wavenumber k. These power spectra include the
    ! matter power spectrum, and the transfer function from the primordial scalar perturbation feild
    ! zeta to the matter power spectrum, and in some models the power spectrum for an additional
    ! early-universe field chi and its associtated trasfer functions relating it to zeta.
    implicit none
    integer i, nchi, i_low, offset_bytes, ierr
    real dummy_real, k_read, pk_read, i_interp

    ! Non-Gaussian parameters
    real, optional :: A_in,B_in,R_in
    real, allocatable :: k_extreme(:), p_extreme(:)
    real A_nonG, B_nonG, R_nonG, kchimin, kchimax, P0_low, n_low, P0_high, n_high

    if(present(A_in))then
        A_nonG=A_in
    else
        A_nonG=1.0e-23
    endif
    if(present(B_in))then
        B_nonG=B_in
    else
        B_nonG=1.0e-10
    endif
    if(present(R_in))then
        R_nonG=R_in
    else
        R_nonG=64.0
    endif

    ! Open file containing interpollation tables for power spectrum as a function of wavenumber 
    open(unit=1, file=filepktab, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
        print *, 'Error opening power spectrum interpollation table file: ', trim(filepktab)
        stop
    endif

    ! Read the first wavelength value from the table
    read(unit=1, fmt=*, iostat=ierr) pktabmin
    if (ierr /= 0) then
        print *, 'The power spectrum interpollation table, ', trim(filepktab), ' contains no data.'
        stop
    endif

    ! Determine the length of the power spectrum table and read in its minimum and maximum values
    npktab = 1
    pktabmax = pktabmin
    do
        ! First read
        read(unit=1, fmt=*, iostat=ierr) dummy_real

        ! End do loop if every row has been read
        if (ierr /= 0) exit

        ! increase the count of rows by 1
        npktab = npktab + 1
        
        ! set the maximum value of k
        pktabmax = dummy_real

    enddo

    ! Since k is in log scale, k[j]-k[j-1] is not a constant, rather log_10(k[j]) - log_10(k[j-1]) is constant.
    pktabminl=log10(pktabmin)

    ! dktab represents this value
    pktabmaxl=log10(pktabmax)

    ! In general dktab = dlog_10(k)/dj (where dj is the difference from one array element j to the next j+1)
    dpktab=(pktabmaxl-pktabminl)/(npktab-1)

    close(1)

    ! Allocate power spectrum and transfer function arrays used in core Peak Patch calculation
    if(.not. allocated(tsav))     allocate(tsav(npktab))
    if(.not. allocated(tsav_in))  allocate(tsav_in(npktab)) 
    if(.not. allocated(pksav))    allocate(pksav(npktab))
    if(.not. allocated(pksav_in)) allocate(pksav_in(npktab))

    ! Allocate chi power spectrum arrays used in chaotic billiards preheating non-Gaussianity model
    if(.not. allocated(pkchisav))    allocate(pkchisav(npktab))
    if(.not. allocated(pkchisav_in)) allocate(pkchisav_in(npktab))

    ! Allocate chi power spectrum and transfer function arrays used in PING model
    if(.not. allocated(pkchixsav))    allocate(pkchixsav(npktab))
    if(.not. allocated(pkchixsav_in)) allocate(pkchixsav_in(npktab))
    if(.not. allocated(T_chi2phi))    allocate(T_chi2phi(npktab))
    if(.not. allocated(T_chi2phi_in)) allocate(T_chi2phi_in(npktab))

    call read_table(filepktab, A_nonG, B_nonG)

  end subroutine read_pktable

  subroutine read_Tktable( alpha_e_in , H_e_in , Tchi2phi_flag_in , Tphi2zeta_flag_in )
    ! Reads \Delta\phi \to \Delta\zeta transfer function used in NonGaussianity models >= 9
    ! Reads in W*\chi^2 \to \Delta\phi transfer function used in NonGaussianity models >= 11
    use, intrinsic :: iso_c_binding
    implicit none
    integer, parameter            :: dp=selected_real_kind(12,200)
    integer                          nrows, ncols, i, j, ios
    double precision, allocatable :: dummy_double(:)

    ! PING non-Gaussianity cosmological parameters
    real,    optional :: alpha_e_in , H_e_in
    real                 a_e_PING   , H_e_PING

    ! Physical constants
    double precision, parameter :: e = 2.7182818284590452353602874713527_dp

    ! Option to also load in chi filter and chi-to-phi transfer function
    logical, optional :: Tchi2phi_flag_in, Tphi2zeta_flag_in
    logical              Tchi2phi_flag   , Tphi2zeta_flag

    ! Assign PING non-Gaussianity model cosmological parameters if not present
    if(present(alpha_e_in)) then
       a_e_PING = e ** alpha_e_in
    else
       a_e_PING = 1.0e-23
    endif
    if(present(H_e_in)) then
       H_e_PING = H_e_in
    else
       H_e_PING = 1.0e-10
    endif

    ! Check if the chi^2-to-phi transfer function is included in the table
    if(present(Tchi2phi_flag_in)) then
       Tchi2phi_flag = Tchi2phi_flag_in
    else
       Tchi2phi_flag = .False.
    endif

    ! Check if the chi filter W is included in the table
    if(present(Tphi2zeta_flag_in)) then
       Tphi2zeta_flag = Tphi2zeta_flag_in
    else
       Tphi2zeta_flag = .False.
    endif

    ! Determine the number of columns
    ncols = 2
    if(Tchi2phi_flag ==.True.) ncols = ncols + 1
    if(Tphi2zeta_flag==.True.) ncols = ncols + 1
    if(.not. allocated(dummy_double)) allocate(dummy_double(ncols))

    ! opens the file tabulating wavenumber and power spectrum and determines its length
    open(unit=2,file=fileTktab)
    nrows = 0
    do
       read(2,*,iostat=ios) (dummy_double(i), i=1,ncols)
       if(ios/=0) exit
       nrows = nrows + 1
    enddo
    close(2)

    ! Allocate size of transfer function array
    if(.not. allocated(T_dphi2dzeta)) allocate( T_dphi2dzeta(nrows,ncols) )
 
    ! Read the transfer function table
    open(unit=2,file=fileTktab)
    do i=1,nrows
       read(2,*) (dummy_double(j), j=1,ncols)
       do j=1,ncols
          T_dphi2dzeta(i,j) = sngl(dummy_double(j))
       enddo
    enddo
    close(2)

    ! ! Column 1 of T_dphi2dzeta is k/aH, so we multiply by aH to get it in the form of a k
    ! T_dphi2dzeta(:,1) = T_dphi2dzeta(:,1) * a_e_PING * H_e_PING

    ! Make sure no values of T(k)=0 for any k above the minumum k for which T(k)<=0
    j = 0
    do i=1,nrows
       if(T_dphi2dzeta(i,2) < 0.0) then
          if(j == 0) j=i
          if(j > 0)  T_dphi2dzeta(i,2) = 0.0
       endif
    enddo

  end subroutine read_Tktable

  ! Power spectrum for transverse inflaton field chi which imprints
  ! intermittent, spatially-localised non-Gaussianity peaks on the
  ! otherwise Gaussian zeta field 
  real function p_k_chi_extreme(wavenumber,A,B,R)
    implicit none
    real, parameter :: pi=4.0*atan(1.0)
    real wavenumber
    real A ! amplitude encoding strength of chi instability
    real B ! high wavenumber amplitude mitigating hard cutoff at k=2pi/R
    real R ! characteristic scale related to horizon size when instability
           ! takes place (during inflation so horizon size is small)

    ! Dimensionless power spectrum (power per ln k)
    p_k_chi_extreme=A*(wavenumber*R)**2.0*(B+exp(-wavenumber**2.0*R**2.0))

    ! Classical power spectrum in the form <|chi_k|^2>
    p_k_chi_extreme=p_k_chi_extreme*(2*pi**2)*wavenumber**-3

    return
  end function p_k_chi_extreme

subroutine read_table(filepktab, A_nonG, B_nonG)
    character(len=512) filepktab, p_chichi_file
    real, allocatable:: k_extreme(:), p_extreme(:)
    real A_nonG, B_nonG, R_nonG 
    real kchimin, kchimax, P0_low, n_low, P0_high, n_high
    real dummy_real, k_read, pk_read, i_interp
    integer i, nchi, i_low, offset_bytes
 
    ! Error handling
    integer :: ios, unit
    character(len=256) errmsg, line
    logical file_exists

  ! Check if the file exists
  inquire(file=filepktab, exist=file_exists)
  if (.not. file_exists) then
    print *, 'Error: File does not exist:', trim(filepktab)
    stop
  end if

  ! Open the file
  open(newunit=unit, file=filepktab, status='old', action='read', &
       iostat=ios, iomsg=errmsg)
  if (ios /= 0) then
    print *, 'Error opening file:', trim(filepktab)
    print *, 'Error code:', ios
    print *, 'Error message:', trim(errmsg)
    stop
  end if

  !print *, 'Successfully opened file:', trim(filepktab)
  !print *, 'Here are pksav_in values:'

  ! First, try to read the entire line as a string
  read(unit, '(A)', iostat=ios, iomsg=errmsg) line
  if (ios /= 0) then
    print *, 'Error reading first line of the file'
    print *, 'Error code:', ios
    print *, 'Error message:', trim(errmsg)
    stop
  end if

  !print *, 'First line of the file:'
  !print *, trim(line)

  ! Now try to parse the line
  read(line, *, iostat=ios, iomsg=errmsg) dummy_real, pksav_in(1), tsav_in(1), pkchisav_in(1)
  if (ios /= 0) then
    print *, 'Error parsing the first line of the file'
    print *, 'Error code:', ios
    print *, 'Error message:', trim(errmsg)
  !else
  !  print *, 'Successfully parsed the first line:'
  !  print *, 'dummy_real =', dummy_real
  !  print *, 'pksav_in(1) =', pksav_in(1)
  !  print *, 'tsav_in(1) =', tsav_in(1)
  !  print *, 'pkchisav_in(1) =', pkchisav_in(1)
  end if

  ! Rewind the file to start reading from the beginning
  rewind(unit)

  ! Now try reading the rest of the file
  if(NonGauss==0)then
    do i = 1, npktab
      read(unit, *, iostat=ios, iomsg=errmsg) dummy_real, pksav_in(i), tsav_in(i), pkchisav_in(i)
      tsav_in(i) = 0.0
      pkchisav_in(i) = 0.0
      if (ios /= 0) then
        if (ios > 0) then
          print *, 'Error reading file at record', i
          print *, 'Error message:', trim(errmsg)
        end if
        exit
      end if
      !print *, pksav_in(i)
    end do
  elseif(NonGauss>0) then
    do i=1,npktab
       !read(1,*) dummy_real,pksav_in(i),tsav_in(i),pkchisav_in(i)
       read(unit, *, iostat=ios, iomsg=errmsg) dummy_real, pksav_in(i), tsav_in(i), pkchisav_in(i)
       if(NonGauss==6) then
          pkchixsav_in(i)=p_k_chi_extreme(dummy_real,A_nonG,B_nonG,R_nonG)
       else
          pkchixsav_in(i)=0.0
       endif
       if (ios /= 0) then
         if (ios > 0) then
           print *, 'Error reading file at record', i
           print *, 'Error message:', trim(errmsg)
         end if
         exit
      end if
    enddo
    if((NonGauss>=7).and.(NonGauss<=11)) then

       ! Overwrites array for preheating P_chi 
       pkchisav_in=0.0

       ! Set the name of the chi power spectrum file
       if(NonGauss<=10) then
          p_chichi_file = 'tables/p_chichi.dat'
       elseif(NonGauss==11) then
          p_chichi_file = 'tables/p_chichi_e.dat'
       endif

       ! Read in chi power spectrum from Fortran-ordered 32-bit floats
       ! in unformatted binary file
       inquire(file=trim(p_chichi_file),size=nchi)
       nchi=nchi/4/2 ! p_chichi.dat has 4 bytes per value, 2 columns
       allocate(k_extreme(nchi))
       allocate(p_extreme(nchi))
       open(unit=2,file=trim(p_chichi_file),access='stream')
       do i=1,nchi
           offset_bytes = (i-1)*4*2
           read(2,pos=1+offset_bytes) k_extreme(i) ! k in Mpc^-1
           read(2,pos=5+offset_bytes) p_extreme(i) ! P_chich in Mpc^3 
       enddo

       if( (p_extreme(1)==0.0).or.(p_extreme(2)==0.0) ) then
           ! Extend the chi power spectrum as zeros if the spectrum read
           ! in has zero values at the low k
           n_low =1.0
           P0_low=0.0
       else
           ! Extend the chi power spectrum as a power law P(k) = P0 * k^n to
           ! cover the same k range as the other power spectra in Peak Patch
           n_low   = log( p_extreme(  2 )/p_extreme(   1  ) ) / &
                     log( k_extreme(  2 )/k_extreme(   1  ) )
           P0_low  = p_extreme(  1 )*exp(-n_low *log(k_extreme(  1 )))
       endif

       if( (p_extreme(nchi)==0.0).or.(p_extreme(nchi-1)==0.0) ) then
           ! Extend the chi power spectrum as zeros if the spectrum read
           ! in has zero values at the high k
           n_low =1.0
           P0_low=0.0
       else
           ! Extend the chi power spectrum as a power law P(k) = P0 * k^n to
           ! cover the same k range as the other power spectra in Peak Patch
           n_high  = log( p_extreme(nchi)/p_extreme(nchi-1) ) / &
                     log( k_extreme(nchi)/k_extreme(nchi-1) )
           P0_high = p_extreme(nchi)*exp(-n_high*log(k_extreme(nchi)))
       endif

       ! Interpolate
       do i=1,npktab
          k_read   = pktabmin * 10.**(real(i-1)*dpktab)
          i_interp = 1 + (nchi-1) * log(k_read         /k_extreme(1)) &
                                  / log(k_extreme(nchi)/k_extreme(1))
          i_low    = floor(i_interp)
          if(i_low<1) then
             ! Extrapolate fit to power law if k < k_extreme(1)
             pkchisav_in(i) = P0_low*k_read**n_low
          elseif(i_low+1>nchi) then
             ! Extrapolate fit to power law if k > k_extreme(nchi)
             pkchisav_in(i) = P0_high*k_read**n_high
          else
              if((p_extreme(i_low+1)==0.0).and.(p_extreme(i_low)==0.0)) then
                  pkchisav_in(i) = 0.0
              elseif((p_extreme(i_low+1)==0.0).and.(p_extreme(i_low)/=0.0)) then
                  ! linear interpolation if high-k value of P is zero
                  pkchisav_in(i) = -p_extreme(i_low) * &
                      (k_read-k_extreme(i_low+1))    / &
                      (k_extreme(i_low+1)-k_extreme(i_low))
              elseif((p_extreme(i_low+1)/=0.0).and.(p_extreme(i_low)==0.0)) then
                  ! linear interpolation if low-k value of P is zero
                  pkchisav_in(i) = p_extreme(i_low+1) * &
                      (k_read-k_extreme(i_low))       / &
                      (k_extreme(i_low+1)-k_extreme(i_low))
              else
                  ! loglog interpolate if k \in [k_extreme(1),k_extreme(nchi)]
                  pkchisav_in(i) = p_extreme(i_low) *                &
                      (k_read/k_extreme(i_low)) ** (                 &
                      log( p_extreme(i_low+1) / p_extreme(i_low) ) / &
                      log( k_extreme(i_low+1) / k_extreme(i_low) ) )
              endif
          endif
          ! Write out chi powerspectrum if you need to check it
          if(pkchisav_in(i)<0.0) write(*,*) k_read, pkchisav_in(i)
       enddo
    endif!(NonGauss==7)
  endif
  close(1)
end subroutine read_table

end module pktable
