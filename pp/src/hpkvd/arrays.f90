module arrays

  use intreal_types

  integer, parameter :: n1         =  N_REPLACE
  integer, parameter :: n2         =  n1
  integer, parameter :: n3         =  n1
  integer, parameter :: np1d       =  n1
  integer, parameter :: np2d       =  n1
  integer, parameter :: np3d       =  n1
  integer, parameter :: nFmax      =  n1*n2*n3
  integer, parameter :: ircmax     =  max(n1,n2,n3)/2

  integer, parameter :: itabmax    =  2001, &
                        np         =  (np1d)**3, &
                        npt        =  (np1d+1)**3
  integer, parameter :: Npv_evmax =   10,   &
                        Npv_evmin  = -10,   &
                        Nevmax     =  200,   &
                        Nfscmax    =  200
  integer, parameter :: nbormax    =  500000
  integer, parameter :: nodemax    =  100000000
  integer, parameter :: ntabmax    =  2000

  integer, parameter :: npkmax=10
  integer, parameter :: npkmaxl=10000000
  integer, parameter :: nclmax=513
  integer, parameter :: nboxmax=100000
  integer, parameter :: nhalomax=163

  ! LOCAL GRID VARIABLES
  real(C_FLOAT),            pointer :: delta(:,:,:),delta_u(:,:,:)
  real(C_FLOAT),            pointer :: lapd(:,:,:),lapd_u(:,:,:)
  real(C_FLOAT),            pointer :: etax(:,:,:), etay(:,:,:), etaz(:,:,:)
  real(C_FLOAT),            pointer :: eta2x(:,:,:), eta2y(:,:,:), eta2z(:,:,:)

  complex(C_FLOAT_COMPLEX), pointer :: deltac(:,:,:),deltac_u(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: lapdc(:,:,:), lapdc_u(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: etaxc(:,:,:), etayc(:,:,:), etazc(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: eta2xc(:,:,:), eta2yc(:,:,:), eta2zc(:,:,:)

  real(C_FLOAT), pointer :: F(:,:,:)

  type(C_PTR)            :: deltap, delta_up, etaxp, etayp, etazp
  type(C_PTR)            :: lapdp,lapd_up,  eta2xp, eta2yp, eta2zp

  ! GLOBAL GRID VARIABLES
  real(C_FLOAT),            pointer :: deltag(:,:,:),lapdg(:,:,:)
  real(C_FLOAT),            pointer :: etaxg(:,:,:), etayg(:,:,:), etazg(:,:,:)
  real(C_FLOAT),            pointer :: eta2xg(:,:,:),eta2yg(:,:,:),eta2zg(:,:,:)

  complex(C_FLOAT_COMPLEX), pointer :: deltagc(:,:,:),lapdgc(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: etaxgc(:,:,:), etaygc(:,:,:),etazgc(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: eta2xgc(:,:,:),eta2ygc(:,:,:),eta2zgc(:,:,:)

  type(C_PTR)                       :: deltagp, etaxgp, etaygp, etazgp
  type(C_PTR)                       :: lapdgp, eta2xgp, eta2ygp,eta2zgp

  ! CLOUD
  integer, allocatable :: irs2(:)
  integer(kind=2), allocatable :: ixsvec(:,:)
  integer(kind=1), allocatable :: mask(:,:,:)

  ! HALOS
  integer, allocatable :: iwas(:,:),lagrange(:,:)

  real, allocatable :: xpk(:,:),ypk(:,:),zpk(:,:)
  real, allocatable :: vx2pk(:,:),vy2pk(:,:),vz2pk(:,:)
  real, allocatable :: vxpk(:,:),vypk(:,:),vzpk(:,:)
  real, allocatable :: strain_bar(:,:,:)
  real, allocatable :: Fcollv(:,:),RTHLv(:,:),vTHvir(:,:)
  real, allocatable :: hp(:,:),F_ev(:,:),F_pv(:,:), F_d2(:,:)
  real, allocatable :: FcollvRf(:,:),F_d2Rf(:,:)
  real, allocatable :: grad(:,:,:),gradrf(:,:,:),gradf(:,:,:)
  real, allocatable :: zform(:,:)

  type(C_PTR) :: outbufferp

  ! BOXES
  integer, allocatable :: boxlist(:)
  real, allocatable :: xbox(:),ybox(:),zbox(:)

  ! INTEGRATION
  real akk(1100),dkk(1100)

contains

  subroutine allocate_halos(num_redshifts,ioutshear)

    implicit none
    
    integer num_redshifts, ioutshear
    
    allocate(iwas(npkmaxl,num_redshifts))
    allocate(lagrange(npkmaxl,num_redshifts))
    allocate(xpk(npkmaxl,num_redshifts))
    allocate(ypk(npkmaxl,num_redshifts))
    allocate(zpk(npkmaxl,num_redshifts))
    allocate(vxpk(npkmaxl,num_redshifts))
    allocate(vypk(npkmaxl,num_redshifts))
    allocate(vzpk(npkmaxl,num_redshifts))
    allocate(vx2pk(npkmaxl,num_redshifts))
    allocate(vy2pk(npkmaxl,num_redshifts))
    allocate(vz2pk(npkmaxl,num_redshifts))
    
    allocate(RTHLv(npkmaxl,num_redshifts))
    allocate(vTHvir(npkmaxl,num_redshifts))
    allocate(hp(npkmaxl,num_redshifts))
    
    allocate(Fcollv(npkmaxl,num_redshifts))

    RTHLv=0

    RTHLv=0
    
    if(ioutshear>0) then

       allocate(F_ev(npkmaxl,num_redshifts))
       allocate(F_pv(npkmaxl,num_redshifts))
       allocate(strain_bar(6,npkmaxl,num_redshifts))
       allocate(F_d2(npkmaxl,num_redshifts))
       allocate(zform(npkmaxl,num_redshifts))
       allocate(grad(3,npkmaxl,num_redshifts))
       allocate(gradf(3,npkmaxl,num_redshifts))

       allocate(FcollvRf(npkmaxl,num_redshifts))
       allocate(F_d2Rf(npkmaxl,num_redshifts))
       allocate(gradrf(3,npkmaxl,num_redshifts))

    endif
    
  end subroutine allocate_halos
  
  subroutine allocate_boxes
    
    allocate(xbox(nboxmax))
    allocate(ybox(nboxmax))
    allocate(zbox(nboxmax))
    allocate(boxlist(nboxmax))

  end subroutine allocate_boxes

end module arrays
