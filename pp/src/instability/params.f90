! params.f90

module params
#include "macros.h"

  integer, parameter :: dl = kind(1.d0)
  integer, parameter :: nfld = NFIELD
  real(dl), parameter :: mpl = 1.e5
  
end module params
