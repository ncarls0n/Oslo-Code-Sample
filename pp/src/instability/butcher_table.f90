  ! butcher_table.f90

module butcher_table
#include "gl_macros.h"
  use params

  implicit none

#if (GL_ORDER==2)
  real(dl), parameter :: a11_gl2 = 0.5_dl
  real(dl), parameter :: b1_gl2 = 1._dl
#elif (GL_ORDER==4)
  real(dl), parameter :: a11_gl4 = 0.25_dl, &
                         a12_gl4 = -0.0386751345948128822_dl, &
                         a21_gl4 = 0.53867513459481288_dl, &
                         a22_gl4 = 0.25_dl
  real(dl), parameter :: b1_gl4 = 0.5_dl, &
                         b2_gl4 = 0.5_dl
#elif (GL_ORDER==6)
  real(dl), parameter :: a11_gl6 = 0.138888888888888_dl, &
                         a12_gl6 = -0.0359766675249389_dl, &
                         a13_gl6 = 0.00978944401530832_dl, &
                         a21_gl6 = 0.300263194980964_dl, &
                         a22_gl6 = 0.222222222222222_dl, &
                         a23_gl6 = -0.022485417203086_dl, &
                         a31_gl6 = 0.267988333762469_dl, &
                         a32_gl6 = 0.480421111969383_dl, &
                         a33_gl6 = 0.138888888888888_dl
  real(dl), parameter :: b1_gl6 = 0.277777777777777_dl, &
                         b2_gl6 = 0.444444444444444_dl, &
                         b3_gl6 = 0.277777777777777_dl
#elif (GL_ORDER==8)
  real(dl), parameter :: a11_gl8 = 0.08696371128436346434_dl, &
                         a12_gl8 = -0.02660418008499879331_dl, &
                         a13_gl8 = 0.01262746268940472451_dl, &
                         a14_gl8 = -0.00355514968579568317_dl, &
                         a21_gl8 = 0.18811811749986807163_dl, &
                         a22_gl8 = 0.16303628871563653565_dl, &
                         a23_gl8 = -0.02788042860247089521_dl, &
                         a24_gl8 = 0.00673550059453815551_dl, &
                         a31_gl8 = 0.16719192197418877317_dl, &
                         a32_gl8 = 0.35395300603374396651_dl, &
                         a33_gl8 = 0.16303628871563653565_dl, &
                         a34_gl8 = -0.01419069493114114295_dl, &
                         a41_gl8 = 0.17748257225452261185_dl, &
                         a42_gl8 = 0.31344511474186834679_dl, &
                         a43_gl8 = 0.35267675751627186461_dl, &
                         a44_gl8 = 0.08696371128436346434_dl
  real(dl), parameter :: b1_gl8 = 0.17392742256872692868_dl, &
                         b2_gl8 = 0.32607257743127307130_dl, &
                         b3_gl8 = 0.32607257743127307130_dl, &
                         b4_gl8 = 0.17392742256872692868_dl
#elif (GL_ORDER==10)
  real(dl), parameter :: a11_gl10 = , &
                         a12_gl10 = , &
                         a13_gl10 = , &
                         a14_gl10 = , &
                         a15_gl10 = , &
                         a21_gl10 = , &
                         a22_gl10 = , &
                         a23_gl10 = , &
                         a24_gl10 = , &
                         a25_gl10 = , &
                         a31_gl10 = , &
                         a32_gl10 = , &
                         a33_gl10 = , &
                         a34_gl10 = , &
                         a35_gl10 = , &
                         a41_gl10 = , &
                         a42_gl10 = , &
                         a43_gl10 = , &
                         a44_gl10 = , &
                         a45_gl10 = , &
                         a51_gl10 = , &
                         a52_gl10 = , &
                         a53_gl10 = , &
                         a54_gl10 = , &
                         a55_gl10 = 
  real(dl), parameter :: b1_gl10 = , &
                         b2_gl10 = , &
                         b3_gl10 = , &
                         b4_gl10 = , &
                         b5_gl10 = 
#endif
end module butcher_table
