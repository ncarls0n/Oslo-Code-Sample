#define GL_ORDER 8

#if (GL_ORDER == 2)
#define GL_SET_A reshape((/a11_gl2/),(/1,1/))
#define GL_SET_B (/b1_gl2/)
#elif (GL_ORDER == 4)
#define GL_SET_A reshape((/a11_gl4, a21_gl4, \
			  a12_gl4, a22_gl4/),(/2,2/))
#define GL_SET_B (/b1_gl4, b2_gl4/)
#elif (GL_ORDER == 6)
#define GL_SET_A reshape((/a11_gl6, a21_gl6, a31_gl6, \
			  a12_gl6, a22_gl6, a32_gl6, \
			  a13_gl6, a23_gl6, a33_gl6/),(/3,3/))
#define GL_SET_B (/b1_gl6, b2_gl6, b3_gl6/)
#elif (GL_ORDER == 8)
#define GL_SET_A reshape((/a11_gl8, a21_gl8, a31_gl8, a41_gl8, \
			  a12_gl8, a22_gl8, a32_gl8, a42_gl8, \
			  a13_gl8, a23_gl8, a33_gl8, a43_gl8, \
			  a14_gl8, a24_gl8, a34_gl8, a44_gl8/),(/4,4/))
#define GL_SET_B (/b1_gl8, b2_gl8, b3_gl8, b4_gl8/)
#endif

#define RECURSEFIX 1
#define CONVTEST 0
#define CONVTEST_NEWT 0
