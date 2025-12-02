#define METRIC_PERT 0

#define NFIELD 2
#define SYS_DIM_BG 2+2*nfld
#define SYS_DIM_PERT SYS_DIM_BG+nfld*(2*nfld+1)

#define SYS_BG_I 1:2+2*nfld

#define ALPHA_I 1
#define HUB_I 2
#define PHI_I 3:2+nfld
#define DPHI_I 3+nfld:2+2*nfld

#define L_I SYS_DIM_BG+1:SYS_DIM_BG+nfld*(nfld+1)/2
#define DL_I SYS_DIM_BG+nfld*(nfld+1)/2+1:SYS_DIM_BG+nfld*(nfld+1)
#define CORR_XDX_I SYS_DIM_BG+nfld*(nfld+1)+1:SYS_DIM_BG+nfld*(3*nfld+1)/2
#define CORR_DXDX_I SYS_DIM_BG+nfld*(3*nfld+1)/2+1:SYS_DIM_BG+nfld*(2*nfld+1)

#define CORR_FF_I 1:2*(nfld-1)+1:2,1:2*(nfld-1)+1:2
#define CORR_FDF_I 1:2*(nfld-1)+1:2,2:2*(nfld-1)+2:2
#define CORR_DFDF_I 2:2*(nfld-1)+2:2,2:2*(nfld-1)+2:2
