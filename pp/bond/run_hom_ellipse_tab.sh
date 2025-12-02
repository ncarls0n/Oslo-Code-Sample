#!/bin/sh
disk=/Users/bond
computer=_MAC.x
#computer=_DEC

Univer=pyk_LCDM_n096_s08
#LCDMgas_b1.1
#sCDM_b1.545
#oCDM_b1.1
#LCDM_b1.1

case $Univer in
pyk_LCDM_n096_s08)
powfile=/${disk}/bond/Dropbox/bdrop/pk_src/pkp_transfer_n096_s08.030;
Ulist=0.25,0.05,0.0,0.0,0.70,0.0,0.7;;
hCDM_b1.2)
densICfile=${disk}/decnu/hcdmOMnu02Nmnu2z3.030 ;
Ulist=0.75,0.050,0.20,0.0,0.00,0.0,0.5;;
sCDM_b1.5)
densICfile=${disk}/decnu/cdmb_81n1.030 ;
Ulist=0.95,0.050,0.0,0.0,0.00,0.0,0.5;;
sCDM_b1.545)  echo here here here ;
densICfile=${disk}/decnu/cdmb_81n1.030 ;
Ulist=0.95,0.050,0.0,0.0,0.00,0.0,0.5;;
sCDM_b1.5L)
densICfile=${disk}/decnu/cdmb_81n1.030 ;
Ulist=0.9489,0.050,0.0,0.0,0.0011,0.0,0.5;;
LCDM_b1.1)
densICfile=${disk}/decnu/Gam02_n.030 ;
Ulist=0.2849,0.050,0.0,0.0,0.6651,0.0,0.7;;
LCDMgas_b1.1)
densICfile=${disk}/decnu/utfn_LCDMgas.030 ;
Ulist=0.2592,0.00408,0.0,0.0,0.70,0.0,0.7;;
oCDM_b1.1)
densICfile=${disk}/decnu/Gam0217_n.030 ;
Ulist=0.3189,0.050,0.0,0.0,0.0,0.6311,0.7;;
esac
echo Ulist $Ulist

echo running run_hom_ellipse_tab${computer}
run_hom_ellipse_tab${computer} << EOF
`echo ${Ulist}`  ! Omx OmB Omhdm fhdmclus Omvac Omcurv h 0.05 0.0 0.95 0 0.8 
3    iwant_evpvtab /1 READ, 2 MAKE, 3 MAKE-READ/
1           table check / 1 TABNAME 0 DUMMY on next line ${disk}/sphx/Peaks/zvir1F3.tab 
fc1inv_${Ulist}.tab 
9 1 1.0 25.0 20.0 0.05 0.20 0.7 0.9999  Nfsc,ifc1inv,Fsmin,Fsmax,zinit_fac_ell,devtab,dpv_evtab,evtabmax,pv_evtabmax
3.9 0.0 0.0             Frho,ev,pv 
1.9 0.0 0.0             Frho,ev,pv
1.9 0.19 0.05             Frho,ev,pv
1.9 0.1 0.01             Frho,ev,pv
2.9 0.4 0.1             Frho,ev,pv
2.2 0.5 -0.3             Frho,ev,pv
-1.9 -0.4 0.3             Frho,ev,pv
1            do you want rd for elliptic
3   iwant_evmap/ (1) turn/ vir z vs e_v (2)z vs p_v for fixed e_v (3) table/
1 want an e_v-p_v table?/dum if iwant_evmap=0/
4   iforce_strat /0 no bg,1 sbg,3 bg+NLstrain,4 stbg+Lstrain,5 Lstrain,6 SW b_i/
0.005        ! tfac /fraction of local 1-axis Hubble time for dt / 
101          ! zinit
2            ! ivir_strat /1 a_jeq=fcoll_3 a_b3,2 a_jeq=fcoll_j a_bj/
0.18         ! fcoll_3
100         ! critical overdensity
500          ! nout - not needed for table building only
1.9 0.1 0.1     Fbar,de_vtab,dp_vtab OR Fbar,e_vmax,de_v,p_vbar OR Fbar,e_vbar,dp_v
1            'hardcopy? /1 usual output, 2 for plot or tab/' TABFILE or DUMMY
deldel
.3 .1         input ev,pv /ev < 0 stop /
.2 .05         input ev,pv /ev < 0 stop /
-1 -1 
EOF
exit



