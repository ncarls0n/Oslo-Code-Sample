#!/bin/sh
disk=/Users/bond/Dropbox/bdrop/pk_src
computer=_MAC.x

Univer=pyk_LCDM_n096_s08

echo Univer $Univer

case $Univer in
pyk_LCDM_n096_s08)
densICfile=${disk}/pkp_transfer_n096_s08.030 ;
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
oCDM_b1.1)
densICfile=${disk}/decnu/Gam0217_n.030 ;
Ulist=0.3189,0.050,0.0,0.0,0.0,0.6311,0.7;;
oCDM_b1.1L)
densICfile=${disk}/decnu/Gam02_n.030 ;
Ulist=0.2849,0.050,0.0,0.0,0.0,0.6651,0.7 ;;
esac
echo Ulist $Ulist

run_hom_sphere_tab${computer} << EOF
`echo ${Ulist}`  ! Omx OmB Omhdm,fhdmclus, Omvac Omcurv h 0.05 0.0 0.95 0 0.8 
1            ! make an F-z table
1            ! icheck_tab
0.005        ! tfac /fraction of local 1-axis Hubble time for dt / 
101          ! zinit_fac
1            ! ivir_strat 
0.01         ! fcoll_3
2000         ! critical overdensity
50          ! nout - not needed for table building only
1.686        ! Frho
2.11         ! Frho
2.686        ! Frho
3.86         ! Frho
4.58164930        ! Frho
5.86         ! Frho
16.86        ! Frho
-1.11         ! Frho
20 0.8 0.001 ! Frho_max, Frho_min, dlogFrho
0.005        ! BEGIN CHECK PART / tfac / 
101          ! zinit
1            ! ivir_strat 
0.01         ! fcoll_3
2000         ! critical overdensity
50          ! nout - not needed for table building only
1.686        ! Frho
16.86        ! Frho
2.11         ! Frho
2.686        ! Frho
3.86         ! Frho
4.86         ! Frho
5.86         ! Frho
30           ! Frho
-1           ! Frho
EOF
exit



