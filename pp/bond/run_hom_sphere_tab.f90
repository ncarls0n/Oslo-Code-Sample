PROGRAM run_hom_sphere
  call run_hom_sphere_tab
END PROGRAM run_hom_sphere
!C THE PROGRAM RUNS hom_sphere_tab AND CHECKS IT

subroutine run_hom_sphere_tab
  USE intreal_types

  common/evalues_sph/aLam_3,Frho
  common/params_sph/iwant_evmap,tfac,zinit,nstepmax,nyp,nout,dcrit,Fbar
  common/equil_sph/fcoll_3,a_3eq,lvirv,ivir_strat,iforce_strat
  common/hard/ihard
  common/Dlin_params/D_ainfty,D_a0,HD_Hainfty,HD_Ha0,chininfty,chin0
  character(LEN=80) zevoutfile
  external get_derivs_sph

  call read_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)
  call Dlinear_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)

  !c      if(iEdeS.ne.1) then
  !c         call make_F_zvir_table
  !c      endif

  write(*,*) 'z,a,chi :  chi check'
  dz=0.2
  z=0.1
313 continue
  z=z+dz
  a=1.0/z
  chi=chifn(a)
  write(*,*) z,a,chi
  if(z.lt.3) goto 313
  z=100
  a=1.0/z
  chi=chifn(a)
  write(*,*) z,a,chi

  write(*,*) ' do you want to make an F-z table? '
  read(*,*) iwant_Fac_table
  if(iwant_Fac_table.eq.1) then
     write(*,*) ' icheck_tab ' 
     read(*,*) icheck_tab 
     call make_F_zvir_table
     write(*,*) 'TABLE MADE'
  else
     icheck_tab=0
  endif
  if(iwant_Fac_table.eq.1.and.icheck_tab.ne.1) return

  nstepmax=50000
  ny=2
  nyp=ny
  write(*,*) ' tfac / fraction of local 1-axis Hubble time for dt/'
  read(*,*) tfac
  write(*,*) ' zinit_fac NOTE zinit=min(zinit_fac*Frho,zinit_fac)'
  read(*,*) zinit_fac
  iwant_evmap=0
  !c      write(*,*) ' iforce_strat / 0 no bgnd, 1 std bgnd, '
  !c      write(*,*) '.... 3 bgnd + nonlinear strain, 4 std bgnd ',
  !c     $     '+ linear strain, '
  !c      read(*,*) 
  iforce_strat=0
  write(*,*) ' ivir_strat / 1 a_jeq=fcoll_3 a_b3 / '
  read(*,*) ivir_strat
  write(*,*) ' fcoll_3 '
  read(*,*) fcoll_3

  write(*,*) 'critical overdensity /<= 0 stops/ '
  read(*,*) dcrit
  if(dcrit.lt.0.0) return
  if(icheck_tab.eq.1) then 
     iwant_evmap=2
  else
     iwant_evmap=0
  endif
  write(*,*) ' nout / output every nout steps / '
  read(*,*) nout 
110 continue
  if(iwant_evmap.eq.0.or.iwant_evmap.eq.2) then
     write(*,*) 'Frho / Frho.le.0 stop/'
     read(*,*) Frho
     write(*,*) 'Frho ', Frho
     if(Frho.lt.0.0) then
        if(iwant_Fac_table.eq.0) return
        zvir3=1.0
        dzvir3=0.2
        write(*,*) 'zvir3,Frho,Frhoac,FrhoDc,D_a,delta '
        delta_approx_now=178.0/Omnr**0.6
        write(*,*) ' delta_approx_now ',delta_approx_now
555     continue
        a_b=1.0/zvir3
        Dlin=Dlinear(a_b,chi,HD_Ha,D_a)
        Frho_tab=Frho_tabfn(zvir3)
        Frhoac_tab=Frho_tab/zvir3
        FrhoDc_tab=Frhoac_tab*D_a
        a_3vir_a_ba=a_3vir_tabfn(zvir3)*zvir3
        delta=1.0/a_3vir_a_ba**3
        !C THIS IS NOT THE APPROX  delta_approx=178.0*D_ainfty/D_a
        write(*,*) zvir3,Frho_tab,Frhoac_tab,FrhoDc_tab,D_a,delta 
        !c            write(*,*) ' delta,delta_approx ',delta,delta_approx
        zvir3=zvir3+dzvir3
        if(zvir3.le.20.0) goto 555
        return
     endif
  endif
  zinit=max(zinit_fac,zinit_fac*Frho)
!c is this max ok or shoule it be min
  call evolve_sphere(zvir3,zturn3,a_3ta,a_3vir_a_bvir,E_Mta_H02rpk2,EK_Mvir_H02rpk2_corr)
  if(iwant_Fac_table.eq.1) then 
     a_3vir=a_3vir_a_bvir/zvir3
     zvir3_tab=zvir3_tabfn(Frho)
     Frho_tab=Frho_tabfn(zvir3)
     a_3vir_tab=a_3vir_tabfn(zvir3)
     write(*,*) 'tab check',Frho,zvir3,Frho_tab,zvir3_tab,a_3vir,a_3vir_tab
  endif
  goto 110
end subroutine run_hom_sphere_tab

subroutine evolve_sphere(zvir3,zturn3,a_3ta,a_3vir_a_bvir,E_Mta_H02rpk2,EK_Mvir_H02rpk2_corr)
  USE intreal_types
  PARAMETER (ny=2)
  PARAMETER (one_third=1.0/3.0)
  dimension y(ny),dy(ny)
  common/params_sph/iwant_evmap,tfac,zinit,nstepmax,nyp,nout,dcrit,Fbar
  common/hard/ihard
  common/evalues_sph/aLam_3,Frho
  common/equil_sph/fcoll_3,a_3eq,lvirv,ivir_strat,iforce_strat
  common/univer/Omt,Omnr,Omvac,Omcurv,Rvac_nr,Rcurv_nr,Ha_b_nr
  external get_derivs_sph
  nyp=ny
  a_3eq=0.0
  lvirv=0
  t=1.0/zinit
  a3p=0.0            
  tp=0.0
  lvir_3p=0
  v3p2=0.0
  v3p=0.0
  ldcrit=0
  zturn3=-1
  zturn3p=-1
  delta1_ta=0.0
  delta1_tap=0.0
  zvir3=-1
  zdcrit=-1
  pi=4.0*atan(1.0)
  sq3=sqrt(3.0)
  delta1_ta0=9.0*pi*pi/16.0
  !c      delta1_ira_const=(2.0-1.831)/2.0
  delta1_ira_cnst=0.0845
  aLam_3 = Frho/3.0
  !C ROUGH ESTIMATE OF COLLAPSE zzc FROM aLam_1
  zzc_170sph=3.0*(aLam_3)/1.606
  zzc_Fsph=3.0*(aLam_3)/1.686

  call ic_set_sph(zinit,aLam_3,t,nyp,y)
  istep=0
  if(iwant_evmap.eq.0) then
     write(*,*)'  zz   a_b   a_3   Frho 1+delta vir:3'
  endif
111 continue
  istep=istep+1
  a_3=y(1)
  dtstep=tfac*a_3*sqrt(a_3)*Ha_b_nr
  !C NOT FUNNY RK4 OF PRESS REQUIRES FIRST CALL TO derivs
  call get_derivs_sph(t,y,dy)
  call rk4_sph(y,dy,nyp,t,dtstep,y,get_derivs_sph)
  t=t+dtstep
  delta1_e=t**3/(y(1)**3)
  if(a3p.ne.-1) then
     if(y(1).lt.a3p) then
        !c LINEAR INTERPOLATION: v3_interp=v3p+(v3-v3p)*(a_b_interp-a_bp)/(a_b-a_bp)
        !C    SET THIS TO ZERO
        a3=y(1)
        v3=y(2)
        pfac=-v3p/(v3-v3p)
        a_3ta=a3p+(a3-a3p)*pfac
        a_bta=tp+pfac*(t-tp)
        zturn3=1.0d0/a_bta
        delta1_ta=(a_bta/a_3ta)**3
        write(*,*) 'v3p,v3,pfac,a3p,a3,a_3ta,a_bp,a_b,a_bta'
        write(*,*) v3p,v3,pfac,a3p,a3,a_3ta,tp,t,a_bta
        if(omvac.eq.0.0) then
           a_3vir=0.5*a_3ta
        else
           call get_cubic_soln(a_3ta,a_3vir)
        endif
        a_3vir_a_3ta=a_3vir/a_3ta
        E_Mta_H02rpk2=-0.5*(Omnr/a_3ta+Omvac*a_3ta**2)*3.0/5.0
        EK_Mvir_H02rpk2_corr=0.5*(3.0*Omvac*a_3vir**2)*3.0/5.0
        write(*,*) 'a_3vir,a_3vir_a_3ta = ',a_3vir,a_3vir_a_3ta
        write(*,*) 'E_Mta_H02rpk2 = ',E_Mta_H02rpk2            
        write(*,*) 'EK_Mvir_H02rpk2_corr = ',EK_Mvir_H02rpk2_corr 
        !C ETA_LLPR=Lambda/(4 pi G rho_ta) = 2OMvac/OMnr/(delta1_ta*zturn3**3)
        eta_LLPR=OMvac/OMnr/delta1_ta/zturn3**3
        delta1_ta_ira=delta1_ta0*((1.0-delta1_ira_cnst*eta_LLPR)/(1.0-eta_LLPR))**sq3
        a3p=-1
        if(eta_LLPR.ne.0) then
           x_ira=2.0/eta_LLPR
        else 
           x_ira=0.0
        endif
     else 
        a3p=y(1)
        v3p2=v3p
        v3p=y(2)
        tp=t
        delta1_tap=delta1_ta
        delta1_ta=delta1_e
        zturn3p=zturn3
        zturn3=1.0d0/t
     endif
  endif
  if(lvirv.ne.lvir_3p) then
     zvir3=1.0/t
     lvir_3p=1
     delta1_e3=(t/y(1))**3
     delta1_vir=(t/a_3vir)**3
     a_3vir_a_bvir=a_3vir*zvir3
     write(*,*) 'delta1_ta,delta1_vir =',delta1_ta,delta1_vir
     write(*,*) 'a_3vir_a_bvir,a_3ta/a_bta =',a_3vir_a_bvir,a_3ta/a_bta
  endif
  if(iwant_evmap.eq.0) then
     if((istep.eq.(istep/nout)*nout).or.lvirv.eq.1.or.istep.ge.nstepmax) then
        zz=1.0/t
        a_b=1.0/zz
        delta1_e=t**3/(y(1)**3)
        Dlin=Dlinear(a_b,chi,HD_Ha,D_a)
        Frhoz=Frho*Dlin
        write(*,11) istep,zz,a_b,y(1),Frhoz,delta1_e,lvirv
11      format(i5,5(1pe11.3),2x,2i2)
     endif
  endif
  if(lvirv.eq.1.or.istep.ge.nstepmax.or.t.ge.2.0) then
     if(iwant_evmap.eq.2) return
     if(iwant_evmap.eq.0) write(*,*) 'last axis (1) is virialized '
     write(*,*) '  '
     write(*,67) Frho,aLam_3
67   format(' F,Lam123: ',6f9.3)
     write(*,66) zturn3,zvir3,zdcrit
66   format('zturn,vir,dcr: ',3(f9.3))
     write(*,68) zzc_170sph,zzc_Fsph,delta1_e3,dcrit
68   format('zc170sph,zcF: ',2(F9.3),' dvir,dcr: ',2(F9.3))
     write(*,*) 'zturn3,zturn3p ', zturn3,zturn3p 
     write(*,*) 'v3p,v3p2 at calculated turnaround ',v3p,v3p2 
     write(*,*) 'delta1_ta,delta1_tap,delta1_ta_ira ',delta1_ta,delta1_tap,delta1_ta_ira
     write(*,*) 'x_ira,eta_LLPR ',x_ira,eta_LLPR 
     if(ihard.eq.1.and.iwant_evmap.ne.3) then
        write(10,67) Frho,aLam_3
        write(10,66) zturn3,zvir3,zdcrit
        write(10,68) zzc_3sph,zzc_Fsph,delta1_e3,dcrit
     elseif(ihard.eq.2.and.iwant_evmap.ne.3) then
        write(10,69) Frho,zvir3,zdcrit,zzc_Fsph
69      format(8F9.3)
     endif
     return
  endif
  goto 111
end subroutine evolve_sphere


