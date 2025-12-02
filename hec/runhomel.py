#!/usr/bin/env python
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Usage 
"""
Starting from the peak patch repository

cd src/modules/HomogeneousEllipsoid/
python3 runhomel.py
"""

# Homogeneous Ellipsoid initial parameters
Frho = 1.686    # overdensity
e_v  = .1       # ellipticity, e_v >= 0
p_v  = -.05     # prolateness, -e_v <= p_v <= e_v
z    = (10.,0.) # (initial redshift , final redshift)
                # where initial redshift <~ 1100.0
                # and final redshift < initial redshift

# Compilers
machine = subprocess.check_output('hostname')
machines = { b'homes-MacBook-Air.local\n' :'NateMacBookAir',
             b'Nates-MacBook-Pro.local\n' :'NateMacBookPro',
             b'nia-login01.scinet.local\n':'Niagara'       ,
             b'sheep\n'                   :'Sheep'          }
fortrancompilers = { 'NateMacBookAir':'gfortran',
                     'NateMacBookPro':'gcc'     ,
                     'Niagara'       :'ifort'   ,
                     'Sheep'         :'gfortran' }
if machine in machines: 
    fortrancompiler = fortrancompilers[machines[machine]]
else: 
    fortrancompiler = 'gfortran'

# Copiler options
fortranoptions = '-O4 -w'
    # -Oi is optimization with -O1 being minimum and -O4 being maximum
    # -w flag inhibits all warnings from the compiler

# Directories that HomogeneousEllipsoid.f90 makes calls to
dir_ = {'homel':'./'}
dir_['src']             = dir_['homel']+'../../'
dir_['External']        = dir_['src']+'modules/External/'
dir_['GlobalVariables'] = dir_['src']+'modules/GlobalVariables/'
dir_['Solvers']         = dir_['src']+'modules/Solvers/'
dir_['cosmology']       = dir_['src']+'cosmology/'

# HomogeneousEllipsoid dependencies
f90s = { 'intreal_types'        : dir_['External']        ,
         'params'               : dir_['GlobalVariables'] ,
         'cosmoparams'          : dir_['GlobalVariables'] ,
         'input_parameters'     : dir_['GlobalVariables'] ,
         'Solvers'              : dir_['Solvers']         ,
         'Dlin_params'          : dir_['cosmology']       ,
         'psubs_Dlinear'        : dir_['cosmology']       ,
         'HomogeneousEllipsoid' : dir_['homel']           , 
         'runhomel'             : dir_['homel']            }

# All fortran scripts to compile and link
scripts_f90 = ''
for i in f90s:
    scripts_f90 += f90s[i]+i+'.f90 '

# Comment out calls to mpi_finalize(ierr) in HomogeneousEllipsoid.f90
os.system('''
    cd {0}
    touch temp.f90
    sed 's/call mpi_finalize(ierr)/!call mpi_finalize(ierr)/' \\
        HomogeneousEllipsoid.f90 > temp.f90
    rm -f HomogeneousEllipsoid.f90
    mv temp.f90 HomogeneousEllipsoid.f90
    '''.format( dir_['homel'] ))
    # This feature is used for parallel Peak Patch runs on Niagara, but we
    # don't need it for single homel runs.

# Change mode to avoid issues with Peak Patch .f90 scripts being executable
os.system('cd {0};chmod -x {1}'.format( dir_['homel'], scripts_f90 ))

# Make logfile
filename = 'runhomel.log'
os.system('cd {0};rm -f {1};touch {1}'.format( dir_['homel'], filename ))

# Compile and Link .f90 files
exe_homel = 'runhomel' # Executable file name
os.system(f'''
    cd {dir_['homel']}
    {fortrancompiler} {fortranoptions} {scripts_f90} -o {exe_homel}
    ''')

""" SIMULATION PARAMETERS """
# Cosmological parameters
params = {
    # homel parameters, declared in params.f90   
    'iwant_evmap' : 4,    # 1) turn z_vir vs e_v, 2) z vs p_v for fixed e_v, 3) table
    'nstepmax'    : 10000,# Stops neverending looping, using same value as hpkvd
    'iwant_rd'    : 1,    # 1) Use Carlson's elliptic integrals, else) don't use them
    'tfac'        : 0.01, # fraction of local 1-axis Hubble time for dt
    'zinit'       : z[0], # initial redshift of homel simulation
    'dcrit'       : 200., # Delta_critical, critical overdensity
    'e_vmax'      : 0.,   # e_vmax and de_v only used if iwant_evmap==1
    'de_v'        : 0.,   #
    'p_vbar'      : 0.,   # p_vbar, dp_v and e_vbar only used if
    'dp_v'        : 0.,   #     iwant_evmap==2 
    'e_vbar'      : 0.,   # 
    'Fbar'        : 0.,   # Fbar only used if iwant_evmap==1 or 2
    'fcoll_3'     : .171, # f_{coll,i} radial freeze-out factors for each
    'fcoll_2'     : .171, #     axis of the ellipsoid (see note below).
    'fcoll_1'     : .01,  # 
    'ivir_strat'  : 2,    # 1) a_jeq=fcoll_3 a_b3, 2) a_jeq=fcoll_j a_bj
    'iforce_strat': 4,    # 0) no bg, 1) sbg, 3) bg+NLstrain, 4) stbg+Lstrain, 5) Lstrain, 6) SW b_i
    # The radial freeze-out factors f_{coll,i} relates to the background
    # scale factor a, the scale factor a_i for each axis (and corresponding
    # semiaxes R_i) of the ellipsoid found at a filter scale R_f by 
    #           a_{coll,i} = R_{coll,i}/R_f = a f_{coll,i}
    # The values for the freeze-out are chosen to best recover the halo
    # mass function. The first two axes halt when they reach collapse con-
    # sistent with a final dcrit^{-1/3}=200^{-1/3}=.171, and the final axis
    # collapses to near-completion .01~0 (but stops at a sufficiently large
    # value that it does not introduce unwanted numerical effects in
    # solving ODEs), see Stein, Alvarez & Bond 2018.

    # Cosmological parameters, declared in cosmoparams.f90
    'Omb'  : 0.0493, # Omega_baryon, baryonic matter density fraction
    'Omx'  : 0.2645, # Omega_x, cold DM density fraction
    'Omvac': 0.6862, # Omega_Lambda, DE density fraction
    'Omcur': 0.0,    # Omega_k, curvature parameter
    'h'    : 0.6735, # little h (dimensionless Hubble constant), h = H_0/(100 km s^-1 Mpc^-1)
    #    'ns'    : 0.9649, # n_s, spectral index
    #    'sigma8': sigma(z=0,r=8 Mpc)
    # Values from Planck 2018 results. VI Cosmological parameters
    # Table 2, column TT,TE,EE+lowE+lensing

    # i/o parameters, declared in input_parameters.f90
    'ihard':1, # used in formatting output, assigned 1 in hpkvd.f90, so I did the same here
    }
params_to_read=str(len(filename))+'\n'+filename+'\n'
for i in params:
    params_to_read+=str(params[i])+'\n'

# Run the executable runhomel
os.system('''
    cd {0}
    ./{1} << EOF
    {2}
    {3} {4} {5}
    EOF'''.format( dir_['homel'], exe_homel, params_to_read,
        Frho, e_v, p_v ))

# Read runhomel.log
raw = np.loadtxt(filename, delimiter=',')
redshift = raw[:,0]**-1-1

# Eigenvalues of traceless strain tensor
lam1prime = -np.log( raw[:,1]*(raw[:,1]*raw[:,2]*raw[:,3])**(-1/3) )
lam2prime = -np.log( raw[:,2]*(raw[:,1]*raw[:,2]*raw[:,3])**(-1/3) )
lam3prime = -np.log( raw[:,3]*(raw[:,1]*raw[:,2]*raw[:,3])**(-1/3) )

# Plot of a_i(z) and delta(z)
figa,axa1 = plt.subplots()
axa1.plot(redshift, raw[:,1], label=r'$a_1(z)$', ls='-', color='b')
axa1.plot(redshift, raw[:,2], label=r'$a_2(z)$', ls='--', color='b')
axa1.plot(redshift, raw[:,3], label=r'$a_3(z)$', ls=':', color='b')
axa1.invert_xaxis()
axa1.set_xlabel(r'redshift $z$')
axa1.set_ylabel(r'principle axis scale factor $a_i$', color='b')
axa1.tick_params(axis='y',colors='b')
#axa1.set_yscale('log')
#axa1.set_xscale('log')
axa2 = axa1.twinx()
axa2.plot(redshift, raw[:,4], label=r'$1+\delta_c(z)$', color='r')
axa2.set_ylabel(r'overdensity $\delta$', color='r')
axa2.set_yscale('log')
axa2.tick_params(axis='y',colors='r')
axa1.legend(loc='upper left')
axa1.set_title(r'Inititial conditions $\delta=$'+str(Frho)+r', $e_v=$'+str(e_v)+r', $p_v=$'+str(p_v))
# # Scale factor alone
# axa2.plot(redshift, raw[-1,4]*raw[:,0]**3)

# Plot of e_v(z) and p_v(z)
figb,axb1 = plt.subplots()
axb1.plot(redshift, raw[:,5], label=r'$e_v(z)$')
axb1.plot(redshift, raw[:,6], label=r'$p_v(z)$')
axb1.plot(redshift, np.log( raw[:,1]/raw[:,3] ) / (2*raw[:,4]), label='test e')
axb1.plot(redshift, np.log( raw[:,2]**2/raw[:,1]/raw[:,3] ) / (2*raw[:,4]), label='test p')
axb1.invert_xaxis()
axb1.set_xlabel(r'redshift $z$')
axb1.set_ylabel(r'eigenvalue')
axb1.legend()


figd,axd1 = plt.subplots()
axd1.plot(redshift, raw[:,7], label=r'$\lambda_1(z)$')
axd1.plot(redshift, raw[:,8], label=r'$\lambda_2(z)$')
axd1.plot(redshift, raw[:,9], label=r'$\lambda_3(z)$')
axd1.plot(redshift, raw[:,7]+raw[:,8]+raw[:,9], label=r'$\delta(z)=\sum_i \lambda_i(z)$')
axd1.plot(redshift, raw[:,4], label=r'$\delta(z)$', ls='--')
axd1.set_yscale('log')
axd1.invert_xaxis()
axd1.legend()
plt.show()




#plt.figure(7)
#plt.plot( raw[1:,0], raw[1:,0]-raw[:-1,0], ls='none', marker='.' )
#plt.show()

# Once execution is complete, clean up stray .mod files
scripts_mod = ''
for i in f90s:
    scripts_mod += i+'.mod '
os.system('''
    cd {0}
    rm -f {1}
    '''.format( dir_['homel'], scripts_mod ))

# Restore calls to mpi_finalize(ierr) in HomogeneousEllipsoid.f90
os.system('''
    cd {0}
    touch temp.f90
    sed 's/!call mpi_finalize(ierr)/call mpi_finalize(ierr)/' \\
        HomogeneousEllipsoid.f90 > temp.f90
    rm -f HomogeneousEllipsoid.f90
    mv temp.f90 HomogeneousEllipsoid.f90
    '''.format( dir_['homel'] )) # So that we don't mess up Peak Patch
