import numpy as np
import matplotlib.pyplot as plt
#plt.rcParams['text.usetex']=True
import sys,os
import pykDict

###########################################################################
#                                                                         #
#                        COMPARE CHI POWER SPECTRA                        #
#                                                                         #
# This script compares primordial chi power spectra for different runs.   #
#                                                                         #
# USAGE:                                                                  #
#     python3 <...>P_chi_compare.py <run-1> <run-2> ...                   #
#                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pphome = os.environ["PP_DIR"]

fig,ax = plt.subplots()

for i in range(1,len(sys.argv)):

    run_dir    = sys.argv[i]
    param_file = run_dir+'/param/param.params'
    pk_file    = run_dir+'/tables/p_chichi.dat'
    hpkvdparam = run_dir+'/hpkvd_params.bin'

    # Read from parameterfile into dictionary
    params=pykDict.pykDict()
    params.read_from_file(param_file)

    # Peak Patch geometry
    n = params['next']
    s = params['boxsisze']

    # Early universe parameters
    m_phi   = params['m_phi']
    m_chi   = params['m_chi']
    phi_w   = params['phi_w']
    m_tach  = params['m_tach']
    a_e     = params['a_e']
    a_prime = np.fromfile(hpkvdparam,dtype=np.float32,count=41)[-1]
    l_code  = 2.625907e-52

    # Read in power spectrum
    count = int( os.path.getsize(pk_file)/4/2 )

    pk = np.fromfile(pk_file,dtype=np.float32,count=count)
    pk = np.reshape(pk,(int(count/2),2))

    # # Convert k to early universe units (comment out to display k in LSS units)
    # pk[:,0] /= a_e/a_prime/l_code

    # Actually plot the things
    ax.plot( pk[:,0] , pk[:,1] , label='mlambda = '+str(m_tach) , ls='none', marker='.' )

ax.set_xscale('log');ax.set_yscale('log')
fig.set_size_inches(14,12)
fig.savefig( '{0}/lambda_chi_compare.png'.format(run_dir) , dpi=100 )
ax.legend()
plt.show()
