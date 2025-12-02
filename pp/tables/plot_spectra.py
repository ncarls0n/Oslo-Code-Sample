import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex']=True
from scipy.optimize import curve_fit

# The cosmology of peak patch is described by power spectra loaded into the
# simulaiton from tabulated data in one of the .dat files in the working
# directory for this script (i.e. peak-patch/tables/). The most up-do-date
# power spectrum is planck2018_powerspectrum.dat
# 
# These data are in an array with dimension 500x4 where the 4 columns are
# respectively wavenumber k, powerspectrum in fourier space P(k), transfer
# function multiplied by k squared k^2 T (where the transfer function
# relates scalar potential Phi to curvature perturbations zeta
#     Phi(k,t) \propto T(k,t_0) zeta(k,t)
# ) and the last column represents the chi power spectrum (where chi is...)

# Power spectra file
filename = 'planck2018_powerspectrum.dat' #'bbps_power_spec_NG.dat'

# Spectral index
ns=0.9650 # Note that using ns=0.9649 as in parameter file gives A slight
          # slope, and 0.9649995690890052 gives best fit for zero slope

# Read in powerspectrum data
data = np.loadtxt(filename, delimiter=' ')
k    = data[:,0] # k
Pk   = data[:,1] # P_k = | \delta_k |^2
k2T  = data[:,2] # k^2 T, where T is the transfer function
Pchi = data[:,3] # P_\chi (scaling uncertain)

print(len(k), k[0], k[-1])

# Transfer function
T = k2T*k**-2

# Since Pk \propto T^2 k^{n_s}, the following should be constant
A    = Pk * T**-2 * k**-ns #A    = Pk * k2T**-2 * k**(4-ns)
Abar = np.mean(A)
Pbar = Abar * T**2 * k**ns

Pzeta = k2T**-2 * Pk

# Plot of k vs element j
'''
fig0,ax0 = plt.subplots()
ax0.plot(np.arange(1,len(k)+1,1), k)
ax0.set_yscale('log')
ax0.set_xlabel(r'Fortran array index $j$')
ax0.set_ylabel(r'wavenumber $k$')
'''

# Plot A and Abar to show that it is roughly constant
'''fig1,ax1 = plt.subplots()
ax1.plot(k, A, label=(r'$k^{{-n_s}}T^{{-2}}\mathcal{{P}}(k),~n_s={:.3f}$'.format(ns)) )
ax1.plot(k, 0*k+Abar, label=r'$\bar{A}=$'+str(round(Abar*1e8,4))+r'$\cdot10^{-8}$' )
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel(r'$k$ $[h$ Mpc$^{-1}]$')
ax1.set_ylabel(r'$\mathcal{P}(k)$ $[(h^{-1}$ Mpc$)^3]$')
ax1.legend()'''

# Plot P(k) and T(k) (up to factors of 2pi)
fig2,ax2 = plt.subplots()

# Matter power spectrum
ax2.plot(k, #(2*np.pi)**(3/2)*
    Pk, label=r'$P_{\delta\delta}(k)$')
#ax2.plot(k, T, label=r'$T(k)$' )

# Zeta power spectrum
############ax2.plot(k, Pzeta, label=r'$P_{\zeta\zeta}(k) = \frac{2\pi^2}{k^3} A_s \left(\frac{k}{k_0} \right)^{n_s-1}$')#, $A_s=2.1\times10^{-9}$, $k_0=0.05$ Mpc, $n_s=0.9649$')

# Read peak patch initial conditions power
#pkp_ics = np.loadtxt('bears.txt')
#k_ics   = pkp_ics[:,0]
#Pk_ics  = pkp_ics[:,1]
#prefact = 2*np.pi**2/k_ics**3
#ax2.plot( k_ics, prefact*(48/128)**3*Pk_ics, label='bears')

'''
s_box  = 4000 # Mpc
n_eff  = 360
s_latt = s_box/n_eff # Mpc

zetag = np.loadtxt('zetag.txt')
ax2.plot( zetag[:,0], s_latt**-3*zetag[:,1], label='zetag', ls='--' )

#zetang = np.loadtxt('zetang_fnl5e3.txt')
#ax2.plot( zetang[:,0], zetang[:,1], label=r'zetang, $f_{{NL}}=5\cdot10^3$', ls='--' )

deltag = np.loadtxt('rhog.txt')
ax2.plot( deltag[:,0], s_latt**-3*deltag[:,1], label=r'deltag' )

#deltang = np.loadtxt('delta_fnl1e3.txt')
#ax2.plot( deltang[:,0], deltang[:,1], label=r'deltang, $f_{{NL}}=1\cdot10^3$', ls='--' )

chi = np.loadtxt('chi.txt')
ax2.plot( chi[:,0], s_latt**-3*chi[:,1], label=r'chi')
'''


#def P_0(n,s_latt):
#    return (2*np.pi)**-3 * n**-3 * s_latt**3
#
#sb  = 128
#n   = 32
#nb  = 10
#ne  = n-2*nb
#sl  = sb/ne
#pk  = np.loadtxt('Fvec_128Mpc_n128_nb40_nt1_partial_power.dat')
#ax2.plot( pk[:,0], pk[:,1], label='partial' )
#
#sb1 = 128
#n1  = 128
#nb1 = 40
#ne1 = n1-2*nb1
#sl1 = sb1/ne1
#pk1 = np.loadtxt('Fvec_128Mpc_n128_nb40_nt1_power.dat')
#ax2.plot( pk1[:,0], P_0(n1,sl1) * pk1[:,1], label='sn128nb40' )

'''
sb2 = 256
n2  = 128
nb2 = 40
ne2 = n2-2*nb2
sl2 = sb2/ne2
pk2 = np.loadtxt('Fvec_256Mpc_n128_nb40_nt1_power.dat')
ax2.plot( pk2[:,0], P_0(n2,sl2) * pk2[:,1], label='s256n128nb40' )

sb3 = 400
n3  = 128
nb3 = 40
ne3 = n3-2*nb3
sl3 = sb3/ne3
pk3 = np.loadtxt('Fvec_400Mpc_n128_nb40_nt1_power.dat')
ax2.plot( pk3[:,0], P_0(n3,sl3) * pk3[:,1], label='s400n128nb40' )

sb4 = 400
n4  = 236
nb4 = 20
ne4 = n4-2*nb4
sl4 = sb4/ne4
pk4 = np.loadtxt('Fvec_400Mpc_n236_nb20_nt1_power.dat')
ax2.plot( pk4[:,0], P_0(n4,sl4) * pk4[:,1], label='s400n236nb20' )

sb5 = 400
n5  = 400
nb5 = 40
ne5 = n5-2*nb5
sl5 = sb5/ne5
pk5 = np.loadtxt('Fvec_400Mpc_n400_nb40_nt1_power.dat')
ax2.plot( pk5[:,0], P_0(n5,sl5) * pk5[:,1], label='sn400nb40' )
'''

#sba = 4000
#na  = 500
#nba = 5
#nea = na-2*nba
#sla = sba/nea
#pka = np.loadtxt('Fvec_4000Mpc_n236_nb20_nt10_partial_power.dat')
#ax2.plot( pka[:,0], pka[:,1], label='partial' )
#
#sb6 = 4000
#n6  = 2000
#nb6 = 20
#ne6 = n6-2*nb6
#sl6 = sb6/ne6
#pk6 = np.loadtxt('Fvec_4000Mpc_n236_nb20_nt10_power.dat')
#ax2.plot( pk6[:,0], P_0(n6,sl6) * pk6[:,1], label='s4000n1960,nb20')

# pk4_= np.loadtxt('Fvec_400Mpc_n400_nb40_nt1_97051_power.dat')
# ax2.plot( pk4_[:,0], P_0(n4,sl4) * pk4_[:,1], label='sn400nb40, alternate seed' )

p = np.loadtxt('Fvec_4000Mpc_n236_nb20_nt10_ng3_fnl10_power.dat')
ax2.plot( p[:,0], p[:,1], label=r'Model 3, $f_{\rm{NL}}=10$',ls=':' )

p = np.loadtxt('Fvec_4000Mpc_n236_nb20_nt10_ng4_fnl1e3_power.dat')
ax2.plot( p[:,0], p[:,1], label=r'Model 4, $f_{\rm{NL}}=10^3$',ls='-.' )

p = np.loadtxt('Fvec_4000Mpc_n236_nb20_nt10_ng6_fnl7.5e2_power.dat')
ax2.plot( p[:,0], p[:,1], label=r'Model 6, $f_{\rm{NL}}=7.5\times10^2$', ls='--' )

ax2.plot( p[:,0], (.8111/.8298)**2*p[:,1], label=r'Model 6, $f_{\rm{NL}}=7.5\times10^2$ normalised t    o $\sigma_8=0.8111$', ls='--' )

p = np.loadtxt('Fvec_4000Mpc_n236_nb20_nt10_ng6_fnl1e4_power.dat')
ax2.plot( p[:,0], p[:,1], label=r'Model 6, $f_{\rm{NL}}=10^4$', ls='-')

ax2.plot( p[:,0], (.8111/1.5251)**2*p[:,1], label=r'Model 6, $f_{\rm{NL}}=10^4$ normalised to $\sigma_8=0.8111$', ls='-' )

'''
# Test from Fortran
test_sbox =400.
test_n    =236
test_nbuff=20
test_slatt=test_sbox/(test_n-2*test_nbuff)
p_test=np.loadtxt('Fvec_400Mpc_n236_nb20_nt1_power.dat')
ax2.plot( p_test[:,0], (2*np.pi)**-1.5*test_n**-3*test_slatt**-3*2**1.5*p_test[:,1], label='fortran s400n236' )



test2_sbox =128.
test2_n    =128
test2_nbuff=40
test2_slatt=test2_sbox/(test2_n-2*test2_nbuff)
p_test2=np.loadtxt('Fvec_128Mpc_n128_nb40_nt1_power.dat')
ax2.plot( p_test2[:,0], (2*np.pi)**-1.5*test2_n**-3*test2_slatt**-3*p_test2[:,1], label='fortran sn128')



test3_sbox =400.
test3_n    =400
test3_nbuff=40
test3_slatt=test3_sbox/(test3_n-2*test3_nbuff)
p_test3=np.loadtxt('Fvec_400Mpc_n400_nb40_nt1_power.dat')
ax2.plot( p_test3[:,0], (2*np.pi)**-1.5*test3_n**-3*test3_slatt**-3*2**-1.5*p_test3[:,1], label='fortran sn400')



# Test from python
# p_testpy=np.loadtxt('Fvec_400_power_from_py.txt')
# ax2.plot( p_testpy[:,0], test_n**-3*test_slatt**-3*p_testpy[:,1], label='python test')
'''


## Chaotic billiards chi power
#ax2.plot(k, Pchi, label=r'$\mathcal{P}_{\chi\chi}(k)$')
#ax2.plot(k, Pbar, label=r'$\langle\mathcal{P}_k\rangle = \bar{A} T^2 k^{n_s}$' )

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlabel(r'$k$ $[$Mpc$^{-1}]$')
ax2.set_ylabel(r'$P(k)$ $[$Mpc$^3]$')



# # P_m(k) as presented in Planck collab. papers
# k_planck = np.array([ 1e-4,4e-4,1e-3,3e-3,6e-3 ,9e-3,1.8e-2,3e-2,6e-2,1e-1,3e-1 ,6e-1,1 ,3,6 ])
# P_planck = np.array([ 400 ,1500,4000,1e4 ,1.7e4,2e4 ,2.2e4 ,2e4 ,1e4 ,6e3 ,8.2e2,2e2 ,62,5,1 ])
# #ax2.plot( k_planck, P_planck, label=r'Theory $P_m(k)$ original')
# ax2.plot( k_planck/.6736, P_planck*.6736**3, ls='none', marker='.', label=r'Theory $P_m(k)$')
# #ax2.plot( k, Pk*(2*np.pi)**2, label=r'$(2\pi)^2\mathcal{P}(k)$')
# #ax2.plot( k, Pk*(2*np.pi)**3, label=r'$(2\pi)^3\mathcal{P}(k)$')
# #ax2.plot( k, Pk*(2*np.pi)**1.5, label=r'$(2\pi)^{3/2}\mathcal{P}(k)$')
# #ax2.plot( k_planck, (k_planck-k_planck[0])**3+P_planck[0], label=r'$\propto k^3$')
# #ax2.plot( k_planck[:-1], (k_planck[:-1]-k_planck[-1])**-3+P_planck[-1], label=r'$\propto k^{-3}$')


# Fitting high- and low-k for each P(k) model
def linear_in_loglog(x,a,n):
    return a*x**n
'''
# Fitting matter power spectrum P_{\delta\delta}(k) from table
low_k  = curve_fit(linear_in_loglog, k[:50] , Pk[:50] ,
    np.array([5e4,ns]) )
high_k = curve_fit(linear_in_loglog, k[-25:], Pk[-25:],
    np.array([.5,-3]) )
ax2.plot(k[:200] , linear_in_loglog(k[:200] , *low_k[0] ),
    ls='--',label=r'$({:.3f}\cdot10^4)k^{{n_s-{:.3f}}}$'.format(
        low_k[0][0]*1e-4 , ns-low_k[0][1]
        )
    )
ax2.plot(k[-150:], linear_in_loglog(k[-150:], *high_k[0]),
    ls='--',label=r'$({:.3f})k^{{{:.3f}}}$'.format(
        high_k[0][0] , high_k[0][1]
        )
    )
'''
# My P(k) crude model using sigmoid
#ax2.plot(k, low_k[0][0]*k**ns / (1+np.exp((k-.1)*1e5)) + k**-3/(1+np.exp(-(k-.1)*1e2)) )

# P_m(k) from Planck papers
'''low_k_planck  = curve_fit(linear_in_loglog, k_planck[:3], P_planck[:3],
    np.array([5e4,ns]) )
high_k_planck = curve_fit(linear_in_loglog, k_planck[-3:], P_planck[-3:],
    np.array([.5,-3]) )
ax2.plot(k_planck[:5], linear_in_loglog(k_planck[:5] , *low_k_planck[0] ),
    ls='--',label=r'$({:.3f}\cdot10^6)k^{{n_s-{:.3f}}}$'.format(
        low_k_planck[0][0]*1e-6 , ns-low_k[0][1]
        )
    )
ax2.plot(k_planck[-5:], linear_in_loglog(k_planck[-5:] , *high_k_planck[0] ),
    ls='--',label=r'$({:.3f})k^{{{:.3f}}}$'.format(
        high_k_planck[0][0] , high_k_planck[0][1]
        )
    )'''

# Fit for P_\chi
'''p_chi_fit = curve_fit(linear_in_loglog, k, Pchi, np.array([2e-14,-3]))
ax2.plot(k, linear_in_loglog(k, *p_chi_fit[0]),ls='--',
    label=r'$({:.3f}\cdot10^{{-14}})k^{{{:.3f}}}$'.format(
        p_chi_fit[0][0]*1e14, p_chi_fit[0][1]))'''

# P_\chi model for spatially-localized intermittent non-Gaussianity
As=1.e-10
Bs=1.e-10
R =64

def Pchichi(k,a,b,r):
    return ( 2*np.pi**2/k**3 *
        a*(k*r)**2*( b + np.exp(-(k*r)**2) ) )

#############ax2.plot(k, Pchichi(k,As,Bs,R), label=
#############    r'$P_{{\chi_E\chi_E}}(k)=\frac{{2\pi^2}}{{k^3}}A(kR)^2\left[B+e^{{-(kR)^2}}\right]$')#, $A=B=10^{{-10}}$, $R=64$ Mpc' )

#As=4.096e-10
#R =0.04
#ax2.plot(k, Pchichi(k,As,R), label=
#    r'$P_{{\chi\chi}}(k)=\frac{{2\pi^2}}{{k^3}}A_s(kR)^2e^{{-(kR)^2}}$, $A_s={:.3f}\times10^{{-10}}$, $R=2a_{{\rm{{latt}}}}={:.2f}$'.format(As*1e10,R) )

## Add Gaussian bump
#Pchi_sliNGbump = Pchi_sliNG * (1+1.0e2*np.exp( -(k-1.0e-2)**2/(2*(1.0e-3)**2) ))
#ax2.plot(k, Pchi_sliNGbump, label=r'with bump')

ax2.legend(framealpha=0.0)

'''
# Plot P(k)
file2 = 'bbps_power_spec_NG.dat'
data2 = np.loadtxt(file2, delimiter=' ')
k2    = data2[:,0] # k
Pk2   = data2[:,1] # P_k = | \delta_k |^2
k2T2  = data2[:,2] # k^2 T, where T is the transfer function
Pchi2 = data2[:,3] # P_\chi (scaling uncertain)
T2    = k2T2*k2**-2

fig3,ax3 = plt.subplots()
ax3.plot(k, Pk, label=r'$\mathcal{P}_{Planck2018}(k)$' )   
ax3.plot(k, T, label=r'$T_{Planck2018}(k)$' ) 
ax3.plot(k2, Pk2, label=r'$\mathcal{P}_{BBPS}(k)$' )   
ax3.plot(k2, T2, label=r'$T_{BBPS}(k)$' ) 
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlabel(r'$k$')
ax3.set_ylabel(r'$\mathcal{P}(k)$')
ax3.legend()
'''
plt.show()
