import numpy as np
import matplotlib.pyplot as plt

# Load in Peak Patch power spectrum and transfer function
filename       = 'planck18_intermittent.dat'
data1          = np.loadtxt(filename, delimiter=' ')
k              = data1[:,0]
P_ff_peakpatch = data1[:,1]
P_ff           = (2*np.pi)**3 * P_ff_peakpatch
T_z2f_peakpatch= data1[:,2]
T_z2f          = (2*np.pi)**1.5 * T_z2f_peakpatch

# Typical linear matter power
data2 = np.array([
    [ 0.00010597966732950987, 447.6349881823656   ],
    [ 0.0001900779549468776,  789.4608456694704   ],
    [ 0.00037004049919527255, 1494.6447730252912  ],
    [ 0.000705772428652144,   2750.5848799739883  ],
    [ 0.0012529185772987977,  4715.321331600419   ],
    [ 0.0025938603320724645,  8927.270326414578   ],
    [ 0.004113801417068824,   12908.713346520148  ],
    [ 0.006867445985095601,   18143.743574961307  ],
    [ 0.010670626309758322,   22445.504092346997  ],
    [ 0.016580001661406787,   24095.18054446415   ],
    [ 0.024225513813452774,   21817.703472887944  ],
    [ 0.03162277660168379,    18665.826637911658  ],
    [ 0.04127879423736665,    14666.430727487737  ],
    [ 0.07873039552182777,    8198.939848152986   ],
    [ 0.10489903991675777,    5281.9295923507825  ],
    [ 0.30770217338405625,    823.7797050262151   ],
    [ 0.34442287541866695,    665.8994010007127   ],
    [ 0.4732320871580583,     351.72337516743244  ],
    [ 0.5406767386313648,     268.6322260242656   ],
    [ 0.7660780113838219,     123.12580302417484  ],
    [ 1.2529185772987976,     40.72454706622886   ],
    [ 2.4391603746171002,     8.08346451044848    ],
    [ 5.369951044068961,      1.1415498111326394  ]
    ])
h=0.6735
k_lmp    = data2[:,0] * h
P_ff_lmp = data2[:,1] / h**3

# Compare P(k)
fig1,ax1 = plt.subplots()
#ax1.plot( k,     P_ff,     label=filename, ls='-'  )
#ax1.plot( k_lmp, P_ff_lmp, label='sketchy', ls=':' )

# Primordial power spectrum
A_s  = 2.100e-9
n_s  = 0.9649
k_o  = 0.05
P_zz = (2*np.pi**2) * k**-3 * A_s * (k/k_o)**(n_s-1)
ax1.plot( k, k**3/(2*np.pi**2)*P_zz, label=r'$P_{\zeta\zeta}(k)$', ls='-' )

Pzz = T_z2f**-2 * P_ff

ax1.plot( k , k**3/(2*np.pi**2)*Pzz, label='Tk', ls='--' )

# Labels etc.
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel(r'$k ~ [\mathrm{Mpc}^{-1}]$')
ax1.set_ylabel(r'$P(k) ~ [\mathrm{Mpc}^3]$')
#ax1.set_xlim([1e-5,1e3])
#ax1.set_ylim([1e-6,1e6])
ax1.legend()

fig,ax = plt.subplots()
ax.plot( k , Pzz / P_zz )
print( np.mean( Pzz / P_zz ) )

plt.show()
