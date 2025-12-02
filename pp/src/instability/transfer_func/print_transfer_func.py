import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 16})

# Load in the Delta phi(H_e) to Delta zeta(H_f) transfer function
data = np.loadtxt( 'transfer_func_log_spaced_k.out' )

# Find where T(k) dips below zero
flag = False
j    = -1
while flag == False:
    j += 1
    if data[j,1] <= 0:
        flag = True

# Set up a figure
fig,ax = plt.subplots()

# Plot the transfer function from the file
ax.plot( data[:,0], data[:,1], label=r'$T(k)$ read from file' )

# Plot the transfer function that's actually used
ax.plot( data[:,0], data[:,1] * np.heaviside( data[j,0] - data[:,0], 0.0 ), label=r'$T(k)$ used by Peak Patch', ls='--' )

# Labelling the plot
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel(r'$k/aH$')
ax.set_ylabel(r'$T_{\Delta\phi(H_e) \to \Delta\zeta_{ng}(H_f)}(k)$')
plt.legend()
plt.show()
