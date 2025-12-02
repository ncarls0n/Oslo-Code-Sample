import numpy as np
import sys

# USAGE:
# 
# Geven a Peak Patch halo catalogue /tmp/run_merge.pksc.13579 we can create a
# .npz file by running
# 
#     python pksc3npz.py /tmp/run_merge.pksc.13579 /tmp/run_merge.npz
# 
# To read this in a python script, you would then run use
# 
#     with np.load( '/tmp/run_merge.npz' ) as data:
#         M        = data['M']
#         x, y, z  = data['x'],  data['y'],  data['z']
#         vx,vy,vx = data['vx'], data['vy'], data['vz']
# 
# etc. where M is an array with halo masses, x,y,z are arrays with halo final
# (Eulerian) positions, and vx,vy,vz are arrays of halo velocities.

# Read input/output files from command line
file_in  = open( sys.argv[1] , 'rb' )
file_out = sys.argv[2]

# Read .pksc file header
Non      = np.fromfile(file_in, dtype=np.int32, count=1)[0]
RTHmax   = np.fromfile(file_in, dtype=np.float32, count=1)[0]
redshift = np.fromfile(file_in, dtype=np.float32, count=1)[0]

# Read .pksc file catalogue
peakdata = np.fromfile(file_in, dtype=np.float32, count=int(Non*11)).reshape((Non,11))

# Eulerian position of halos
xpk = peakdata[:,0]
ypk = peakdata[:,1]
zpk = peakdata[:,2]

# Halo velocities
vxpk = peakdata[:,3]
vypk = peakdata[:,4]
vzpk = peakdata[:,5]

# Halos found at filter scalae
RTH = peakdata[:,6]

# Halo comoving distance
xLpk = peakdata[:,7]
yLpk = peakdata[:,8]
zLpk = peakdata[:,9]

# Field density at halo peak
Fpk  = peakdata[:,10]

# Halo mass
h      = .6735
omegam = .3138
rho    = 2.775e11*omegam*h**2
M      = RTH**3 *4/3 * np.pi * rho

# Save catalogue as .npz
np.savez( file_out, Nhalo=Non, M=M, x=xpk, y=ypk, z=zpk, xL=xLpk, yL=yLpk,
          zL=zLpk, vx=vxpk, vy=vypk, vz=vzpk )
