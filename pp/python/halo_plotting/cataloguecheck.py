import numpy as np
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

#This simple script plots all the variables in a merged peak patch halo catalogue, in order to confirm the properties make sense

if len(sys.argv) != 2:
    print '\nusage: python cataloguecheck.py <merged peakfile> '
    sys.exit(2)

infile   = open(sys.argv[1],"rb")

#Cosmo parameters to get mass from R_TH
omegam  = 0.25
h       = 0.7

rho      = 2.775e11*omegam*h**2
deltavir = 200

outnum  = 11

#Load in peak patch catalogue
Non    = np.fromfile(infile,dtype=np.int32,count=1)[0]
RTHmax = np.fromfile(infile,dtype=np.float32,count=1)[0]
zin    = np.fromfile(infile,dtype=np.float32,count=1)[0]
print "\nNumber of halos to read in = ",Non

npkdata=outnum*Non
peakdata = np.fromfile(infile, dtype=np.float32, count=npkdata)
peakdata = np.reshape(peakdata,(Non,outnum))
    
x      = peakdata[:,0]
y      = peakdata[:,1]
z      = peakdata[:,2]
vx     = peakdata[:,3]
vy     = peakdata[:,4]
vz     = peakdata[:,5]
Rth    = peakdata[:,6] 
xL     = peakdata[:,7]
yL     = peakdata[:,8]
zL     = peakdata[:,9]
deltah = peakdata[:,10]

M = 4*np.pi/3*Rth**3*rho
#MAKE FIGURE
f = plt.figure(figsize=(12, 8)) 
gs = gridspec.GridSpec(4, 3) 
ax00 = f.add_subplot(gs[0,0])
ax01 = f.add_subplot(gs[0,1])
ax02 = f.add_subplot(gs[0,2])
ax10 = f.add_subplot(gs[1,0])
ax11 = f.add_subplot(gs[1,1])
ax12 = f.add_subplot(gs[1,2])
ax20 = f.add_subplot(gs[2,0])
ax21 = f.add_subplot(gs[2,1])
ax22 = f.add_subplot(gs[2,2])
ax30 = f.add_subplot(gs[3,0])
ax31 = f.add_subplot(gs[3,1])
ax32 = f.add_subplot(gs[3,2])

plt.subplots_adjust(wspace=0.4,  hspace=0.6)


ax00.plot(xL[abs(zL)<25.],yL[abs(zL)<25.],".",alpha=0.25)
ax01.plot(xL[abs(yL)<25.],zL[abs(yL)<25.],".",alpha=0.25)
ax02.plot(yL[abs(xL)<25.],zL[abs(xL)<25.],".",alpha=0.25)

ax10.plot(x[abs(z)<25.],y[abs(z)<25.],".",alpha=0.25)
ax11.plot(x[abs(y)<25.],z[abs(y)<25.],".",alpha=0.25)
ax12.plot(y[abs(x)<25.],z[abs(x)<25.],".",alpha=0.25)

ax20.plot(vx,vy,".",alpha=0.25)
ax21.plot(vx,vz,".",alpha=0.25)
ax22.plot(vy,vz,".",alpha=0.25)

ax30.hist(M,bins=np.logspace(13,15.5,100))
ax30.set_xscale('log')
ax30.set_yscale('log')
ax31.semilogx(M,deltah,".",alpha=0.25)

ax00.set_xlabel("x Lag [Mpc]")
ax01.set_xlabel("x Lag [Mpc]")
ax02.set_xlabel("y Lag [Mpc]")
ax00.set_ylabel("y Lag [Mpc]")
ax01.set_ylabel("z Lag [Mpc]")
ax02.set_ylabel("z Lag [Mpc]")

ax10.set_xlabel("x Eul [Mpc]")
ax11.set_xlabel("x Eul [Mpc]")
ax12.set_xlabel("y Eul [Mpc]")
ax10.set_ylabel("y Eul [Mpc]")
ax11.set_ylabel("z Eul [Mpc]")
ax12.set_ylabel("z Eul [Mpc]")

ax20.set_xlabel("vx [km/s]")
ax21.set_xlabel("vx [km/s]")
ax22.set_xlabel("vy [km/s]")
ax20.set_ylabel("vy [km/s]")
ax21.set_ylabel("vz [km/s]")
ax22.set_ylabel("vz [km/s]")

ax30.set_xlabel("M [Msun]")
ax30.set_ylabel("Number")
ax31.set_xlabel("M [Msun]")
ax31.set_ylabel("delta halo")

plt.show()
