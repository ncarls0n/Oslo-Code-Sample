import numpy as np
import pylab as py
import sys
import matplotlib as m

if len(sys.argv) != 8:
    print '\nusage: python patchplot.py <peakfile> <boxsize [Mpc]> <slicewidth [Mpc]> <outfile> <z> <rankflag> <iwrap>\n'
    sys.exit(2)

infile   = open(sys.argv[1],"rb")
bx       = float(sys.argv[2])
slic     = float(sys.argv[3])
outfile  = sys.argv[4]
rankflag = int(sys.argv[6])
iwrap    = int(sys.argv[7])

if iwrap==1: wrap=True 
if iwrap==0: wrap=False 

L = sys.argv[2]
W = sys.argv[3]
Z = sys.argv[5]

simtype='Peak Patch Halos'

#Cosmo parameters to get mass from R_TH
omegam  = 0.32
h       = 0.67

rho      = 2.775e11*omegam*h**2
deltavir = 200

minmass = 1.e10
outnum  = 11

#Load in peak patch catalogue
Non    = np.fromfile(infile,dtype=np.int32,count=1)[0]
RTHmax = np.fromfile(infile,dtype=np.float32,count=1)[0]
zin    = np.fromfile(infile,dtype=np.float32,count=1)[0]
print "\nNumber of halos to read in = ",Non

npkdata=outnum*Non
peakdata = np.fromfile(infile, dtype=np.float32, count=npkdata)
peakdata = np.reshape(peakdata,(Non,outnum))
    
x    = peakdata[:,0]+bx/2
y    = peakdata[:,1]+bx/2
z    = peakdata[:,2]+bx/2
vx   = peakdata[:,3]
vy   = peakdata[:,4]
vz   = peakdata[:,5]
Rth  = peakdata[:,6] 
#xL   = peakdata[:,7]+bx/2
#yL   = peakdata[:,8]+bx/2
#zL   = peakdata[:,9]+bx/2
#Rth is  ~R_{TH,lagrangian} of the halo. Where M_200 = 4*pi/3 * Rth**3 * rho 
#Therefore, R_vir = Rth/(deltavir)**(1./3)

print "\nxmin, xmax:     ",x.min(),x.max()
print "ymin, ymax:     ",y.min(),y.max()
print "zmin, zmax:     ",z.min(),z.max()
print "RTHmin, RTHmax: ",Rth.min(),Rth.max()

#take thin z slice centered at box center
dm=[(abs(z-bx/2)<slic/2)]
x=x[dm]
y=y[dm]
z=z[dm]
vx=vx[dm]
vy=vy[dm]
vz=vz[dm]
Rth=Rth[dm]

#take x,y from 0 to bx
if wrap:
    x[x<0] += bx
    y[y<0] += bx
    z[z<0] += bx
    
    x[x>bx] -= bx
    y[y>bx] -= bx
    z[z>bx] -= bx


dm=[ (x>0.) & (x<bx) & (y>0.) & (y<bx)]
x=x[dm]
y=y[dm]
z=z[dm]
vx=vx[dm]
vy=vy[dm]
vz=vz[dm]
Rth=Rth[dm]


#Cut small halos below minmass
M   = 4*np.pi/3.*Rth**3*rho
dm  = [M>minmass]
x   = x[dm]
y   = y[dm]
vx  = vx[dm]
vy  = vy[dm]
Rth = Rth[dm]

#Get Rvir to plot from Rth
Fpk  = deltavir**(1.0/3)
Rvir  = Rth/Fpk

#keep only a certain number of halos
if rankflag<0:
    nhalos = -rankflag
    dm = np.argsort(Rvir)[-nhalos:]
    x=x[dm]
    y=y[dm]
    Rvir=Rvir[dm]
    M=M[dm]

print "\nAfter mass and position cuts:"
print "Number of halos remaining = ",len(x)
print 'min, max mass:',np.min(M),np.max(M)
print 'min, max    x:',np.min(x),np.max(x)
print 'min, max    x:',np.min(y),np.max(y)

#Multiply Rvir by factor to make more visually aesthetic
#Rvir *= 2

fb=0.0

#Plot
fig=py.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1)

print "\nPlotting %d halos\n" % len(Rth)
for i in range(len(Rvir)-1):
    ax.add_patch(py.Circle((x[i],y[i]),Rvir[i], edgecolor='b',facecolor='none',lw=1))

#titletext=simtype+"\n"+L+" x "+L+" x "+W+" Mpc/h\n"+"z = "+Z
titletext=L+" x "+L+" x "+W+" Mpc\n"+"z = "+Z
font = {'horizontalalignment' :'left',
        'fontsize'            : 26,
        'multialignment'      :'left'}

ax.set_xlabel("x [Mpc]",fontsize=26)
ax.set_ylabel("y [Mpc]",fontsize=26)
ax.tick_params(axis='both',which='major',labelsize=24)
t=py.title(titletext,font)
t.set_x(0.0)
py.xlim(bx*fb,bx*(1-fb))
py.ylim(bx*fb,bx*(1-fb))

py.savefig(outfile+'.pdf')
