#!/usr/bin/env python

import numpy as np
import sys
from scipy import stats
import scipy as sp
from scipy.special import *
from scipy import integrate
from scipy.interpolate import *
import os 
#import ntpath
import glob

verbose = True
print(len(sys.argv))
if len(sys.argv) < 5:
    print('\nusage: pksc2npz.py <Mmin [Msun/h]> <zmax (<0 for no cut)> <pattern> <fileout> <write positions [=0]> <write velocities [=0]> <write lag positions [=0]> <ni> <nf>\n')
    sys.exit(2)

Mmin = float(sys.argv[1])
zmax = float(sys.argv[2])
fpat = sys.argv[3]
fout = sys.argv[4]
ipos = int(sys.argv[5])
ivel = int(sys.argv[6])
ilag = int(sys.argv[7])
ni   = int(sys.argv[8])
nf   = int(sys.argv[9])

print( Mmin, zmax)
wpos = False
wvel = False
wlag = False
if ipos > 0 : wpos = True
if ivel > 0 : wvel = True
if ilag > 0 : wlag = True

h      = 0.7
omegam = 0.25
rho    = 2.775e11*omegam*h**2

Mmin /= h # convert to Msun
Rmin  = (3./4/np.pi/rho*Mmin)**(1./3)

z=np.linspace(0,4,1000)

def hubble(z):
    return h*100*np.sqrt(omegam*(1+z)**3+1-omegam)

def drdz(z):
    return 3e5 / hubble(z)

r  = np.cumsum(drdz(z)*(z[1]-z[0]))
r -= r[0]
z_to_r   = sp.interpolate.interp1d(z,r)
r_to_z   = sp.interpolate.interp1d(r,z)

if zmax >= 0: 
    rmax = z_to_r(zmax)
else:
    rmax = -1

print( 'rmax = ',rmax, Rmin)

RTH32M = 4*np.pi/3*rho

def get_peak_data(filename):

    outnum=11

    pkfile   = open(filename,"rb")
    Non      = np.fromfile(pkfile, dtype=np.int32, count=1)[0]
    RTHMAXin = np.fromfile(pkfile, dtype=np.float32, count=1)[0]
    redshift = np.fromfile(pkfile, dtype=np.float32, count=1)[0]
    npkdata  = outnum*Non

    if verbose: print( 'reading')
    peakdata = np.fromfile(pkfile, dtype=np.float32, count=npkdata)

    if verbose: print( 'reshaping')
    peakdata = np.reshape(peakdata,(Non,outnum))
    pkfile.close()


    xpk = peakdata[:,0]
    ypk = peakdata[:,1]
    zpk = peakdata[:,2]
    
    vxpk = peakdata[:,3]
    vypk = peakdata[:,4]
    vzpk = peakdata[:,5]

    RTH = peakdata[:,6]

    xLpk = peakdata[:,7]
    yLpk = peakdata[:,8]
    zLpk = peakdata[:,9]
    
    print( RTH.min(),RTH.max())

    # mass cut
    #if verbose: print( 'indexing RTH')
    #dm = [RTH > Rmin]

    #if verbose: print( 'cutting RTH')
    #RTH = RTH[dm]

    #if verbose: print( 'cutting positions')
    #xpk = xpk[dm]
    #ypk = ypk[dm]
    #zpk = zpk[dm]

    #xLpk = xLpk[dm]
    #yLpk = yLpk[dm]
    #zLpk = zLpk[dm]

    #if wvel:
    #    vxpk     = vxpk[dm]
    #    vypk     = vypk[dm]
    #    vzpk     = vzpk[dm]


    '''
    if rmax > 0:
        if verbose: print( 'calculating radii')
        r = np.sqrt(xpk**2+ypk**2+zpk**2)
    
        if verbose: print( 'indexing redshift')
        dm = [r<rmax]

        if verbose: print( 'cutting mass')
        RTH = RTH[dm]

        if verbose: print( 'cutting radii')
        r  = r[dm]

        if wpos:
            print( 'xpk len',len(xpk))
            xpk      = xpk[dm]
            ypk      = ypk[dm]
            zpk      = zpk[dm]
            xpkL     = xpkL[dm]
            ypkL     = ypkL[dm]
            zpkL     = zpkL[dm]
            print( 'xpk len after',len(xpk))
        if wvel:
            print( 'vxpk len',len(vxpk))
            vxpk     = vxpk[dm]
            vypk     = vypk[dm]
            vzpk     = vzpk[dm]
            print( 'vxpk len after',len(vxpk))

    if verbose: print( 'calculating mass from RTH')
    M = RTH**3 * RTH32M

    if wpos and wvel and wlag:
        return xpk,ypk,zpk,xLpk,yLpk,zLpk,vxpk,vypk,vzpk,M
    elif wpos and wvel:
        return xpk,ypk,zpk,vxpk,vypk,vzpk,M
    elif wpos:
        return xpk,ypk,zpk,M
    elif wvel:
        return vxpk,vypk,vzpk,M
    else:
        return M,r
    '''
    return xpk,ypk,zpk,xLpk,yLpk,zLpk,vxpk,vypk,vzpk,M

M = []
Nhalo = []

if wpos: 
    x = []
    y = []
    z = []
if wlag:
    xL = []
    yL = []
    zL = []
if wvel: 
    vx = []
    vy = []
    vz = []
if (not wpos) and (not wvel):
    r = []

nfiles = len(glob.glob(fpat))

ndone = 0
for file in glob.glob(fpat):

    ndone += 1
    if ndone < ni or ndone > nf: continue
    if wpos and wvel and wlag:
        xpk,ypk,zpk,xLpk,yLpk,zLpk,vxpk,vypk,vzpk,Mpk = get_peak_data(file)
    elif wpos and wvel:
        xpk,ypk,zpk,vxpk,vypk,vzpk,Mpk = get_peak_data(file)
    elif wpos:
        xpk,ypk,zpk,Mpk = get_peak_data(file)
    elif wvel:
        vxpk,vypk,vzpk,Mpk = get_peak_data(file) 
    else:
        Mpk,rpk = get_peak_data(file)

    Npk = len(Mpk)
    print( file[-5:]+" contains", Npk, "halos with M >",Mmin,"and z <",zmax,ndone,'/',nfiles)

    Nhalo.append([Npk])
    M.append(       Mpk)

    if wpos:
         x.append( xpk)
         y.append( ypk)
         z.append( zpk)
    if wlag:
         xL.append(xLpk)
         yL.append(yLpk)
         zL.append(zLpk)
    if wvel:
        vx.append(vxpk)
        vy.append(vypk)
        vz.append(vzpk)

    if (not wpos) and (not wvel):
        r.append(rpk)
    

M = np.asarray(        M)
Nhalo = np.asarray(Nhalo)

if wpos:
     x = np.asarray(   x)
     y = np.asarray(   y)
     z = np.asarray(   z)

if wlag:
     xL = np.asarray(  xL)
     yL = np.asarray(  yL)
     zL = np.asarray(  zL)

if wvel:
    vx = np.asarray(  vx)
    vy = np.asarray(  vy)
    vz = np.asarray(  vz)

if (not wpos) and (not wvel):
    r = np.asarray(    r)

if wpos and wvel and wlag:
    np.savez(fout,Nhalo=Nhalo,M=M,x=x,y=y,z=z,xL=xL,yL=yL,zL=zL,vx=vx,vy=vy,vz=vz)
elif wpos and wvel:
    np.savez(fout,Nhalo=Nhalo,M=M,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz)
elif wpos:
    np.savez(fout,Nhalo=Nhalo,M=M,x=x,y=y,z=z)
elif wvel:
    np.savez(fout,Nhalo=Nhalo,M=M,chi=r,vx=vx,vy=vy,vz=vz)
else:
    np.savez(fout,Nhalo=Nhalo,M=M,chi=r)
