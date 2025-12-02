import numpy as np
import sys
from scipy import stats
import scipy
from scipy.special import *
from scipy import integrate
from scipy.interpolate import *
import os 
import ntpath
import glob

# input should be of the form "/path/to/outputdir/<basename>". eg /output/1750Mpc_2048
dirin  = ntpath.dirname(sys.argv[1])+"/"
basein = ntpath.basename(sys.argv[1])
print "\ndirin    = ",dirin
print "basename = ",basein

h      = 0.7
omegam = 0.25
rho    = 2.775e11*omegam*h**2

L    = 1750.
n    = 2000
Mmin = 100*(L/n)**3*rho
Rmin = (3./4/np.pi/rho*Mmin)**(1./3)
print "Cutting at Mmin, Rmin = ", Mmin, Rmin
print ""


#DATA NAMES:
#d2F is the laplacian at the final peak size
#F   is the overdensity at the final peak size
#Rf  is the filtersize it was found
#gradRth is the mean gradient at the final peak size
#gradF is the mean gradient of a random location in the field at the final peak size
def get_peak_data(filename):

    outnum=23

    pkfile   = open(filename,"rb")
    Non      = np.fromfile(pkfile, dtype=np.int32, count=1)
    RTHMAXin = np.fromfile(pkfile, dtype=np.float32, count=1)
    redshift = np.fromfile(pkfile, dtype=np.float32, count=1)
    npkdata  = outnum*Non
    peakdata = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    peakdata = np.reshape(peakdata,(Non,outnum))
    pkfile.close()
    # mass cut
    dm = [(peakdata[:,6] > Rmin)]
    peakdata = peakdata[dm]

#    x        = peakdata[:,0]
#    y        = peakdata[:,1]
#    z        = peakdata[:,2]
#    vx       = peakdata[:,3]
#    vy       = peakdata[:,4]
#    vz       = peakdata[:,5]
    Rth      = peakdata[:,6]
    d2F      = peakdata[:,7]
    F        = peakdata[:,8]
    ev       = peakdata[:,9]
    pv       = peakdata[:,10]
    Rf       = peakdata[:,11]
#    FRf      = peakdata[:,12]
#    d2FRf    = peakdata[:,13]
    gradRth   = np.sqrt(peakdata[:,14]**2+peakdata[:,15]**2+peakdata[:,16]**2)
    gradF    = np.sqrt(peakdata[:,20]**2+peakdata[:,21]**2+peakdata[:,22]**2)


#    return x,y,z,vx,vy,vz,Rth,d2F,F,ev,pv,Rf,FRf,d2FRf,gradRth,gradF
    return Rth,d2F,F,ev,pv,Rf,gradRth,gradF #Rf,FRf,d2FRf,gradRth,gradF


# SIGMA DATA FROM GRID                                                                                                      
pkfile=open("../sigrho_out.13579","rb") #DATAFILE WILL BE CALLED sigrho_<seed>.dat                                                    
nic = int(np.fromfile(pkfile, dtype=np.int32, count=1))
data = np.fromfile(pkfile, dtype=np.float32, count=nic*4)
data = data.reshape(nic,4)

Rthsig = data[:,0]
sig0   = data[:,1]
sig1   = data[:,2]
sig2   = data[:,3]

sig0f  = interp1d(Rthsig,sig0,bounds_error=False)
sig1f  = interp1d(Rthsig,sig1,bounds_error=False)
sig2f  = interp1d(Rthsig,sig2,bounds_error=False)

cstart = 0 
cend   = 0
count  = 0

flist_tot = glob.glob(dirin+basein+"*merge*")
ncut      = 50  #max number of merge catalogues in each .npz file
nout      = int(len(flist_tot))/ncut+1 #number of npz files

if len(flist_tot)<ncut:
    ncut = len(flist_tot)
    nout = 1
    
print "Number of files in, number out = ",len(flist_tot), nout

for i in range(nout):
    flist = flist_tot[ncut*i:ncut*(i+1)]
    count = 0
    for f in flist:
        print f
        if count==0:
#            x,y,z,vx,vy,vz,Rth,d2F,F,ev,pv,Rf,FRf,d2FRf,gradRth,gradF = get_peak_data(dirin+f)
            Rth,d2F,F,ev,pv,Rf,gradRth,gradF = get_peak_data(dirin+f)
            sigma0  = sig0f(Rth)
            sigma1  = sig1f(Rth)
            sigma2  = sig2f(Rth)
            sigma0f = sig0f(Rf)
            sigma1f = sig1f(Rf)
            sigma2f = sig2f(Rf)
        else:
#            xi,yi,zi,vxi,vyi,vzi,Rthi,d2Fi,Fi,evi,pvi,Rfi,FRfi,d2FRfi,gradRthi,gradFi = get_peak_data(dirin+f)
            Rthi,d2Fi,Fi,evi,pvi,Rfi,gradRthi,gradFi = get_peak_data(dirin+f)
            Rth = np.append(Rth,Rthi)
            d2F = np.append(d2F,d2Fi)
            F = np.append(F,Fi)
            ev = np.append(ev,evi)
            pv = np.append(pv,pvi)
            Rf = np.append(Rf,Rfi)
            gradRth = np.append(gradRth,gradRthi)
            gradF = np.append(gradF,gradFi)
            sigma0  = np.append(sigma0,sig0f(Rthi))
            sigma1  = np.append(sigma1,sig1f(Rthi))
            sigma2  = np.append(sigma2,sig2f(Rthi))
            sigma0f = np.append(sigma0f,sig0f(Rfi))
            sigma1f = np.append(sigma1f,sig1f(Rfi))
            sigma2f = np.append(sigma2f,sig2f(Rfi))


        print "finished file ", count
        count+=1


    print ""
    print "min max mean Rf     = " , np.min(Rf), np.max(Rf), np.mean(Rf)
    print "min max mean Rth    = " , np.min(Rth), np.max(Rth), np.mean(Rth)
    print ""
    print "min max mean F      = " , np.min(F), np.max(F), np.mean(F)
    print ""
    print "min max mean d2F    = " , np.min(d2F), np.max(d2F), np.mean(d2F)
    print ""
    print "min max mean gradRf = " , np.min(gradRth), np.max(gradRth), np.mean(gradRth)
    print "min max mean gradF  = " , np.min(gradF), np.max(gradF), np.mean(gradF)
    print ""

    np.savez(basein+"_append_"+str(i),Rth=Rth,d2F=d2F,F=F,
             ev=ev,pv=pv,Rf=Rf,gradRth=gradRth,
             gradF=gradF,sigma0=sigma0,sigma1=sigma1,sigma2=sigma2,
             sigma0f=sigma0f,sigma1f=sigma1f,sigma2f=sigma2f)            


#print ""
#print "min max mean x      = " , np.min(x), np.max(x), np.mean(x)
#print "min max mean y      = " , np.min(y), np.max(y), np.mean(y)
#print "min max mean z      = " , np.min(z), np.max(z), np.mean(z)
#print "min max mean vx     = " , np.min(vx), np.max(vx), np.mean(vx)
#print "min max mean vy     = " , np.min(vy), np.max(vy), np.mean(vy)
#print "min max mean vz     = " , np.min(vz), np.max(vz), np.mean(vz)

