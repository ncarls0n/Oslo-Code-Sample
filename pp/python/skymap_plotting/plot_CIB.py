import healpy as hp
import numpy as np
import matplotlib

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap
import ntpath,sys

cib_file = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.02.05_SBsuite/ng7_m40_test/cenz6500Mpc/cib100.fits'


#Code for plotting fancy healpy and flatsky maps
if len(sys.argv)<4:
    sys.exit("\n\n\tUSAGE: python mapview.py <model> <filein> <flatsky> [models: tsz, ksz, kap, cib, tco, thi, cmb  flatsky: 0=hp.mollview, 1=flatsky imshow, 2=hp.cartview]\n\n")

model   = sys.argv[1]
filein  = sys.argv[2]
flatsky = int(sys.argv[3])

#Set colours and max/mins
if model=='tsz': vmin, vmax, cmap = -7.8 , -3.51 , cm.get_cmap('magma').copy()
if model=='ksz': vmin, vmax, cmap = -3.99, 3.99 , cm.get_cmap('coolwarm').copy()
if model=='kap': vmin, vmax, cmap = 0.012, 0.072, cm.get_cmap('jet').copy() 
#if model=='cib': vmin, vmax, cmap = -3.3 , -2.2 , cm.get_cmap('hot').copy()
if model=='cib': vmin, vmax, cmap = -6. , -1. , cm.get_cmap('hot').copy()
if model=='cmb': vmin, vmax, cmap = -560 ,  560 , ListedColormap(np.loadtxt("./Planck_Parchment_RGB.txt")/255.)

rgbab    = 'w'
rgbaf    = 'k'
#cmap.set_under("w")

#Set colourbar labels 
if model=='tsz': cblab = 'log Compton-y'
if model=='ksz': cblab = r'$\Delta T_{\textrm kSZ}\ [\mu K]$'
if model=='kap': cblab = r'$\kappa$'
if model=='cib': cblab = 'log MJy/sr'
if model=='cmb': cblab = r'$\mu K$'


#Set font sizes etc.
plt.rcParams.update({'axes.labelsize': 36, 'font.size': 36,
                    'xtick.labelsize': 36,
                    'ytick.labelsize': 36, 'axes.linewidth': 4, 
                    'xtick.labelsize': 36})

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('font', **font)


#change font
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Palatino']})

#=====================HEALPIX mollview or cartview=====================================
if (flatsky == 0 or flatsky==2):
    if flatsky==2: cart=True  #hp.cartview for healpix map   
    if flatsky==0: cart=False #hp.mollview for healpix map      
    #Set up figure size
    dpi = 300
    if not cart: figsize_inch = 24, 16
    if cart:     figsize_inch = 12, 13

    fig = plt.figure(figsize=figsize_inch,dpi=dpi)

    #Read map in
    map = hp.read_map(filein)
    if (model=='tsz' or model=='cib'): map = np.log10(map+1e-9)
    if (model=='ksz')                : map = map * 2.725e6
    if (model=='kap')                : map = hp.smoothing(map, fwhm = np.radians(1./60*15)) # 15 arcmin smoothing seems to look good
    if (model=='cib')                : map = hp.smoothing(map, fwhm = np.radians(1./60*15)) # 15 arcmin smoothing
    if (model=='cmb')                : map = map * 1e6
    print("Map min, max, mean = ", np.min(map), np.max(map), np.mean(map))

    if not cart: 
        hp.mollview(map,fig=fig.number, xsize=figsize_inch[0]*dpi,cmap=cmap,cbar=None,title="", min=vmin,max=vmax)
    if cart:
        fov=8.1821*60 # angular field of view in arcmin
        npix=figsize_inch[0]*dpi
        ang_res=fov/npix
        hp.gnomview(map,fig=fig.number, xsize=npix,cmap=cmap,cbar=None,title="", min=vmin,max=vmax, rot=(0,90,0), reso=ang_res )
        #hp.cartview(map,fig=fig.number, xsize=figsize_inch[0]*dpi,cmap=cmap,cbar=None,title="", min=vmin,max=vmax, lonra=[-180,180],latra=[80,90])

    #Add colorbar    
    ax = plt.gca()
    image = ax.get_images()[0]
    if not cart: cbaxes=fig.add_axes([0.25,0.1,0.5,0.03])
    if cart:     cbaxes=fig.add_axes([0.1,0.05,0.8,0.04])
    cb = plt.colorbar(image,cax=cbaxes,orientation="horizontal")#, fraction=0.05)
    cb.set_label(cblab)
  
    #Save figure
    if not cart: 
        fileout = 'Fullsky_'+str(model)+'_'+ntpath.basename(filein)[:-5]+'.pdf'#'Images/Fullsky_'+str(model)+'_'+ntpath.basename(filein)[:-5]+'.pdf'
        plt.savefig(fileout,bbox_inches='tight',pad_inches=0)
    if cart:     
        fileout = 'Cartview_'+str(model)+'_'+ntpath.basename(filein)[:-5]+'.pdf'#'Images/Cartview_'+str(model)+'_'+ntpath.basename(filein)[:-5]+'.pdf'
        plt.savefig(fileout,bbox_inches='tight')




#=====================FLATSKY=====================================
if flatsky == 1:
    #Set up figure
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)

    ax.set_xlabel("Degrees")
    ax.set_ylabel("Degrees")

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    ax.tick_params('both', length=10, width=4, which='major')
    ax.tick_params('both', length=6, width=3, which='minor')

    #Read in data
    read_file=open(filein,"rb")
    #Read in 16 byte header (2 ints, 2 floats)
    npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
    npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
    fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
    fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
    fov  = float(fov/2/np.pi*360) #convert from rad to deg
    print("npix, fov = ", npix, fov)
    
    #Read in map
    map = np.fromfile(read_file,dtype=np.float32,count=npix**2)
    map = np.reshape(map, (npix,npix))

    if (model=='tsz' or model=='cib'): map = np.log10(map+1e-9)

    print("Map min, max, mean = ", np.min(map), np.max(map), np.mean(map))

    #Imshow map
    im = ax.imshow(map, vmin=vmin,vmax=vmax,cmap=cmap,extent=[-fov/2,fov/2,-fov/2,fov/2],origin="upper")
    #Add colorbar
    cbaxes=fig.add_axes([0.925,0.1,0.05,0.8])
    cb = plt.colorbar(im, cax=cbaxes)
    cb.set_label(cblab)

    #Save figure
    fileout = 'Flatsky_'+str(model)+'_'+ntpath.basename(filein)[:-6]+'fov_'+str(fov)+'.pdf'
    #fileout = 'Images/Fullsky_'+str(model)+'_'+ntpath.basename(filein)[:-5]+'_fov_'+str(fov)+'.pdf'
    plt.savefig(fileout, bbox_inches="tight")
    plt.show()


