import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys

run_dir = '/mnt/scratch-lustre/njcarlson/peak-patch-runs/23.01.30_redshift_kernels'
run_name= '1000Mpc_n580_nb40_nt2_'
model   = sys.argv[1]

#Set colourbar labels 
if model=='tsz': ylabel = r'$\frac{dy}{dz}$'#r'$d\log_{10}(y)/dz$'
if model=='ksz': ylabel = r'$\frac{d|\Delta T_{\rm{kSZ}}|}{dz} \ [\mu\rm{K}]$'
if model=='kap': ylabel = r'$\frac{d\kappa}{dz}$'
if model=='cib': ylabel = r'$\frac{dI}{dz}$ [MJy/sr]'#r'$d\log_{10}I/dz$ log MJy/sr'

seed = np.array([ 13579,23456,31317,42297,59841,66600,76543,81111 ]) #,92154,10314 ])
cenz = np.array([   500, 1500, 2500, 3500, 4500, 5500, 6500, 7500 ]) #, 8500, 9500 ])
distance_bins = np.concatenate( (np.array([0.]),cenz+500.) )
z_bins = np.array([ .0,.23854,.51414,.84518,1.2587,1.7960,2.5237,3.5553,5.0984 ]) #,7.5685,11.901 ])

# Set font sizes etc.
plt.rcParams.update({'axes.labelsize': 36, 'font.size': 36#,
                    #'xtick.labelsize': 36,
                    #'ytick.labelsize': 36, 'axes.linewidth': 4, 
                    #'xtick.labelsize': 36
                    })
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('font', **font)
#change font
plt.rc('text', usetex=True)
#plt.rc('font',**{'family':'serif','serif':['Palatino']})

#Set up figure
fig = plt.figure(figsize=(12,12))


ax  = plt.subplot(111)
ax.set_xlabel(r'comoving distance $\chi$ [Mpc]')
ax.set_ylabel(ylabel)
ax.set_yscale('log')

#for axis in ['top','bottom','left','right']:
#    ax.spines[axis].set_linewidth(4)
#ax.tick_params('both', length=10, width=4, which='major')
#ax.tick_params('both', length=6, width=3, which='minor')

dI_dz = np.zeros(8)

M_halo_max = np.zeros(8)

for j in range(len(seed)):

    # Read in data
    filein = '{0}/23.01.29_{1}{2}_cenz{3}Mpc/maps/{1}{4}_{2}_fs.map'.format(run_dir,run_name,seed[j],cenz[j],model)
    read_file=open(filein,"rb")

    # Read in 16 byte header (2 ints, 2 floats)
    npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
    npix = np.fromfile(read_file,dtype=np.int32,count=1)[0]
    fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
    fov  = np.fromfile(read_file,dtype=np.float32,count=1)#[0]
    #fov  = float(fov/2/np.pi*360) #convert from rad to deg
    
    # Fraction of sky covered
    fsky = 4*fov*np.sin(fov) / (4*np.pi)

    # Read in map
    fsmap = np.fromfile(read_file,dtype=np.float32,count=npix**2)
    #fsmap = np.reshape(fsmap, (npix,npix))
    #if (model=='tsz' or model=='cib'): fsmap = np.log10(fsmap+1e-9)
    if model=='ksz': fsmap = np.abs(fsmap)

    # dI/dz averaged over box, extrapolated to full sky
    dI_dz[j] = fsky * np.mean(fsmap) / (z_bins[j+1]-z_bins[j])


    #in_catalogue=open('{0}/23.01.29_{1}{2}_cenz{3}Mpc/output/{1}merge.pksc.{2}'.format(run_dir,run_name,seed[j],cenz[j]))
    #Non     = np.fromfile(in_catalogue,dtype=np.int32,count=1)[0]
    #RTHmax    = np.fromfile(in_catalogue,dtype=np.float32,count=1)[0]
    #zin       = np.fromfile(in_catalogue,dtype=np.float32,count=1)[0]
    #catalogue = np.fromfile(in_catalogue,dtype=np.float32,count=11*Non)
    #catalogue = np.reshape(catalogue,(Non,11))
    #Rth    = catalogue[:,6]
    #del(catalogue)
    #rho=2.775e11*(0.6735)**2 * (0.2645+0.0493)
    #M = 4*np.pi/3*Rth**3*rho
    #M_halo_max[j] = np.max(M)
    #del(M)

ax.hist(distance_bins[:-1], distance_bins, weights=dI_dz)
ax.set_xlim([distance_bins[0],distance_bins[-1]])

print(dI_dz)

ax2 = ax.twiny()
ax2.set_xticks( ax.get_xticks() )
ax2.set_xticklabels(z_bins)
ax2.set_xlabel(r'redshift $z$')

#Save figure
fileout = str(model)+'_redshift_kernel.png'
plt.savefig(fileout, bbox_inches="tight")
#plt.show()




#fig_,ax_ = plt.subplots(figsize=(12,12))
#ax_.plot(cenz,M_halo_max,ls='none',marker='.')
#ax_.set_xlabel('comoving distance [Mpc]')
#ax_.set_ylabel('max halo mass in this cube [Msol]')
#ax_.set_yscale('log')
#plt.savefig('thing.png',bbox_inches="tight")

