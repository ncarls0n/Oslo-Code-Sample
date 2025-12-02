import numpy as np
import sys

# USAGE:
# 
# python3 setup-run.py <n> <s> [ <cluster>='Niagara'
#     <CPUs/node>=40 <rbuff>=36.0 <iLPT>=2 <cat_len>=11 ]
# 
# <n>^3  the simulation resolution (including buffers)
# <s>^3  the Lagrangian volume in Mpc of the simulation cube (excluding
#     buffers).
#
# Optional parameters:
# <cluster>    the computer to run on (Peak Patch is optimized for UofT's
#              Niagara supercomputer, so 'Niagara' is default)
# <CPUs/node>  self-explanatory (Niagara has 2048 nodes each with 40 CPUs)
# <rbuff>      the Lagrangian space thickness of the buffers in Mpc
# <iLPT>       use linear or quadratic Lagrangian perturbation theory
#              (accepts values 1 or 2)
# <cat_len>    the number of 32-bit floating point values per halo in the
#              catalogue files (11 generally, or 33 if outputting shear).

print()

if len(sys.argv)==2 and sys.argv[1] == 'help':
    print('Usage:\n\npython3 setup-run.py <n> <s> [ <cluster>=\'Niagara'+
          '\' <CPUs/node>=40\n\t<rbuff>=36.0 <iLPT>=2 <cat_len>=11 ]\n')
    sys.exit(1)

if len(sys.argv)<3:
    raise IOError('incorrect number of arguments passed at command line.\n'
        +'USAGE: <n> <s> [ <cluster> <CPUs/node> <rbuff> <iLPT> <cat_len> '
        +'] where\n<n>^3  the simulation resolution (including buffers)\n<'
        +'s>^3  the Lagrangian volume in Mpc of the simulation cube (exclu'
        +'ding\n       buffers).\nOptional parameters:\n<cluster>    the c'
        +'omputer to run on (Peak Patch is optimized for UofT\'s\n        '
        +'     Niagara supercomputer, so \'Niagara\' is default)\n<CPUs/no'
        +'de>  self-explanatory (Niagara has 2048 nodes each with 40 CPUs)'
        +'\n<rbuff>      the Lagrangian space thickness of the buffers in '
        +'Mpc\n<iLPT>       use linear or quadratic Lagrangian perturbatio'
        +'n theory\n             (accepts values 1 or 2)\n<cat_len>    the'
        +' number of 32-bit floating point values per halo in the\n       '
        +'      catalogue files (11 generally, or 33 if outputting shear).'
        )

# Size of simulation volume
n = int(sys.argv[1])   # n = (n_mesh-2n_buff)*n_tile + 2n_buff
s = float(sys.argv[2]) # a_latt * (n-2n_buff)
# n^3 is the resolution of the simulation in voxels (including buffers) and
# s^3 is the volume of the simulation in Mpc (excluding buffers)

print(n,s)

# Cluster characteristics
cluster    = 'Niagara' # cluster to run on
rbuff      = 36.       # buffer thickness in Mpc
# Niagara has 2024 nodes, our allocation allows up to 1000, 
# buffer size should be the largest peak radius you expect to find (for
# runs with boxsize < 10 Gpc/h at redshifts z < 0, we expect
# r_buff <~ 40 Mpc/h).

# Read in optional inputs
if len(sys.argv)>3:
    cluster = sys.argv[3]

if cluster.lower()=='niagara':
    nproc      = 40   # number of CPUs/node
    noderam    = 188. # RAM/node in GiB (i.e. 1 GiB = 2^10 MiB = 2^30 B)
    totalnodes = 1000 # nodes available on cluster

# Single-cpu niagara cluster
elif cluster.lower()=='niagara_single':
    nproc      = 1   # number of CPUs/node
    noderam    = 188. # RAM/node in GiB (i.e. 1 GiB = 2^10 MiB = 2^30 B)
    totalnodes = 1000 # nodes available on cluster

# Sunnyvale is a heterogeneous cluster of about 120 nodes with varying
# sizes and speeds. This will probably make it more difficult to get the
# full parallel code running on it, so I haven't done this yet.
# elif cluster.lower()=='sunnyvale':
#     nproc      = ?
#     noderam    = ?
#     totalnodes = ~120

else:
    print('I don\'t recognize cluster "{0}".'.format(sys.argv[3]))
    print('I\'m going to assume you meant "Niagara".')
    cluster = 'Niagara'
    nproc   = 40
    noderam = 188.
    totalnodes = 1000

if len(sys.argv)>4:
    nproc = int(sys.argv[4]) # number of CPUs/node
if len(sys.argv)>5:
    rbuff = float(sys.argv[5]) # buffer thickness in Mpc/h

# Number of fields stored (this determines the RAM load)
if len(sys.argv)>6:
    lpt = int(sys.argv[6])
else:
    lpt = 2
if lpt==1:
    nfields=4
elif lpt==2:
    nfields=7

rlatt  = s/n                  # lattice spacing
nbuffi = int(rbuff/rlatt+0.5) # min buffer rounded to integer # of voxels

# Function for n_mesh: the sidelength in voxels of a cubic simulation tile
# (w/ buffers) as a function of n (sidelength in voxels of full simulation
# cube w/ buffers), nb (n_buff, the thickness in voxels of the buffers) and
# nt (n_tile, the sidelength of the simulation volume in tiles)
ng = lambda n,nb,nt: (n+2*nb*(nt-1))/nt # (n-2n_buff)/n_tile + 2n_buff

# n_buff as an integer based on n_grid, n_buff, n_tile, r_buffer
nb = lambda nmesh,nbuff,ntile,rb: int( rb*(nmesh-2*nbuff)*ntile/s +0.5)

# check if number is prime. FFTW is very slow for primes so don't use these
def is_prime(a):
    return all(a % i for i in range(2, int(a/2)+1))

ngrid = []
nt    = []
nbuff = []
found = 0

# nmax^3 is roughly how many 32-bit floats that we can fit per CPU on the
# selected computer architecture, so a tiling in which the simulation
# volume is divided into cubes of sidelength n_mesh voxels with n_mesh <=
# nmax is to be found
nmax = int(         # [
    (0.5*noderam    # 1/2 RAM in GiB per node
    /nproc          # / number of CPUs per node
    *2**30          # * 2^30 B/GiB
    /4              # / 4 B per 32-b float
    /(nfields+1.25) # / ( 7 fields + 1 working field + 1/4 extra )
    )**(1./3)       # ]^{1/3}
    +0.5)           # rounded to the nearest integer
nmin = 199 # Not much point running tiles smaller than nmesh = 199

if n<nmin:
    raise ValueError( ('Potential parallelization configurations for a box'
        +' with n={ 0},\nboxsize={1} Mpc/h and buffersize={2} Mpc/h requir'
        +'e that integers nmesh,\nnbuff and ntile be found such that\n    '
        +'n = (nmesh-2nbuff)*ntile+2nbuff\nAnd for the {3} cluster {4}'
        +' < nmesh <= {5}, which means that n must be\nless than {4}.'
        ).format( n,s,rbuff,cluster,nmin,nmax ) )

# find possible values for ntile, nbuff, nmesh in order to get desired resolution
for j in range(5):
    nti = 1
    for i in range(50):
        ngridi = ng(n,nbuffi,nti) # estimate n_mesh
        if(ngridi.is_integer()==True and ngridi<nmax and
          ngridi>nmin and is_prime(int(ngridi))==False):
            if(nbuffi>nb(ngridi,nbuffi,nti,rbuff)):
                ngridi = int(ngridi)
                print("possible ntile, nbuff, nmesh = ", nti,nbuffi,ngridi)
                ngrid.append(ngridi)
                nt.append(nti)
                nbuff.append(nbuffi)
                found = 1
        nti += 1
    nbuffi += 1

if(found==0):
    raise ValueError( ('No potential parallelization configurations could '
        +'be found with n={0},\nboxsize={1} Mpc/h and buffersize={2} Mpc/h'
        +'. Parallelization requires that\nintegers nmesh, nbuff and ntile'
        +' be found such that\n    n = (nmesh-2nbuff)*ntile+2nbuff\nwhere'
        +'\n    nbuff >= (nmesh-2nbuff)*ntile*buffersize/boxsize\nand for '
        +'the {3} cluster, {4} < nmesh < {5}. Any of these conditions '
        +'not\nbeing met could cause this error.'
        ).format( n,s,rbuff,cluster,nmin,nmax ) )

# use ntile, nbuff, nmesh in order to set minimum number of nodes,
# while also evenly dividing slabs across total number processors (required by fftw)
npk = 1e6
if len(sys.argv)>7:
    catparams = int(sys.argv[7])
else:
    catparams = 11#33

# Memory of each code section:
def fftw_mem(n,nnodes):
    # Returns memory in GiB/node needed to perform FFTs on nfields (=7
    # for 2LPT) fields of dimension nxnxn 32-bit floating point values
    # spread over nnodes computational nodes.
    return(nfields # 7 fields per tile
        *4         # * number of B per 32-b float
        /2**30     # / number of B per GiB
        *n**3      # * n^3
        /nnodes)   # / nnodes

def grid_mem(ngrid,nproc):
    # Returns the memory in GiB needed to store and work on fields (usually
    # 7 for 2LPT + 1 working field + 1/4 extra space) for a serial Peak
    # Patch run with ngrid^3 32-bit floating point values times the number
    # of CPUs per node.
    return(nproc   # number of CPUs per node
        *ngrid**3  # * number of 32-b floats per parallelization "tile"
        *( nfields # *( 7 fields
           +1      #    + 1 working field
           +.25 )  #    + 1/4 field for extra space? )
        *4         # * number of B per 32-b float
        /2**30)    # / number of B per GiB

def cat_mem(npk,nproc):
    # Returns the memory in GiB needed to store catalogue files (which will
    # have maximum size npk x catparams where catparams is 33 for runs
    # where we output shear data and 11 for runs without shear) times the
    # number of CPUs per node.
    return (nproc   # number of CPUs per node
        *catparams # * number of 32-b floats per peak in catalogue file
        *npk       # * max number of peaks per catalogue file
        *4         # * number of B per 32-b float
        /2**30)    # / number of B per GiB
    # Note that in some legacy versions of Peak Patch, catparams was either
    # 11 or 23, so you may see 23 pop up randomly in some places.

def lat_mem(nbuff,nproc):
    # Returns the memory in GiB needed to perform peak-finding on lattice.
    return(nproc # number of CPUs per node
        *24      # number of parameters per voxel in peak finding
        *4/3*np.pi*nbuff**3 # * volume (in voxels) of largest peaks
        *4       # * number of B per 32-b float
        /2**30)  # / number of B per GiB

def s2c_mem(n,ngrid,nnodes,nproc):
    # Returns the memory in GiB used to do the slab-to-cube transformation
    return((# number of voxels in a tile:
            ngrid**3       # n_mesh^3
            # number of voxels in a slab:
            +n             # sidelength of full cube /w buff in voxels
            /nnodes        # / number of nodes                         
            /nproc         # / number of CPUs per node                 
            *ngrid**2)     # n_mesh^2                                         
        *nproc         # * number of CPUs per node
        *4             # * number of B per 32-b float
        /2**30)        # / number of B per GiB

# try all possible ntile, nbuff, nmesh values, 
# and see which one splits up best for cluster architecture
nnodes_tot = []
check_ram  = False # Switched to True if a configuration with memory less than noderam is found
check_proc = False # Switched to True if a configuration with an integer n/nnodes/nproci found
for i in range(len(nt)):

    nti = nt[i]       # number of tiles
    nbuffi = nbuff[i] # n_buff
    ngridi = ngrid[i] # n_mesh

    # Number of processes is whichever is less, the number of tiles
    # n_tile^3 or the number of CPUs per node nproc
    nproci = min(nti**3,nproc)

    nnodes = 1          # start with one node and count up
    min_reached = False # flag
    while not min_reached and nnodes<=totalnodes:
        fftw = fftw_mem(n,nnodes)
        grid = grid_mem(ngridi,nproci) 
        cat = cat_mem(npk,nproci)
        lat = lat_mem(nbuffi,nproci)
        s2c = s2c_mem(n,ngridi,nnodes,nproci)
        total_mem = fftw+grid+cat+lat+s2c # Total memory in GiB

        if total_mem<noderam:
            check_ram=True
        else:
            check_ram=False
        if (n/nnodes/nproci).is_integer():
            check_proc=True
        else:
            check_proc=False

        if(total_mem<noderam and check_proc):
            print("\nIf using possible values:")
            print("nnodes,ntile,nbuff,ngrid = ",nnodes,nti,nbuffi,ngridi)
            print("Memory needed for FFTW           ~ ", fftw, "GiB")
            print("Memory needed for initial fields ~ ", grid, "GiB")
            print("Memory needed for halo catalogue ~ ", cat, "GiB")
            print("Memory needed for peak finding   ~ ", lat, "GiB")
            print("Memory needed for slab2cube      ~ ", s2c, "GiB")
            print('Total memory per node            ~ ',total_mem,'GiB /',int(noderam),'GiB')
            min_reached = True
            nnodes_tot.append(nnodes)
        nnodes+=1

# If we couldn't find a configuration that works, this reports that
if not check_ram and check_proc:
    raise ValueError( ('No configuration was found for which the memory re'
        +'quired was less than the\n{0} GiB/node allowed by {1}. You lik'
        +'ely need to choose a different value\nof n with more factors so '
        +'that it will be more likely that we can find a\nresoluiton ntile'
        +' for the parallelization of the volume.'
        ).format(noderam,cluster) )

elif check_ram and not check_proc:
    raise ValueError( ('No configuration was found for which n/nnodes/npro'
        +'ci was an integer, which\nis required for the parallelization sc'
        +'heme used by FFTW. The values used\nwere\nn     ={0}\nnnodes={1}'
        +'\nnproci=min(ntile**3,nproc)\nIf n or nnodes are prime numbers, '
        +'that is likely the culprit, as the\nalgorithm tries to find ntil'
        +'e that divides the total volume evenly (note\nthat nproc for {2}'
        +' is {3}).').format(n,nnodes,cluster,nproc) )

elif not check_ram and not check_proc:
    raise ValueError( ('No configuration was found for which the memory re'
        +'quired was less than the\n{0} GiB/node allowed by {1}, or for '
        +'which n/nnodes/nproci was an integer.\nYou likely need to choose'
        +' a different value of n with more factors so that\nit will be mo'
        +'re likely that we can find a resoluiton ntile for the paral-\nle'
        +'lization of the volume.').format(noderam,cluster) )

else:
    # use setup that requires minimum number of nodes
    nnodes_min = min(nnodes_tot)
    ntile_min  = 999
    for i in range(len(nnodes_tot)):
        if nnodes_tot[i]<=nnodes_min:
            nnodes_min = nnodes_tot[i]
            if nt[i]<=ntile_min:
                ntile_min  = nt[i]
                ind = i

    nt     = nt[ind]
    nbuff  = nbuff[ind]
    ngrid  = ngrid[ind]
    nnodes = nnodes_tot[ind]

    print("\n====================================")
    print("Values optimized for number of nodes")
    print("ntile  = ",nt)
    print("nbuff  = ",nbuff)
    print("nmesh  = ",ngrid)
    print("nnodes = ",nnodes)
    print("tpnode = ",min(nt**3,nproc))
    print("====================================\n")
