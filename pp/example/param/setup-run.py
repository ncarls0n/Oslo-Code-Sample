import numpy as np
import sys
print " "

#USAGE: <resolution> <boxsize>
if len(sys.argv)<3:
    sys.exit("USAGE: <ngrid> <boxsize> [cluster:<(niagara, gpc)> nproc per node:<(32, 8)> <rbuff>]\n")

n       = float(sys.argv[1])
bx      = float(sys.argv[2])

cluster = 'niagara'
rbuff   = 36.

if len(sys.argv)>3:
    cluster = sys.argv[3]

if cluster=='niagara':
    nproc   = 32
    noderam = 185.
elif cluster=='gpc':
    nproc   = 8
    noderam = 12.
else:
    print "don't recognize cluster"+sys.argv[3]+", assuming niagara"
    cluster = 'niagara'
    nproc   = 32
    noderam = 185.

if len(sys.argv)>4:
    nproc = float(sys.argv[4]) #number of processors per node
if len(sys.argv)>5:
    rbuff = float(sys.argv[5])


#Number of fields need stored
lpt = 2
if lpt==1:
    nfields=4
if lpt==2:
    nfields=7

rlatt = bx/n
nbuffi = int(rbuff/rlatt+0.5)

ng = lambda n,nb,nt: (n+2*nb*(nt-1))/nt
nb = lambda ngrid,nbuff,ntile,rb: int(rb*(ngrid-2*nbuff)*ntile/bx+0.5)
def is_prime(a): # check if number is prime. FFTW is very slow for primes so don't use these
    return all(a % i for i in xrange(2, a))

ngrid = []
nt = []
nbuff = []
found = 0

# set nmax tile grid memory is maximum half the total available per node
nmax = int( (noderam/2 * 1./nproc * 1024.**3/(4*8.25))**(1./3) + 0.5)
nmin = 255

# find possible values for ntile, nbuff, nmesh in order to get desired resolution
for j in range(5):
    nti = 1
    for i in range(50):
        ngridi = ng(n,nbuffi,nti)
        if( (ngridi.is_integer() == True) and ( ngridi < nmax) and (ngridi > nmin) and (is_prime(int(ngridi))==False) ):
            
            if(nbuffi>nb(ngridi,nbuffi,nti,rbuff)):
                ngridi = int(ngridi)
                print "possible ntile, nbuff, nmesh = ", nti,nbuffi,ngridi
                ngrid.append(ngridi)
                nt.append(nti)
                nbuff.append(nbuffi)
                found = 1
        nti += 1
    nbuffi += 1

if(found==0):
    print("ERROR: too high resolution for "+str(nproc)+" tasks per node")
    sys.exit()


# use ntile, nbuff, nmesh in order to set minimum number of nodes,
# while also evenly dividing slabs across total number processors (required by fftw)
npartmax = 40.
npk      = 1e6

# to get memory of each code section
fftw_mem = lambda n,nnodes: nfields*4*n**3/1024**3/nnodes
grid_mem = lambda ngrid,nproc: nproc*ngrid**3*(nfields+1+1./4)*4./1024**3
cat_mem  = lambda npk,nproc: nproc*92.*npk/1024**3
lat_mem = lambda nbuff,npartmax,nproc: nproc*96.*4/3*np.pi*nbuff**3/1024**3
s2c_mem = lambda n,ngrid,nnodes,nproc : (ngrid**3 + n/nnodes/nproc*ngrid**2)*nproc*4./1024**3 
                            #(recieve_buffer + send_buffer) 


# try all possible ntile, nbuff, nmesh values, 
# and see which one splits up best for cluster architecture

nnodes_tot = []
for i in range(len(nt)):

    nti = nt[i]
    nbuffi = nbuff[i]
    ngridi = ngrid[i]

    nproci = min(nti**3, nproc)

    print "\nIf using possible values:"

    nnodes = 1 
    min_reached = False
    while not min_reached and nnodes<64:

        fftw = fftw_mem(n,nnodes)
        grid = grid_mem(ngridi,nproci) 
        cat = cat_mem(npk,nproci)
        lat = lat_mem(nbuffi,npartmax,nproci)
        s2c = s2c_mem(n,ngridi,nnodes,nproci)

        total_mem = fftw + grid + cat + lat + s2c
#        print fftw, grid, cat, lat, s2c
#        print nnodes, n, nproci, ngridi, n/nnodes/nproci, total_mem
        if( (total_mem < noderam) and ( (n/nnodes/nproci).is_integer() ) ):
            print "\nnnodes, ntile, nbuff, ngrid = ",nnodes, nti, nbuffi, ngridi
            print "fftw_mem,", fftw, "Gb" 
            print "grid_mem,", grid, "Gb"
            print "cat_mem,", cat, "Gb"
            print "lat_mem,", lat, "Gb"
            print "s2c_mem", s2c, "Gb"
            print "\ntotal memory per node ~ ",total_mem, "Gb"
            min_reached = True
            nnodes_tot.append(nnodes)
        nnodes += 1


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

print "\n==============================="
print "RECOMMENDED VALUES"
print "ntile = ",nt
print "nbuff = ",nbuff
print "nmesh = ",ngrid
print "nnodes = ",nnodes
print "===============================\n"
