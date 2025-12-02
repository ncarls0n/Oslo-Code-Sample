import numpy as np
import sys

# Usage
#
# python3 python/find_tiling.py s_box n [n_tile=1 [R_THmax=34]]
# 
# where s_box is the sidelength of the simulation cube in Mpc/h, n is the
# total lattice size (including buffers), n_tile is the tiling (i.e. the
# simulaiton will be split into n_tile^3 equally sized boxes to be run in
# parallel), and R_THmax is the maximum filter radius.

if len(sys.argv)<3 or len(sys.argv)>5:
    print('You\'ve passed the wrong number of arguments.')
    sys.exit(2)

s_box = float(sys.argv[1])
n     = int(  sys.argv[2])

if len(sys.argv)>3:
    n_tile=int(sys.argv[3])
else:
    n_tile=1

if len(sys.argv)>4:
   R_THmax=float(sys.argv[4])
else:
   R_THmax=34.

primes = np.array([])
for i in range(2,n):
    check=False
    for j in range(2,i):
        if j<i and i%j==0:
            check=True
            continue
    if check==False:
        primes=np.concatenate((primes, np.array([i])))

def a_buff(s_box,n,n_tile,n_buff):
    return s_box*n_buff*n_tile/(n-2*n_buff)

def a_latt(s_box,n,n_buff):
    return s_box/(n-2*n_buff)

check = False
for n_buff in range(1,int(n/n_tile)+1):

    a=a_buff(s_box,n,n_tile,n_buff)   
 
    if a>=R_THmax and ((n-2*n_buff)/n_tile)%1==0 and ((n-2*n_buff)/n_tile) not in primes:
        check=True
        print('\ns_box   = {0} Mpc/h\nn       = {1}\nn_tile  = {2}\nn_buff  = {3}\nR_THmax = {4} Mpc/h\nR_buff  = {5} Mpc/h\na_latt  = {6} Mpc/h\n'.format(s_box,n,n_tile,n_buff,R_THmax,a,a_latt(s_box,n,n_buff)))
        break

if check==False:
    print('no solution found')
