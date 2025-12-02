import numpy as np
import sys

# Initial field file `Fvec_...`
file_in = open(sys.argv[1],'rb')

for i in range(10):
    print np.fromfile( file_in, dtype=np.float64, count=1)[0]

