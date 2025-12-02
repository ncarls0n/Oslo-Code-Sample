import numpy as np
import sys
import os

# USAGE
# This python script checks the output of a peak patch run. It can be run
# at command line by running:
#     python3 check_catalogue.py <...>/<short_name>_merge.pksc.<seed>

# Raise an error if incorrect number of arguments at command line
if len(sys.argv) != 2:
    raise ValueError('You\'ve passed me the wrong number of arguments at t'
        +'he command line. Command\nline should take the form\npython3 che'
        +'ck_catalogue.py <...>/<short_name>_merge.pksc.<seed>')

# Open catalogue file
catalogue_file = str(sys.argv[1]).strip()
in_catalogue = open( catalogue_file, 'rb' )

# Read header from catalogue file and print
num_halos = np.fromfile( in_catalogue, dtype=np.int32,   count=1 )[0]
R_THmax   = np.fromfile( in_catalogue, dtype=np.float32, count=1 )[0]
z_thing   = np.fromfile( in_catalogue, dtype=np.float32, count=1 )[0]
print('Number of halos:  ',num_halos)
print('Max filter scale: ',R_THmax)
print('zthing:           ',z_thing)

# Read file size
bytes_per_halo_param = 4 # 4 for 32-bit floating-point numbers
header_length        = 3 # number of e.g. 32-bit floats or integers
cat_file_size   = os.stat(catalogue_file).st_size # bytes
num_halo_params = int(
    (cat_file_size/bytes_per_halo_param-header_length)/num_halos )
print('This catalogue appears to have {0} parameters per halo.'
    .format(num_halo_params))

# Read catalogue body
catalogue = np.fromfile( in_catalogue, dtype=np.float32,
    count=num_halo_params*num_halos )
catalogue = np.reshape( catalogue, (num_halos,num_halo_params) )

print('The mean value for initial position is ({0},{1},{2}) where {3}<x_L^i<{4}'
    .format( str(round( np.mean(catalogue[:,7])  , 2 )) ,
             str(round( np.mean(catalogue[:,8])  , 2 )) ,
             str(round( np.mean(catalogue[:,9])  , 2 )) , 
             str(round( np.min(catalogue[:,7:9]) , 2 )) ,
             str(round( np.max(catalogue[:,7:9]) , 2 )) ))

# Print some samples of catalogue data
print('x','y','z','v_x','v_y','v_z','R_TH','x_L','y_L','z_L','F_pk',sep='\t\t')
for i in range(2):
    for j in range(11):
        if j<10:
            print(catalogue[i,j],end='\t')
        else:
            print(catalogue[i,j])

print('.\n.\n.')

for i in range(num_halos-1,num_halos):
    for j in range(11):
        if j<10:
            print(catalogue[i,j],end='\t')
        else:
            print(catalogue[i,j])

