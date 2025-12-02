import numpy as np

def get1LPTparams(merged_peak_file):
    '''
    '''
    # Pointer to merged, 1LPT peak patch catalogue
    in_1LPT = open( merged_peak_file, 'rb' )
    
    # Reading data from merged peak patch catalogue pointer
    outnum_1LPT    = 11
    Non_1LPT       = np.fromfile( in_1LPT, dtype=np.int32, count=1 )[0]
    RTHmax_1LPT    = np.fromfile( in_1LPT, dtype=np.float32, count=1 )[0]
    zin_1LPT       = np.fromfile( in_1LPT, dtype=np.float32, count=1 )[0]
    catalogue_1LPT = np.fromfile( in_1LPT, dtype=np.float32,
                                  count = outnum_1LPT * Non_1LPT )
    catalogue_1LPT = np.reshape( catalogue_1LPT,  (Non_1LPT,outnum_1LPT) )
    
    x  = catalogue_1LPT[:,0]  # x,y,z: components of the final
    y  = catalogue_1LPT[:,1]  #        (Eulerian) halo position
    z  = catalogue_1LPT[:,2]  #        vector [h^-1 Mpc]
    vx = catalogue_1LPT[:,3]  # vx,vy,vz: components of the final
    vy = catalogue_1LPT[:,4]  #           (Eulerian) halo velocity
    vz = catalogue_1LPT[:,5]  #           vector [km/s/Mpc]

    return x,y,z,vx,vy,vz
