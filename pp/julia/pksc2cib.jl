# pksc2cib.jl
#
# This julia script generates CIB maps
# 
# USAGE:
#
#     cd <run_dir>
#     julia /home/njcarlson/Repositories/peakpatch/julia/pksc2cib.jl output/<catalogue>.pksc.<seed>

using XGPaint
using Healpix
using HDF5

# Peak Patch directory
pp_dir = "/home/njcarlson/Repositories/peakpatch/"

# Python script to make HDF5 dataset from Peak-Patch formatted catalogue file (.pksc)
pksc2h5_script = string(pp_dir,"python/catalogue_tools/pksc2hdf5.py")

# Read halo catalogue file from command line
in_halos = ARGS[1]

# If passed a Peak-Patch formatted halo catalogue
if occursin( ".pksc" , in_halos )

    # Get halo catalogue directory
    catalogue_dir  = dirname(in_halos)
    if catalogue_dir == ""
        catalogue_dir = "."
    end

    # Run pksc2hdf5.py to make HDF5 halo catalogue
    # run(`rm $catalogue_dir/merge.h5`) need to change so that it checks if merg.h5 is present before trying to remove
    run(`python $pksc2h5_script $in_halos $catalogue_dir/merge.h5`)

    # Overwrite in_halos with HDF5 halo catalogue
    in_halos = string(catalogue_dir,"/merge.h5")

end

# Read in halo positions and masses
h5in      = h5open(in_halos,"r")
halo_pos  = transpose(read(h5in["pos"]))
halo_mass = vec(read(h5in["m200m"]))
close(h5in)

# Note: h5 catalogue has groups:
# pos     position of each halo (x,y,z)
# m200m   virial mass (200 times matter density)
# m200c   virial mass (200 times critical density)
# z       redshift of each halo
# ra      right ascension
# dec     declination
# Each can be read similarly to pos and m200m above.

# Set cosmology and CIB model
cosmo = get_cosmology( h     =0.6735f0,
                       Neff  =2.99f0,   # Planck 2018 as opposed to standard model 3.046f0
                       OmegaK=0.0,      # Planck 2018 curvature density fraction
                       OmegaM=0.3138f0, # Planck 2018 matter density fraction
                       OmegaR=nothing,  # Planck 2018 radiation density fraction
                       Tcmb  =2.7255,   # Mean CMB temperature
                       w0    =-1,       # Using -1 for Lambda as opposed to Planck 2018 -1.03 +/- 0.03
                       wa    =0       ) # No coupled Dark Energy
model = CIB_Planck2013{Float32}() # NEED TO CHANGE TO Planck2018...

# Allocate some arrays and file them up for centrals and satellites
@time sources = generate_sources(model, cosmo, halo_pos, vec(halo_mass) );

# Deposit the sources into maps
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
m = HealpixMap{Float64,RingOrder}(model.nside)

for freq in [#"030", "090", "148", "219", "277", "350",
        "100", "143", "217", "353", "545", "857"] # Planck 2018 frequency bins
    @time begin
        XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, model, sources,
            fluxes_cen, fluxes_sat)
        Healpix.saveToFITS(m, "./maps/cib$(freq).fits")
    end
end
