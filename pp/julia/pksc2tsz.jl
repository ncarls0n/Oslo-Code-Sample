# pksc2cib.jl
#
# This julia script generates CIB maps
# 
# USAGE:
#
#     cd <run_dir>
#     julia pksc2cib.jl <catalogue>.pksc.<seed>

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

halo_ang  = Healpix.vec2ang(halo_pos)



#=

# Note: h5 catalogue has groups:
# pos     position of each halo (x,y,z)
# m200m   virial mass (200 times matter density)
# m200c   virial mass (200 times critical density)
# z       redshift of each halo
# ra      right ascension
# dec     declination
# Each can be read similarly to pos and m200m above.

# Set cosmology and pressure profile
Omega_b = 0.0493f0
Omega_m = 0.3138f0
f_b     = Omega_b/Omega_m
cosmo = get_cosmology( h     =0.6735f0,
                       Neff  =2.99f0,   # Planck 2018 as opposed to standard model 3.046f0
                       OmegaK=0.0,      # Planck 2018 curvature density fraction
                       OmegaM=Omega_m,  # Planck 2018 matter density fraction
                       OmegaR=nothing,  # Planck 2018 radiation density fraction
                       Tcmb  =2.7255,   # Mean CMB temperature
                       w0    =-1,       # Using -1 for Lambda as opposed to Planck 2018 -1.03 +/- 0.03
                       wa    =0       ) # No coupled Dark Energy
profile = XGPaint.BattagliaProfile( f_b , cosmo )

# beam stuff (not used in this particular script)_
N_logtheta = 512
rft = RadialFourierTransform(n=N_logtheta, pad=256)

# Making or reading table for table read
model_file::String = "cached_battaglia.jld2"
if isfile(model_file)
    print("Found cached Battaglia profile model. Loading from disk.\n")
    model = load(model_file)
    prof_logθs, prof_redshift, prof_logMs, prof_y = model["prof_logθs"],
        model["prof_redshift"], model["prof_logMs"], model["prof_y"]
else    
    print("Didn't find a cached profile model. Computing and saving.\n")
    logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
    @time prof_logθs, prof_redshift, prof_logMs, prof_y = profile_grid(p;
        N_logθ=N_logθ, logθ_min=logθ_min, logθ_max=logθ_max)
    save(model_file, Dict("prof_logθs"=>prof_logθs, 
        "prof_redshift"=>prof_redshift, "prof_logMs"=>prof_logMs, "prof_y"=>prof_y))
end

# interpolate
itp = Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, prof_logθs, prof_redshift, prof_logMs);

##
function paint_map!(m, p::XGPaint.AbstractProfile, psa, sitp, masses,
                    redshifts, αs, δs, irange; mult=4)
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax = XGPaint.θmax(p, mh * XGPaint.M_sun, z, mult=mult)
        profile_paint!(m, α₀, δ₀, psa, sitp, z, mh, θmax)
    end
end

function chunked_paint!(m, p::XGPaint.AbstractProfile, psa, sitp, masses,
                        redshifts, αs, δs; mult=4)
    m .= 0.0
 85
 86     N_sources = length(masses)
 87     chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
 88     chunks = chunk(N_sources, chunksize);
 89
 90     Threads.@threads for i in 1:Threads.nthreads()
 91         chunk_i = 2i
 92         i1, i2 = chunks[chunk_i]
 93         paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2, mult=mult)
 94     end
 95
 96     Threads.@threads for i in 1:Threads.nthreads()
 97         chunk_i = 2i - 1
 98         i1, i2 = chunks[chunk_i]
 99         paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2, mult=mult)
100     end
101 end
102
103 ##
104 # cut = eachindex(halo_mass)[begin:end-5]
105 m = Enmap(zeros(shape), wcs)
106
107 print("Painting map.\n")
108 @time chunked_paint!(m, p, psa, sitp, halo_mass, redshift, ra, dec, mult=cutoff)
109
110
111 #
112 write_map(
113     "/mnt/raid-cita/mlokken/buzzard/ymaps/ymap_buzzard_standard_bbps_car_zgtrpt03.fits",
114     m)
115
116 #
117 using PyPlot
118 plt.clf()
119 plt.figure()
120 plt.imshow(log10.(m.data'))
121 plt.axis("off")
122 plt.savefig("test_buzzard.png", bbox_inches="tight",pad_inches = 0)
123 plt.gcf()





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
        Healpix.saveToFITS(m, "./cib$(freq).fits")
    end
end

=#


