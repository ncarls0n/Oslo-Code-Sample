Peak Patch julia scripts are set up to run on CITA machines.

To download julia, go to the julia language downloads page https://julialang.org/downloads/ and select the Generic Linux on x86 (glibc), copy the link and run wget on CITA to download. As of Januaru 2023, this command looks like

    wget https://julialang-s3.julialang.org/bin/mac/x64/1.8/julia-1.8.5-mac64.dmg
    tar -xzf julia-1.8.5-linux-x86_64.tar.gz

This will install julia in the current directory. The julia executable is in ./julia-1.8.5/bin/. If you add this to your PATH (i.e. add `export PATH=$PATH:/home/njcarlson/julia-1.8.5/bin` to `~/.bashrc`), then you can run julia with command

    julia

Next we'll need to add some things in the package manager. To access the package manager, when running interactive julia, type `]`. And add the following packages

    add IJulia
    dev https://github.com/xzackli/XGPaint.jl.git
    add Healpix

IJulia will allow you to make Jupyter notebooks. The second downloads Zack Li's XGPaint repository which makes CIB maps. And the final adds HEALPix support so that we can make maps in the HEALPix projection which can then be plotted using the python scripts in `peak-patch/python/plotting/skymap_plotting/.`.
