# Peak Patch MP

Main documentation file for the Peak Patch Package

-----

### Getting started
* [Setting up on Niagara](docs/niagara_setup.md)
* [Setting up on CITA cluster](docs/cita_setup.md)
* [Setting up on CITA starq cluster](docs/starq_setup.md)
* [Setting up with Nix](docs/nix_setup.md)
* [Running with python](docs/running_with_python.md)
* [Running without python](docs/running_without_python.md)
* [Notes on running on Mac (doesn't currently work)](docs/running_on_mac.md)

### Running/Editing the code
* [INI File Documentation](docs/ini_file_docs.md)
* [Debugging PeakPatch](docs/debugging.md)

*Where can it run*: As of June 2024, PeakPatch can run on
SciNet Niagara supercomputer, CITA machines
(interactive ones like kingcrab or as cluster jobs like greenq),
on CITA starq, or through the
[nix](https://nixos.org/) [flake](https://nixos.wiki/wiki/Flakes).

The only thing that you have to change for different machines is the "machine"
parameter in the parameter file, for any supported machine:
niagara, cita, starq. nix does this change automatically, so you don't have to
worry about it at all.

## Directions for designing a peak patch run:

The peak-patch parameter file is mostly self explanatory, and everything you'd ever want to change is commented. But, setting up nmesh, nbuff, and ntile can be difficult without understanding the geometry of how peak-patch tiles a cubic density field. See `http://imgur.com/vL5Vuak` for image of tiling in the case `ntile=2` (or see `peak-patch/vL5Vuak.png` where I've coped it in my personal instance of the Peak Patch repository. This image file is kind of confusing, so I made a more thorough pdf explaining this `peak-patch/Parallelization.pdf`).

So, to set up box it is easiest to use `peak-patch/tools/setup-run.py`:

        python2.7 <~/...>/peak-patch/tools/setup-run.py <ngrid> <boxsize> [<cluster> <number of procesors per node> <buffersize>]

Where `<ngrid>` is an integer denoting the number of lattice points into whilch the simulation volume (including buffers) should be divided, and `<boxsize>` is a number in megaparsecs representing the width of the main simulaiton volume (excluding biffers). One can also include the optional arguments `<cluster>` as either `niagara` (the default) or `gpc`, `<number of procesors per node>` which is pretty self explanitory and set by default to 32 for `<cluster>=niagara` or 8 for `<cluster>=gpc`, and `<buffersize>` is a number in megaparsecs representing the width of the buffer region to be used by the simulation and is set as 36.0 as default.

This will give you all values that will work, and recommend you the best values that you should probably use (unless you're doing a high res run with large buffers). For example

        python2.7 tools/setup-run.py 512 512

If using best values:

        fftw_mem, 2.0 Gb        
        grid_mem, 4.01675608009 Gb
        cat_mem, 0.342726707458 Gb
        lat_mem, 0.147805440499 Gb
        s2c_mem 0.5 Gb
        
        total memory per node ~  7.00728822805 Gb
        
        ===============================
        RECOMMENDED VALUES
        ntile =  2
        nbuff =  39
        nmesh =  295
        nnodes =  1
        ===============================



## Succinct run instructions

Log into Niagara

	ssh user@niagara.scinet.utoronto.ca

Configure

	cd <~/...>/peak-patch
	./configure
        source ~/.bashrc

Load Niagara modules

	module load NiaEnv/2019b
        module load intel/2019u4
        module load intelmpi/2019u4
        module load fftw/3.3.8
        module load gsl/2.5
        module load cfitsio/3.450
        module load intelpython2/2019u5

Figure out what run parameters you should use

	python2.7 <~/...>/peak-patch/tools/setup-run.py <ngrid> <boxsize> [<cluster> <number of procesors per node> <buffersize>]

Create a run (e.g. `<run-name>` in `<~/...>/peak-patch/runs/`

	cp -r example runs/<run-name>

and then plug the values from `setup-run.py` into `<run-name>/param/param.params`). Run Peak Patch

        cd <~/...>/peak-patch/example
        python2.7 ../python/peak-patch.py param/param.params
        ./bin/hpkvd 1

Submit the job to the Niagara queue

        sbatch <short_name>_<seed>.sh

(e.g. `sbatch 512Mpc``_nb14``_13579.sh`) 



## Directions for making maps from a peak patch run

The directory `/peak-patch/src/pks2map` contains all of the mapmaking codes. Maps capabilities so far include:

1. pasted profiles:
   * tSZ   - BBPS profiles
   * kSZ   - BBPS profiles
   * Kappa - BBPS profiles
2. line profiles:
   * CO - Tony Li et al. 2016
   * HI - Villaescusa-Navarro et al. 2014

For some reason, the makefile thinks that one of the fortran scripts we need, `bbps``_profile.f90`, is supposed to be in both `peak-patch/src/pks2cmb` and `peak-patch/src/pks2map`. However, it isn't in the latter. Rather than changing the makefile and possibly introducing other problems, I just copied it, so if `peak-patch/src/pks2map/bbps``_profile.f90` isn't present, make sure to copy it from `pks2cmb`.

To use pasted profiles you first need to create the BBPS tables. To do so (using model 1) on Niagara

	cd <~/...>/peak-patch/runs/<run_name>/src/
	make make_maptable
	../bin/make_maptable maptable1 1 

Then you need to make `pks2map`

	make pks2map

Once made, we can finally produce some maps based on a Peak Patch run using `pks2map` as follows

	Usage: mpirun -np 8 pks2map <filein> <fileout> <tablefile> <profilecode> 
	any modules you may have locally in your login node will not be loaded, so you need to include these in your job scripts		[<zmin> <nside> <scramble> <center> <npix> <fov> <zmax> <chihview> <model> <PSZcut>]  
	profile choices: tsz, ksz, kap, tau, tco, thi
	
	<zmin>      - Minimum redshift for halos in map
	<nside>     - healpy map nside (total pixels = 12*nside**2)
	<scramble>  - 0 for off, 1 to randomly move halo to new position on sky
	<center>    - 0 for off, 1 to centre maps on the most massive halo
	<npix>      - number of pixels on a side for flatsky 
	<fov>       - field of view for flatsky
	<zmax>      - Maximum redshift for halos in map
	<chihview>  - 0 for off, 1 to move observer to this location
	<model>     - use model 1 (AGN feedback from BBPS)
	<PSZcut>    - 0 for off, 1 to cut halos below polynomial fit to 
	              Planck 2015 selection function  
	
	For pasted profiles `<tablefile>=maptable1` that you just created
	For Tco             `<tablefile>=/peak-patch/tables/sfr_behroozi.dat`
	for Thi             `<tablefile>=anything`

Example to make 512 CO maps from a 560Mpc box at 2.4 < z < 2.8:

	mpirun -np 8 ./bin/pks2map peakpatchrun_merge.pksc.13579 mapsout/CO 
	sfr_behroozi.dat tco 0.0 8 0 0 1024 4.8 3.0 0 1 0



## Directions for modifying the code

1) Pull the latest version

2) Create a new branch, see https://confluence.atlassian.com/display/BITBUCKET/Branching+a+Repository#BranchingaRepository-BranchaMercurialrepo)

3) Make your changes

4) Once you have finished and have a result you feel is ready to become part of the main trunk, close the branch by issuing a pull request, see https://confluence.atlassian.com/display/BITBUCKET/Branching+a+Repository#BranchingaRepository-BranchaMercurialrepo) and https://confluence.atlassian.com/display/BITBUCKET/Work+with+pull+requests



## Additional Notes

When you configure Peak Patch and set `~/.bashrc` as the source, it should add some modifications to `.bashrc` that will be critical to the function of `peak-patch.py`:

	# modifications by PykPatchAutoConf
	export PP_DIR=<~/...>/peak-patch
	export PYTHONPATH=<~/...>/peak-patch/python
	export PATH=$PATH:<~/...>/peak-patch/bin:<~/...>/peak-patch/python

Note that `PP_DIR` is the directory in which your local clone of the Peak Patch repository is saved. `peak-patch.py` interfaces with this using the built in `os` module. When `os` is first imported, it takes a snapshot of the system and saves a few key directories as strings in `os.environ`, a mapping object similar to a dictionary and consisting of a series of key strings (i.e. `'USER'`) that map to a second set of strings (i.e. `'Nate'`, in my case). From the `.bashrc` file, python figures out where the Peak Patch repository is saved. This is equivalent to the python command
```python
os.environ['PP_DIR'] = '<~/...>/peak-patch'
```
If, when running the example, you get a key error, this probably means that python was unable to locate `'PP_DIR'` when the `os` module was imported. To fix this, just run

	source ~/.bashrc

again, and it should fix this. If you then get a directory-not-found error, this means that you probably moved your Peak Patch repository. The quickest fix is going to be to go into `.bashrc` and make sure `PP_DIR`, `PYTHONPATH` and `PATH` are all pointing to the correct directories (as above).

### To clean up a run's output

When you're debugging Peak Patch, or tuning parameters, something you will find yourself doing a lot of is having to erase all the files from a test run so that you can run it again. If you don't do this, it can cause issues, so best practice is to always start with a clean run directory `<run_dir>` containing a single directory `param` with a single file, the parameter file `param.params`. The shell script `peak-patch/example/eraseoutput.sh` does this. To run it, simply

	cd <...>/peak-patch
	./example/eraseoutput.sh <...>/<run_dir>

### Moving data from SciNet to CITA

The `rsync` command allows you to quickly syncornize files (and directories if recursive flag is used) in a destination directory with a source. This is an additive synchronization in the sense that source files that are absent in the destination will be added, but destination files that are absent in the source will not be removed from the destination. No changes are made to the source. `rsync` can copy files locally or between local and a remote host (but not between two remote hosts):

	rsync [options] <source> <destination>
	rsync [options] <user>@host:<source> <destination>
	rsunc [options] <source> <user>@host:<destination>

Useful options: `-r` for recursive; `-a` for archive (essentially fast, almost complete recursion); `-v` for verbose (meaning more information is printed during the transfer); -P for long transfers that may be interupted (equivalent to `--partial --progress`). These can be combined, *e.g.* `-avP`.

Once you've ssh'd into Niagara, you need to ssh again into one of the data mover nodes, and then use `rsync` to copy it over to CITA.

        ssh <scinet_username>@niagara.scinet.utoronto.ca
        ssh nia-dm1
        rsync -avP <your_scinet_stuff> <your_CITA_username>@nuexport00.cita.utoronto.ca:/home/<your_CITA_username>/<your_cita_path>

(or alternatively you can use the second Niagara data mover `nia-dm2`, and/or the other CITA data mover node `groundhog.cita.utoronto.ca`). You can only move data 4 times in a 2-minute window using either the login nodes (for transfers < 10GB) or the data mover nodes (for transfers > 10GB), so bundling transfers is encouraged. 

The `nia-dm1` is in the home directory, so as an example, the `rsync` command might look like

        rsync -avP $SCRATCH/peak-patch/. njcarlson@nuexport00.cita.utoronto.ca:/home/njcarlson/Documents/peak-patch

Note that `rsync -avP things/. <CITA login>:<CITA dir>/stuff` will copy all the contents of the directory `things/` from SciNet and copy them into the directory `stuff/` on CITA, whereas `rsync -avP things <CITA login>:<CITA dir>/stuff` will make a subdirectory `stuff/things/` on CITA. You can use this distinction to change the name of a SciNet directory when you sync it (*e.g.* if directory `stuff` doesn't exist on CITA, you can use the former case and it will make a directory `stuff` with the same contents that `things` has on Niagara).

For more info, see the CITA Wiki's [data transfer page](https://wiki.cita.utoronto.ca/index.php/Data_Transfer_between_CITA_and_SciNet) and the SciNet documentation on [Data Management](https://docs.scinet.utoronto.ca/index.php/Data_Management#Moving_data).

### CITA
All of CITA's workstations can be accessed by SSH using a VPN. On Mac, we use the application Tunnelblick to establish a VPN connection. There is a guide to setting up Tunnelblick on the [CITA Wiki](https://wiki.cita.utoronto.ca/index.php/CITA_Remote_Access#SSH_Keys).

Once Tunnelblick is installed, you can SSH into a CITA machine by clicking on the Tunnelblick icon in the menu bar and selecting "Connect CITAvpn". Then SSH in 

        ssh user@machine.cita.utoronto.ca

Where `user` is your CITA username (*i.e.* for me it's njcarlson) and `machine` is a CITA computer. CITA has a number of shared computers, you can also ssh into your personal CITA desktop which mostly have animal names (mine for instance is sheep)

        ssh njcarlson@sheep.cita.utoronto.ca

