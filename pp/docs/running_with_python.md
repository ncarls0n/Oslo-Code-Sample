Return to the [Main README](../readme.md)
# Running Peak Patch-Websky with python

as of 5 December 2020, edited by Nathan Carlson

Once you set up PeakPatch on your machine (be it Niagara, CITA, Nix, etc.), here is how you will run it with python:

### Run Peak Patch (Python)

We are now ready to run peak patch. The repository contains a sample run `<~/...>/peak-patch/example` in which we can execute a Peak Patch run using the parameter file `example/param/param.params` and the coresponding template jobscript (i.e. `<~/...>/peak-patch/templates/ppatchrun_niagara-single.sh` for small runs on Niagara). You'll want to make sure that line ~27 in `/param/param.params` reflects your current machine (i.e. `machine = 'niagara'`). Once that's good, run

	cd <~/...>/peak-patch/example
	python2.7 ../python/peak-patch.py param/param.params
	./bin/hpkvd 1

The python script will take about 30 seconds on Niagara. After it runs you should have a directory set up to run peak-patch. Next we run the executable hpkvd (for Hierarchical PeaK VoiD) to generate collapse tables. Collapse tables are saved in the file `HomelTab.dat` (for HOMogeneous ELlipsoidal collapse TABle).

You can now submit the job to the queue with

	sbatch <short_name>_<seed>.sh

(e.g. `sbatch 512Mpc_nb14_13579.sh`) 

#### 3.1) Check your output

Once Peak Patch is finished, the halo catalogues can be found in `<run-name>/output` (pre-merging/exclusion catalogues have the form `*_raw.pksc.*` and final catalogues have the form `*_merge.pksc.*`). These are binary files, which means you can't read them using a typical editor like VIM, but you can check them using the scripts in `peak-patch/tools/check_catalogue`, i.e.

	cd <...>/<run_dir>/output
	python3 <...>/peak-patch/tools/check_catalogue/check_catalogue.py ./*_merge.pksc.*

#### 3.2) Oh no! Something's gone wrong

If the job fails for some reason, you can try running it yourself. This should only be done for small runs on Niagara because the `$SCRATCH` directory has RAM limits. To run yourself without submitting to the queue, in the same directory where you ran the last command `python2.7...` (it should be `/peak-patch/example` if you're following the directions) run 

	mpirun -np 8 ./bin/hpkvd 

You will see peak-patch running, and once it is complete you will have one "raw" file in `/output/` called

	<run_name>_raw.pksc.<seed>

To merge this file

	mpirun -np 8 ./bin/merge_pkvd

This will return a final peak patch catalogue in output/ called

	output/<run_name>_merge.pksc.<seed>
	        
These are the unmerged and merged halo catalogs, respectively. The merged halo catalog is the final product and is the one you generally want to use

	OUTPUT FILES:
	This final halo catalog is a binary file with a 12 byte header,
		header = (int Nhalo, float RTHmax, float redshiftbox)
	where Nhalo is the total number of halos found, RTHmax is the radius of
	the largest halo in the simulation, and redshiftbox is the redshift of 
	the box, which will be negative for lightcone runs.
	
	This is then followed by a list of Nhalo*11 4 byte floats that represent
		(x_halo,y_halo,z_halo,vx_halo,vy_halo,vz_halo,Rth_halo,xL_halo,yL_halo,zL_halo,delta_halo)
        
	Sample python code to load a datafile:
	pkfile       = open(filein,"rb") # opening in mode 'rb' means 'read-only' + 'binary' (binary as opposed to text, meaning the file is not human readable as openned, but can be used by the computer)
	Nhalo        = np.fromfile(pkfile, dtype=np.int32, count=1)[0]
	RTHMAXin     = np.fromfile(pkfile, dtype=np.float32, count=1)[0]
	redshiftbox  = np.fromfile(pkfile, dtype=np.float32, count=1)[0]
	print "Nhalo = ", Nhalo
	
	nfloats_per_halo = 11     
	npkdata          = nfloats_per_halo * Nhalo

	peakdata = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
	peakdata = np.reshape(peakdata,(Nhalo,nfloats_per_halo))
	
	xpos   = peakdata[:,0]
	ypos   = peakdata[:,1]
	zpos   = peakdata[:,2]
	vx     = peakdata[:,3]
	vy     = peakdata[:,4]
	vz     = peakdata[:,5]
	Rth    = peakdata[:,6]
	xpos_L = peakdata[:,7]
	ypos_L = peakdata[:,8]
	zpos_L = peakdata[:,9]
	deltah = peakdata[:,10]
	
        Omega_M = 0.31
	h       = 0.68
	rho     = 2.775e11 * Omega_M * h**2
	M       = 4.0/3*np.pi * Rth**3 * rho

### 4) Make a plot of the output to see how it looks

	python2.7 ../python/halo_plotting/cataloguecheck.py output/<run_name>_merge.pksc.<seed>

Niagara isn't set up to do graphical stuff, so you'll want to copy the output directory (i.e. `<~/...>/peak-patch/example/output`) back to your local machine. I would recommend you set up another git repository to commit your runs to and then clone that repository to your local machine as directly copying from Niagara to a local machine is made difficult by SciNet's firewalls.


