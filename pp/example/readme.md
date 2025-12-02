# Succinct run instructions

The bare minimum:

        ssh -Y user@niagara.scinet.utoronto.ca
        cd $SCRATCH/peak-patch
        ./setup.sh <run name>

Where the run you want to do is in the directory `<run name>` (e.g. `./example`; to work the run directory must have a subdirectory with the parameter file `<run name>/param/param.params`). That is the *most* succinct summary of run instructions on Niagara. In the subsequent sections I provide a little more detail and go over what each line in the shell script `setup.sh` does.

## Log into Niagara

        ssh -Y user@niagara.scinet.utoronto.ca

## Configure

Change directory into `$SCRATCH` (which has more room for preliminary stuff than `$HOME`) and then configure peak patch.

        cd $SCRATCH/peak-patch
        ./configure
        source ~/.bashrc

## Load Niagara modules

        module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4 fftw/3.3.8 gsl/2.5 cfitsio/3.450 intelpython2/2019u5

## Design a run

Figure out what run parameters you should use given `<boxsize>` (sidelength of real-space volume in Mpc) and `<ngrid>` (sidelength of real-space volume + buffers in cells)

        python2.7 <~/...>/peak-patch/tools/setup-run.py <ngrid> <boxsize> [<cluster> <number of procesors per node> <buffersize>]

Create a run (e.g. `<run-name>` in `<~/...>/peak-patch/runs/`

        cp -r example runs/<run-name>

and then plug the values from `setup-run.py` into `<run-name>/param/param.params`). 

## Run Peak Patch

        cd runs/<run name>
        python2.7 ../../python/peak-patch.py param/param.params
        ./bin/hpkvd 1

## Submit the job to the Niagara queue

        sbatch <short_name>_<seed>.sh

## Play with the data

Push and pull the repository over to your home computer where you can run graphics and then run

	cd <~/...>/peak-patch/runs/<run-name>
	python2.7 ../../python/halo_plotting/cataloguecheck.py output/<run_name>_merge.pksc.<seed>

## More python plotting

As it stands, the main python script for plotting peak patches is `peak-patch/python/halo_plotting/plot-great-attract.py`. To run this code, you first need to get 1LPT data out of Peak Patch. My current, inelegant solution is to run a second simmulation but tell it to only do 1LPT calculations. If you've already run `peak-patch/runs/<run_name>`, then what you need to do on Niagara is

	cd peak-patch/runs
	mkdir <run_name>_1LPT
	mkdir <run_name>_1LPT/param
	cp <run_name>/param/param.params <run_name>_1LPT/param/.

Then in `<run_name>_1LPT/param/param.params`, at line 49, change `ilpt=2` to `iltp=1` and run this second simulation volume. Once it has finished,

	cd peak-patch/runs/<run_name>
	cp -r ../<run_name>_1LPT/output output/1LPT

Then push to gitlab, pull onto your local device, and run the python script

	cd <~/...>/peak-patch/
	python3.8 python/halo_plotting/plot-great-attract.py runs/<run_name>

