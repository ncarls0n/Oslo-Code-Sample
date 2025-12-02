Return to the [Main README](../readme.md)
# Setting up Peak Patch on SciNet Niagara Cluster

## Installation and Quick Start 

as of 5 December 2020, edited by Nathan Carlson

At the moment, I've only got Peak Patch running on the Niagara supercomputer which is mainteined by SciNet (a high performance computing group affiliated with UofT and Ontario Hospitals). I've made some notes about running Peak Patch on macOS, you can find them at the end of this document.

To get started using SciNet, you first need to set up a Compute Canada account, which is done through your supervisor. To see a detailed walkthrough of this, and of Niagara in general, see the [Niagara Docs](https://docs.scinet.utoronto.ca/index.php/Niagara_Quickstart). I'll give an abrieviated version here that is sufficient for getting Peak Patch going.

### 0) Clone the package

Once you're set up with Compute Canada and have a SciNet username (I'll use `username` in the following) ssh into one of Niagara's login nodes

	ssh -Y username@niagara.scinet.utoronto.ca

You have two personal directories on Niagara under your group's allotment

	$HOME=/home/g/groupname/username
	$SCRATCH=/scratch/h/groupname/username

`$HOME` is backed up, you should clone `peak-patch` into this directory. `$SCRATCH` isn't backed up, anything in there that isn't used for two months will get purged. It has more space than `$HOME` so that's where you should run configuration-type stuff for setting up your jobs. In other words, copy `peak-patch` into `$SCRATCH` and continue following the instructions.

To clone the package, `cd` to the directory where you want to save the repository (I'll call this path `<~/...>/`)

	cd <~/...>/
	git clone git@gitlab.com:marceloalvarez/peak-patch.git

This saves the Peak Patch repository on your local machine as a directory `<~/...>/peak-patch/`. Note that you'll need to have either "developer" or "maintainer" privileges on the GitLab project to clone the repository, so if you're just a guest, you'll need to contact George or Marcelo to get that fixed.

You'll probably also want to have a local copy of Peak Patch on your machine for convenience. I've also made a copy of the repository in my own GitLab `natecarlson/peak-patch.git` so I can make personal changes.

For actual runs, you should copy your clone of Peak Patch over from `$HOME` to `$SCRATCH`.

### 1) Configure the peak patch code

Next we need to configure the code. To do so, in the peak patch directory on Niagara   	   

	./configure
	source ~/.bashrc

Peak Patch requires a number of modules in addition to Niagara's environment module `NiaEnv/2019b` (which loads a select number of modules by default). To load an additional module `<module-name>/<module-version>` use

        module load <module-name>/<module-version>

The relevant modules to load for NiaEnv2019b are: `intel/2019u4`, `intelmpi/2019u4`, `fftw/3.3.8`, `gsl/2.5`, `cfitsio/3.450` and `intelpython2/2019u5`. To load these run

To load these run the following on Niagara

	module load intel/2019u4
	module load intelmpi/2019u4
	module load fftw/3.3.8
	module load gsl/2.5
	module load cfitsio/3.450
	module load intelpython2/2019u5

(As an alternative to the final load command, you could also use the shorter: `module load conda2` however it's probably safer to specify the exact version rather than using a pointer like `conda2` so that you know exactly what was running in the event you need to debug something later on.) 

You can also experiment with `NiaEnv/2018a`, which may be faster, in this case the modules to load are: `intel/2018.2`, `intelmpi/2018.2`, `fftw-mpi/3.3.7`, `gsl/2.4`, `cfitsio/3.430` and `anaconda2/5.1.0`.

If it isn't 2020 still and you're doing this it's gonna be a good idea to check that these are the right module versions. You can do so by running e.g. `module avail intel` and it will tell you which version of the `intel` module are available. You can probably just use the newest one.

Note that we can run `module list` to see what you've got running. It should be those six plus NiaEnv/2019b.

You can make sure they are up to date using `module spider ...`.

### 2) Setup for your device compilers

Compilers should be setup and ready to go for a run on Niagara. In the directory `<~/...>/peak-patch/src` you should see a file `Make.mach.niagara`. If you have problems with your compilers, this is where they are specified.

Next we'll look at the shell script templates in
	cd <~/...>/peak-patch/templates
here you should see `ppatchrun_niagara-single.sh` and `ppatchrun_niagara-multi.sh` (the latter is for larger runs). When you pass a job to the scheduler to be run on one of Niagara's compute nodes (via `sbatch`, see below), any modules you may have locally in your login node will not be loaded, so you need to include these in your job scripts. Peak patch compiles job scripts from these templates. Lines 11 through 17 should simply be `module load NiaEnv/2019b` followed by the six `module load` commands specified in the previous section.

### 3) Running PeakPatch

Now you are all set to run PeakPatch! The way that it mostly has been done so far is [with python](running_with_python.md),
but you also can try doing it [without python](running_without_python.md), i.e. directly compiling the code.
