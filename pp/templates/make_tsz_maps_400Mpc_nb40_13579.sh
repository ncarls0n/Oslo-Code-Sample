#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=6:00:00
#SBATCH --job-name mapmaking
#SBATCH --output=mapmaking_%j.txt
#SBATCH --mail-type=END

# Load Niagara modules
module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4 fftw/3.3.8 gsl/2.5 \
    cfitsio/3.450 intelpython2/2019u

# Go to SLURM submit directory
cd $SLURM_SUBMIT_DIR

dir_bin=../bin
ppfile=../output/400Mpc_n400_nb40_nt1_merge.pksc.13579
tabfile=../tables/bbps_1.tab
profile=tsz
zmin=0.0
nside=512
scramble=0
center=0
npix=512
fov=10.0
zmax=1.245
chihview=0
model=1
PSZcut=0
base=400Mpc_n400_nb40_nt1_tsz_13579

stdout=$base.stdout
stderr=$base.stderr

if [ -f $stdout ] ; then rm -f $stdout; fi
if [ -f $stderr ] ; then rm -f $stderr; fi
if [ ! -d $base ] ; then mkdir $base; fi

# Set OpenMP threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun -np 1 ${dir_bin}/pks2cmb $ppfile $base/hp $tabfile $profile $zmin \
    $nside $scramble $center $npix $fov $zmax $chihview $model $PSZcut \
    2>>$stderr 1>>$stdout


# # Use GNU Parallel to run a single map MPI run for each Websky profile
# # specified in the parameter file for this run
# cat maprun_list | parallel --jobs 4 --sshloginfile \
#     $SLURM_JOB_NODELIST --workdir $SLURM_SUBMIT_DIR
