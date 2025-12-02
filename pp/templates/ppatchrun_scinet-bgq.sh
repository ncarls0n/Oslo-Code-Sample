#!/bin/bash
# @ job_name           = SNAME_REPLACE
# @ job_type           = bluegene
# @ comment            = "BGQ Job SUBBLOCK "
# @ error              = $(job_name).$(Host).$(jobid).err
# @ output             = $(job_name).$(Host).$(jobid).out
# @ bg_size            = NNODES_REPLACE
# @ wall_clock_limit   = TLIMIT_REPLACE
# @ bg_connectivity    = Torus
# @ queue

seed=SEED_REPLACE
lname=LNAME_REPLACE

ntasks=NTASKS_REPLACE
nthreads=OMPNT_REPLACE
rpn=TPNODE_REPLACE

while : 
do
  stdout=logfiles/$lname\_$seed.stdout
  stderr=logfiles/$lname\_$seed.stderr

  runjob --np $ntasks --envs OMP_NUM_THREADS=$nthreads --ranks-per-node=$rpn --cwd=$PWD : ./bin/hpkvd $seed 2> $stderr 1> $stdout

  ((seed+=2))

done



