Return to the [Main README](../readme.md)
# Setting up Peak Patch for CITA starq

Repeat the steps with [Niagara](niagara_setup.md), except load the following modules:
```
module load openmpi/4.1.6-gcc-ucx fftw/3.3.10-openmpi-ucx gsl/2.7.1 cfitsio/4.0.0 python/3.10.2
```

Also make sure that your submission script includes compilation of PeakPatch with the job submission.
This is because as of August 2024, compiling code on machines that are not starq would lead to errors
when you run it on starq.

Keep in mind that for now we are compiling the code with gcc. Future optimizations might require
rewriting Makefile to work with AMD compilers instead, situationally. I am honestly not sure what the difference
between gcc/AMD compilers is, you should consult CITA wiki.

### Running PeakPatch

Now you are all set to run PeakPatch! The way that it mostly has been done so far is [with python](running_with_python.md),
but you also can try doing it [without python](running_without_python.md), i.e. directly compiling the code.
