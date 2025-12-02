Return to the [Main README](../readme.md)
# Setting up Peak Patch for CITA clusters

Repeat the steps with [Niagara](niagara_setup.md), except load the following modules:
```
module load openmpi/4.1.6-gcc-ucx fftw/3.3.10-openmpi-ucx gsl/2.7.1 cfitsio/4.0.0 python/3.10.2
```

The code should work well on both interactive (kingcrab, calamari, ...) and Sunnyvale
(workq, greenq, ...) clusters. Tested on kingcrab.

Keep in mind that for now we are compiling the code with gcc. Future optimizations might require
rewriting Makefile to work with intel compilers instead, situationally.

### Running PeakPatch

Now you are all set to run PeakPatch! The way that it mostly has been done so far is [with python](running_with_python.md),
but you also can try doing it [without python](running_without_python.md), i.e. directly compiling the code.
