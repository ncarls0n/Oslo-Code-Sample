Return to the [Main README](../readme.md)

# Running PeakPatch without python

With last edits (August 2024), it became possible to partially run peakpatch without python.
You might need this if you want to use peakpatch as a module for other programs (i.e MUSIC), 
want more control over which parts of the program are run, editing the source and don't want
to recompile EVERYTHING every time you make a change, etc.

Running without python is best with `parameters.ini` config file, rather than `param.params`.

Here are the steps:
1. (Optional) Clone the peakpatch source directory into `<your_rundir>/src`
2. Create all necessary directories in the directory of the run. These directories include: bin, logfiles, output.
3. Make sure to copy tables with all their contents
4. Copy the parameter file (parameters.ini) from the example directory `<pp_source>/example/param>` and edit it as needed. It is advised to put it into `param` directory.
5. In src directory, compile PeakPatch with:
```
make hpkvd CONFIG_FILE=<path/to/parameters.ini>
```
6. Compile other parts of peakpatch (`merge_pkvd`, `pks2map`, etc.) by replacing `hpkvd` in command above
7. Generate collapse tables from your run directory by `./bin/hpkvd 1 <seed> <path/to/parameters.ini>`. 
8. (Optional) Check that tables make sense with `od -f HomelTab.dat` (see that the file isn't too small or doesn't contain values that don't make sense, i.e. very large or very small)
9. Now, finally run peakpatch from your run directory by `./bin/hpkvd 0 <seed> <path/to/parameters.ini>`

## The need for CONFIG_FILE
You might have noticed that makefile requires `CONFIG_FILE` as an argument to run. 
This is because some parameters are supposed to be hardcoded in, namely `n1, n2, n3`.
If you are running the code with `param.params` file, using `peak-patch.py`, there is no need in specifying
`CONFIG_FILE`, since the code was edited when you copied it over. 
If you are running with `parameters.ini`, it is required to specify this file, since otherwise you will get
a bunch of errors during compilation.
