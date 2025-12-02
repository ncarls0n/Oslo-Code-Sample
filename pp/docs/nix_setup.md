Return to the [Main README](../readme.md)
# Setting up PeakPatch with Nix

You can run PeakPatch through [nix](https://nixos.org/) [flake](https://nixos.wiki/wiki/Flakes).

Running through nix flake is recommended for personal use,
including on Linux, MacOS, or Windows via WSL.
For information about running through the flake, look into [flake.nix](./flake.nix).

The Nix flake is for use with NixOS or nix package manager. 
Nix package manager is available on Linux,
MacOS, or (if you are into this kind of stuff) Windows.
                                                         
### Purpose
With nix, you don't have to worry about dependencies, libraries, activation scripts, etc. If it works 
on one machine, it works on all of them (as long as they have nix packages)
                                                         
### Usage
To enter the PeakPatch environment, defined in this flake, after you git cloned this repository and enabled
flakes (as of August 2024, they are experimental), just run:
```
nix develop
```
All the packages will be downloaded and you will automatically enter the PeakPatch Nix environment. You
do not have to do anything else.
                                                         
### Troubleshooting
If the command fail, that is probably because you didn't enable flakes. In that case, run the command with
temporarily enabling them:
```
nix develop --experimental-features 'nix-command flakes'
```

### Limitations
Currently any healpy-related code is not supported, since healpy is not packaged in nixpkgs. I tried packaging it
myself, but didn't finish (too many errors). The code is in the flake.nix file, if you want to finish it.

Also, keep in mind that for now we are compiling the code with gcc. 
This is because intel compilers are not available in nixpkgs. 
Future optimizations might require
rewriting Makefile to work with intel compilers instead, situationally, as well as packaging the intel
compilers into nix (which I don't think is a big issue, but is an issue non the less?)
gcc might be slower than intel, plus it is more strict at compile-time (which is a good thing, since it's better
to have the code running on a stricter compiler, than the other way around).
                                                         
### Credits
The flake was packaged by Vasilii Pustovoit in April 2024

## Running PeakPatch with Nix

Now you are all set to run PeakPatch! The way that it mostly has been done so far is [with python](running_with_python.md),
but you also can try doing it [without python](running_without_python.md), i.e. directly compiling the code.

Once you are nix environment though, there are a few useful aliases that you can use, instead of constantly calling the scripts.
You can see them by calling `pphelp`. Namely, the commands I usually use, look like this:
* `ppclean` runs the cleanup script, remove everything except parameter files
* `ppcopy` copies all the peakpatch source code. This basically calls for `python peak-patch.py ./param/param.params`
* `ppcopyini` does the same as above, but with ini paramterfile instead: `python peak-patch.py ./param/parameters.ini`
* `ppcpep` copies the example param directory to your current one. 
* `hpkvdtest` runs hpkvd with sample seed 13579
* `hpkvdtestini` does the same as above, but with ini paramterfile instead

So if you want to test out if peakpatch works, just create a new directory outside of peakpatch git directory,
and run:
```
ppcpep; ppcopy; hpkvdtest 
```

After that you will see if any errors appear and debug appropriately.
