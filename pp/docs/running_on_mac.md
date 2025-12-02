Return to the [Main README](../readme.md)
# Notes on Running Peak Patch on MacOS

Currently (Aug 2024), running on Mac without external tools doesn't work. 
If you want to run it on Mac, it is advised to do so [with nix](nix_setup.md).

Below is what has been done to make it work on Mac:

### Notes on running peak patch on macOS computers

At the moment (which is to say December 2020), I haven't gotten Peak Patch to run on a Mac. This seems mostly to be because of compatibility issues with modules meant to be loaded on Niagara vs those I have on my system. I may return to this in the future, although running on Niagara is pretty easy so maybe I won't get around to it.

What follows are some notes I have put together outlining how you'd set up a run on macOS.

(Note, all this stuff is written in ancient Python 2.7, so make sure you run it in that and not your default, which in my case would be python3.8. You can check which version of python you're running in terminal with `python --version`.)

#### -1) Make sure terminal is running an appropriate shell language

Around 2019, Apple switched the defaul shell script language used in terminal from bash to Z shell (or zsh). The `configure` script in Peak Patch is designed to be run in bash, and will fail if you try to run it in zsh, so first you'll need to revert your terminal back to bash.

So, to switch terminal back to bash

	chsh -s /bin/bash

And then restart terminal or open a new tab and it should be in bash. (Note that one can change the default shell to zsh in similar fashion: `zsh -s /bin/zsh` so you might want to do that whenever you're done with a peak patch session.)

(Note: later if I try to get this running in zsh shell, we'd presumably want to set the source to the the `.zshrc` file instead of the `.bashrc` one, but zsh doesn't automatically generate one, so first we'd have to make it.

	cd ~
	touch .zshrc
	cd <~/...>/peak-patch
	source ~/.zshrc

If you try this without further changes, peak patch will not work, so, you know, don't do it.)

Generally speaking, you'll probably want to replace `~/.bashrc` with `~/.cshrc` if you're using C shell or `~/.zshrc` if you're using Z shell, but note that other changes will have to be made if you want to run zsh as your shell language.

### 2) Configuration files and job scripts

Go to the peak patch directory `/src` and, if your machine configuration file does not already exist (e.g. Make.mach.cita for cita machines), then copy one of the existing files and edit it, e.g.

	cd <~/...>/peak-patch/src
	cp Make.mach.darwin Make.mach.your-machine-name

If you're running Peak Patch on macOS, you'll want to start with the darwin machine file (Darwin is a higher level programming language used by macOS).

Now edit those files to reflect the Fortran and C compilers and optimization flags you want to use. There are 2 files as one is for merge, which is not threaded and breaks on large runs if you compile with omp.

Next you should copy the shell script templates

	cd <~/...>/peak-patch/templates
	cp ppatchrun_darwin-single.sh ppatchrun_your-machine-name-single.sh

