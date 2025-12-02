# Peak Patch Notebooks

This directory has a series of useful Jupyter notebooks for manipulating Peak Patch and WebSky, mostly interfacing with the `peakpatchtools` module. Below I give succinct instructions for running the Jupyter notebooks in this directory remotely. This is intended as a quick reference if you're already set up to do so, if you haven't see the following section with detailed instructions.

## Succinct instructions to run notebooks in this directory remotely

To run the jupyter notebooks in this direcotry remotely, first run command

	ssh -N -f -L localhost:8884:localhost:8889 njcarlson@sheep.cita.utoronto.ca

Then in your browser enter

	localhost:8884

You may run into the issue of the `localhost` being used, in which case you'll want to change `8884` in both of the above.

## Detailed instructions to run notebooks in this directory remotely

If you haven't done this before, first in remote host, open terminal and go to this directory

	ssh username@workstation.cita.utoronto.ca
	cd <...>/peakpatch/python/notebooks
	ml python
	jupyter notebook --no-browser --port=8889

Where `<...>` will be used throughout this readme to refer to the path to your branch of Peak Patch on this workstation. Next in your local machine, open terminal and enter

	ssh -N -f -L localhost:8884:localhost:8889 username@workstation.cita.utoronto.ca

You'll be prompted to enter a password. If this all works then you can proceed to access the Jupyter notebooks by entering the following into your internet browser's search bar

	localhost:8884

You may run into the issue of the `localhost` being used, in which case you'll want to change `8884` in both of the above. If you haven't logged in before, you'll be asked to enter a password or token. You can get the token by entering the following command in the ssh tunnel to your workstation:

	jupyter notebook list

This will give you some output that looks like the following:

	Currently running servers:
	http://localhost:8889/?token=<your token> :: <...>/peakpatch/python/notebooks

Where I've replaced the actual token (which will be a random sequence of ~48 letters and numbers) with `<your token>`. You can copy the thing that is in place of `<your token>` into the web browser to set up a password, and then access your notebooks with more ease in the future.

There's a useful tutorial on doing this [here](http://amber-md.github.io/pytraj/latest/tutorials/remote_jupyter_notebook) if you run into any trouble.

