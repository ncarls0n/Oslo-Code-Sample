Return to the [Main README](../readme.md)
# Debugging PeakPatch

IMPORTANT! Debuggers __DON'T__ work with PeakPatch! 
This is because running with a debugger requires changing optimizaiton flag from `-O3` to `-O0`, and this breaks 
PeakPatch. 
I suspect that this is because there are too many variables that are not declared properly/are not declared at all,
and when you try running the code with `-O0`, it gets confused. 
It also didn't help that PeakPatch was originally written with intel compilers in mind, which are very much not strict,
so rewriting it for gcc was a bit of a chore.
Making PeakPatch properly debuggable would require a MAJOR effort, 
so for now the best way to debug the code is with 'print' statements.
