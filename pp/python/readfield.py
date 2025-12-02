import numpy as np
import matplotlib.pyplot as plt
from peakpatchtools import PeakPatch
import peakpatchtools as pkpt

run = PeakPatch('/scratch/m/murray/vasissua/PeakPatch/pp_runs/run1')

field_data = run.add_field()
rhog_data = run.rhog


output_file_path = '../rhog00.txt'
np.savetxt(output_file_path, rhog_data[0][0], fmt='%f')
#print(np.size(rhog_data[0][0][0]))
print(rhog_data[0][0][0])
print(rhog_data[1])
print(len(rhog_data[2][0]))
