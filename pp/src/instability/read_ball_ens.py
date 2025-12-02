import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt('ball_ens.out', usecols=(4,3))
plt.plot(data[:,0],data[:,1])

def quadratic(x,a,b,c): return a*x**2+b*x+c
def squareoff(x,a,c): return a*x**2+c
def simplequad(x,a): return a*x**2

a,b = curve_fit(quadratic, data[:,0],data[:,1])
plt.plot(data[:,0],quadratic(data[:,0],a[0],a[1],a[2]))
print('ΔΦ=({0})χ^2+({1})χ+({2})'.format(a[0],a[1],a[2]))
#ζ
a,b = curve_fit(squareoff, data[:,0],data[:,1])
plt.plot(data[:,0],squareoff(data[:,0],a[0],a[1]))
print('ΔΦ=({0})χ^2+({1})'.format(a[0],a[1]))

a,b = curve_fit(simplequad, data[:,0],data[:,1]-data[0,1])
plt.plot(data[:,0],simplequad(data[:,0],a[0]))
print('ΔΦ -〈ΔΦ〉= ({0})χ^2'.format(a[0]))

plt.show()

