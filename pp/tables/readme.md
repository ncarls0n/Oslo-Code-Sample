# Tables

This directory contains tabulated values of power spectra and transfer functions used for generating cosmological fields from Gaussian white noise. 

For a roughly homogeneous and isotropic field, the two-point correlation function $\xi(r)$ describes how closely two points in the field are expected to be correlated if they are separated by a distance $r$. Its Fourier transform is the power spectral density (often called power spectrum by physicists) $\mathcal{F}[\xi(r)](\mathbf{k}) = P(k)$. Convolving a Gaussian white-noise field $N(\mathbf{x})$ with the correlation function gives a realization of a field with that correlation

$$ \tilde{\delta}(\mathbf{k}) = P(k) \tilde{N}(\mathbf{k}) $$

$\tilde{N}(\mathbf{k})$ is the Fourier-space Gaussian white noise field.

These spectra can be made using the python script `peak-patch/tools/powerspectrum_create.py`. See it for more details.

This directory contains tables (in order typically displayed by VIM):

1. Chatoric billiards zeta(chi) tables:
   1. `FNL_spike_w3_piv12.dat`
   1. `FNL_spike_w3_piv12_2e-4.dat`
   1. `FNL_spike_w3_piv14.dat`
   Each table is a discretisation of the power spectrum of non-Gaussian structures, formatted into two columns:
   1. column 1: wavenumber $k$
   1. column 2: power spectrum $\langle | \chi |^2 \rangle$ of non-Gaussian preheating features

1. An approximation of the colour map used by the Planck collaboration
   `Planck_Parchment_RGB.txt`
   (Adapted from script by Andrea Zonca, see their blog post for more information: http://zonca.github.io/2013/09/Planck-CMB-map-at-high-resolution.html)

1. BBPS power spectra:
   1. `bbps_1.tab`
   1. `bbps_model_1.tab`
   1. `bbps_power_spec.dat`
   1. `bbps_power_spec_NG.dat`
   These are different versions of the NFW-like gas pressure profiles that we ascribe to interiors of Peak Patch dark matter halos when we're making high resolution sky maps in WebSky from the relatively low-resolution Peak Patch dark matter halo catalogues. These gas pressure profiles are referred to as BBPS (or Battaglia et al.) after the authors of the series of papers in which they were introduced (arXiv:1109.3709, 1109.3711, 1209.4082, 1405.3346). The exact forms of the profiles used by WebSky are outlined in the WebSky1.0 paper (arXiv:2001.08787). The profile used for a given Peak Patch run is specified in the parameter file by variable `tabfile_map`. Typically `bbps_1.tab` is used.
   Formatted as a binary data file of 32-bit values with a header:
   1. `nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt` where `nchit` is the number of comoving halo distance bins in the table (an integer), `nmht` is the number of halo mass bins in the table (an integer), `nrt` is the number of halo raius bins (an integer), and the subsequent entries are maximum and minimum values Mpc/h and Msol as floating point numbers.
   These tables can be made and read using subroutines `makemaptable` and `loadmaptable` in `src/pks2map/maptable.f90`.

1. Harmonic tSZ Compton parameter map power spectrum $C_\ell^{yy}$ in dimensionles $(\Delta T/T)^2$ units from Table 4 of Boillet et al 2018 arXiv:1712.00788.
   `boillet2018_tsz_cl.txt`
   The columns are
   1. $\ell_\text{eff}$
   1. $10^{12} D_\ell$
   1. $\sigma_\ell$
   1. $10^{12} D_\ell$ best fit

1. Don't remember what this does...
   1. `deltaN-LUT-1.875`
   1. `deltaN-LUT-1.875_adjust.py`
   1. `deltaN-LUT-1.875_s4000n2000.dat`



1. Planck Collaboration 2018 $P_{\zeta\zeta}$ and transfer function tables:
   `planck18_intermittent.dat`
   `planck2018_powerspectrum.dat`


1. Star Formation Rate table used in field component for pks2map
   `sfr_behroozi.dat`


1. Python power spectrum plotting scripts
   `compare_spectra.py`
   `plot_spectra.py`

1. Harmonic tSZ Compton parameter map power spectrum C_l^yy in dimensionless $(\Delta T/T)^2$ units from Table 2 of Planck 2015 Results XXII.
   `planck2015_tsz_cl.txt`
   The columns are:
   1. $\ell_\text{min}$
   1. $\ell_\text{max}$
   1. $\ell_\text{eff}$
   1. $\ell(\ell+1)C_\ell/2\pi$ power spectrum in variance per $ln(\ell)$ form in units of $[10^{12} y^2]$
   1. $\sigma_\text{stat}$ statistical uncertainty in units of $[10^{12} y^2]$
   1. $\sigma_\text{fg}$ foreground uncertainty in units of $[10^{12} y^2]$
   1. Best fit for $\ell(\ell+1)C_\ell/2\pi$ in units of $[10^{12} y^2]$

