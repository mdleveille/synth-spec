"File to open .fits files for wavelength and flux data and plot \
flux vs. wavelength Using synthetic spectra from PHOENIX library"

# Load packages
from astropy.io import fits
import matplotlib.pyplot as plt 

# open fits  (use full pathname of files)
## the wave file is the same for all spectra
wave_hdu = fits.open('/Users/mleveille/synth_spec/flux_data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
## this file will need to change to plot other flux profiles. There must be a better way to do this. Maybe link to file on the web?
flux_hdu = fits.open('/Users/mleveille/synth_spec/flux_data/lte02300-0.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

# extract data from fits files
wave_data = wave_hdu[0].data         # in angstroms
flux_data = flux_hdu[0].data         # in erg/s/cm^2/cm

# convert units
flux_data = flux_data*1e-8           # in erg/s/cm^2/angstrom

# plot data
plt.plot(wave_data,flux_data)
plt.xlabel(r'Wavelength, [$\AA$]')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ylabel(r'Flux, [$ergs \cdot s^{-1} \cdot cm_{surface}^{-2}$]')
plt.title(r'Synthetic Spectra for $[Fe/H] = 0$ and $[\alpha/Fe] = 0$')
plt.show()
