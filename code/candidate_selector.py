"Read flux data .fits files and save 0.5 M_sun M-dwarf and 1.3 BSS"

# load packages and values
from astropy.io import fits
from astropy import units as u 
import numpy as np
from astropy.constants import M_sun
from astropy.constants import L_sun
import matplotlib.pyplot as plt 
from pysynphot import observation
from pysynphot import spectrum

# load wavelength grid	
wave_hdu = fits.open('/Users/michael/synth_spec/flux_data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
wave_data = wave_hdu[0].data         # in angstroms

# manipulate parameters
d = 1                 # distance at which we calculate fluxin parsecs
xmin = 5180           # in angstroms
xmax = 5190           # in angstroms

# set constants and create lists
M_sun = M_sun.value
L_sun_W = L_sun.value
L_sun = L_sun_W*1e7

# Create empty arrays to be filled later
Mdwarfs = []
Flux_Mdwarfs = []
BSSs = []
Flux_BSSs = []
suns = []

# genfromtxt works better than loadtxt; it keeps the filename as a string
files = np.genfromtxt('/Users/mleveille/synth_spec/flux_data/flux_data.txt',dtype=str)

# this loop reads in all the fits files and evaluates their masses 
##(the first 'if' statement excludes a file that is missing the keyword 'PHXMASS')
for i in np.arange(0,len(files)-1):
	if files[i] != 'lte05300-0.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits':
		hdu = fits.open("/Users/mleveille/synth_spec/flux_data/"+str(files[i]))
		M = hdu[0].header['PHXMASS']/1000 # in kg

		# selects ~0.5 M_sun M-dwarf candidate files
		if 0.49 < M/M_sun < 0.51:
			#print('%s with M = %f' % (files[i][3:13], M/M_sun))
			Mdwarfs.append(files[i])

		# selects ~1.3 M_sun BSS candidate files
		if 1.29 < M/M_sun < 1.31:
			#print('%s with M = %f' % (files[i][3:13], M/M_sun))
			BSSs.append(files[i])

		if 0.9 < M/M_sun < 1.1:
			#print('%s with M = %f' % (files[i][3:13], M/M_sun))
			suns.append(files[i])

		hdu.close()

# this loop calculates the flux detected at a distance, d, for the M-dwarfs
for i in np.arange(0,len(Mdwarfs)):
	# loads relevant numbers from fits headers
	hdu = fits.open("/Users/mleveille/synth_spec/flux_data/"+str(Mdwarfs[i]))
	F_stellar_surface = hdu[0].data*1e-8                                         # in erg/s/cm^2/angstrom
	R_eff = hdu[0].header['PHXREFF']                                             # in cm
	
	# calculate the flux a distance, d, from the star (FOR M-DWARFS)
	d_cm = u.pc.to(u.cm,d)
	F_earth = F_stellar_surface*(R_eff/d_cm)**2
	Flux_Mdwarfs.append(F_earth)
	hdu.close()

# this loop calculates the flux detected at a distance, d, for the BSSs
for i in np.arange(0,len(BSSs)):
	# loads relevant numbers from fits headers
	hdu = fits.open("/Users/mleveille/synth_spec/flux_data/"+str(BSSs[i]))
	F_stellar_surface = hdu[0].data*1e-8                                         # in erg/s/cm^2/angstrom
	R_eff = hdu[0].header['PHXREFF']                                             # in cm
	
	# calculate the flux at a distance, d, from the star (FOR BSSs)
	d_cm = u.pc.to(u.cm,d)
	F_earth = F_stellar_surface*(R_eff/d_cm)**2
	Flux_BSSs.append(F_earth)
	hdu.close()

# # loop to print values for stellar effective temperature and luminosity of candidate stars
# ## to determine which stars are on the main sequence
# for i in np.arange(0,len(Mdwarfs)):
# 	hdu = fits.open("/Users/mleveille/synth_spec/flux_data/"+str(Mdwarfs[i]))
# 	T = hdu[0].header['PHXTEFF']
# 	L = hdu[0].header['PHXLUM']
# 	print('For Mdwarf %s, T_eff= %dK Luminosity= %.2f L_sun' % (i,T,L/L_sun))
#	hdu.close()

# for i in np.arange(0,len(BSSs)):
# 	hdu = fits.open("/Users/mleveille/synth_spec/flux_data/"+str(BSSs[i]))
# 	T = hdu[0].header['PHXTEFF']
# 	L = hdu[0].header['PHXLUM']
# 	print('For BSS %s, T_eff= %dK Luminosity= %.2f L_sun' % (i,T,L/L_sun))
#	hdu.close()

# ********** BEST MAIN SEQUENCE FITS ARE:
MS_BSS    = BSSs[5]
MS_Mdwarf = Mdwarfs[3]
solar_analog = suns[50]

# input the desired file for plotting in flux_file
flux_file = MS_BSS                                                               # CHANGE THIS TO PLOT DIFFERENT STARS
plot_hdu = fits.open('/Users/mleveille/synth_spec/flux_data/'+str(flux_file))
flux_data = plot_hdu[0].data

# resample data to lower resolution from 100,000 to 40,000
## create lo-res wavelength grid with 
lo_res_wav = np.arange(500,50000,0.0375)     

def rebin_spec(wave, specin, wavnew):                                            # wavnew is the new wave grid
	spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)                  # configures wave/spectrum for pysynphot
	f = np.ones(len(wave))                                                       # preserve all flux at all wavelengths (no throughput curve)
	filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')          # apply throughput curve 'filter' (not necessary for us?)
	obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
	return obs.binflux

new_binflux = rebin_spec(wave_data,flux_data,lo_res_wav)

# make a plot
plt.plot(wave_data,plot_hdu[0].data,'bo')
plt.plot(lo_res_wav,new_binflux,'ro')
plt.xlabel(r'Wavelength, [$\AA$]')
plt.xlim(xmin,xmax)
plt.ylabel(r'Flux, [$erg / s / cm^{2} / \AA$]')
plt.title(r'$T_{eff} = %dK, M = %.2f M_{\odot}, log(g) = %.2f$' % (plot_hdu[0].header['PHXTEFF'],plot_hdu[0].header['PHXMASS']/1000/M_sun,plot_hdu[0].header['PHXLOGG']))
plt.show()



































