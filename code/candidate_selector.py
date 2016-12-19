"Read flux data .fits files and save 0.5 M_sun M-dwarf and 1.3 BSS"

# load packages and values
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy import units as u 
from astropy.constants import M_sun
from astropy.constants import L_sun
from pysynphot import observation                                                # must cite pysynphot; see notebook page 10
from pysynphot import spectrum                                                   # must cite pysynphot; see notebook page 10

# load wavelength grid	
wave_hdu  = fits.open('../flux_data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
wave_data = wave_hdu[0].data                                                     # in angstroms

# manipulate parameters
d = 300                                                                            # distance at which we calculate flux (in parsecs)
xmin = 15830                                                                     # (for plotting) in angstroms
xmax = 15950                                                                     # (for plotting) in angstroms

# set constants and create lists
M_sun   = M_sun.value
L_sun_W = L_sun.value
L_sun   = L_sun_W*1e7

# Create empty arrays to be filled later
Mdwarfs      = []
BSSs         = []
suns         = []
Flux_Mdwarfs = []
Flux_BSSs    = []

# read in a list of all 841 fits files; candidate stars will be selected from this list 
## genfromtxt works better than loadtxt; it keeps the filename as a string
files = np.genfromtxt('../flux_data/flux_data.txt',dtype=str)

# this loop opens all the fits files and evaluates their masses then makes lists of candidate stars
##the first 'if' statement excludes a file that is missing the keyword 'PHXMASS'
for i in np.arange(0,len(files)-1):
	if files[i] != 'lte05300-0.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits':
		hdu = fits.open("../flux_data/"+str(files[i]))
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

# save list of Mdwarfs to text file (Mdwarf_list.txt in star_arrays directory)
file = open('../star_arrays/Mdwarf_list.txt', 'w')
for i in np.arange(0,len(Mdwarfs)-1):
	file.write("%s \n" % Mdwarfs[i])
file.close()

# save list of BSSs to text file (BSSs_list.txt in star_arrays directory)
file = open('../star_arrays/BSSs_list.txt', 'w')
for i in np.arange(0,len(BSSs)-1):
	file.write("%s \n" % BSSs[i])
file.close()

# save list of suns to text file (suns_list.txt in star_arrays directory)
file = open('../star_arrays/suns_list.txt', 'w')
for i in np.arange(0,len(suns)-1):
	file.write("%s \n" % suns[i])
file.close()

# # FOR HR DIAGRAM MATCHING (to determine what's on the main sequence)
# ## loop to print values for stellar effective temperature and luminosity of candidate stars
# ### to determine which stars are on the main sequence
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
MS_BSS       = BSSs[5]
MS_Mdwarf    = Mdwarfs[3]
solar_analog = suns[50]

# make flux arrays for best candidate stars
MS_Mdwarf_hdu       = fits.open('../flux_data/'+str(MS_Mdwarf))
MS_Mdwarf_flux_data = MS_Mdwarf_hdu[0].data*1e-8                                 # in erg/s/cm^2/angstrom

MS_BSS_hdu       = fits.open('../flux_data/'+str(MS_BSS))
MS_BSS_flux_data = MS_BSS_hdu[0].data*1e-8                                       # in erg/s/cm^2/angstrom

solar_analog_hdu       = fits.open('../flux_data/'+str(solar_analog))
solar_analog_flux_data = solar_analog_hdu[0].data*1e-8                           # in erg/s/cm^2/angstrom

# # general for any file (need to define file = )
# plot_hdu = fits.open('../flux_data/'+str(file))
# flux_data = plot_hdu[0].data*1e-8                                              # in erg/s/cm^2/angstrom

# rebin data to lower resolution from R = 100,000 to R = 40,000
## create lo-res wavelength grid with binsize=0.0375nm
lo_res_wav = np.arange(xmin,xmax,0.375)                                        # see notebook pg. 12 for calculation of this binsize 

def rebin_spec(wave, specin, wavnew):                                            # wavnew is the new wave grid
	spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)                  # configures wave/spectrum for pysynphot
	f = np.ones(len(wave))                                                       # preserve all flux at all wavelengths (no throughput curve)
	filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')          # apply throughput curve 'filter' (not necessary for us?)
	obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
	return obs.binflux

MS_Mdwarf_binflux_data    = rebin_spec(wave_data,MS_Mdwarf_flux_data,lo_res_wav)         # this is the re-binned flux array 
MS_BSS_binflux_data       = rebin_spec(wave_data,MS_BSS_flux_data,lo_res_wav)            # this is the re-binned flux array 
solar_analog_binflux_data = rebin_spec(wave_data,solar_analog_flux_data,lo_res_wav)

# calculate flux at distance d, from stellar surface
Mdwarf_Reff       = MS_Mdwarf_hdu[0].header['PHXREFF']                                 # effective stellar radius of the best-fit Mdwarf in cm
BSS_Reff          = MS_BSS_hdu[0].header['PHXREFF']                                    # effective stellar radius of the best-fit BSS in cm
solar_analog_Reff = solar_analog_hdu[0].header['PHXREFF']

d_cm = u.pc.to(u.cm,d)

MS_Mdwarf_final_flux    = MS_Mdwarf_binflux_data*(Mdwarf_Reff/d_cm)**2
MS_BSS_final_flux       = MS_BSS_binflux_data*(BSS_Reff/d_cm)**2
solar_analog_final_flux = solar_analog_binflux_data*(solar_analog_Reff/d_cm)**2

# make a plot
## define x-axis range
#xmin = 14995                                                                     # in angstroms
#xmax = 15005                                                                     # in angstroms

plt.clf()

plt.figure(1)
plt.plot(lo_res_wav,MS_BSS_final_flux,'r-')
plt.xlabel(r'Wavelength, [$\AA$]')
plt.xlim(xmin,xmax)
plt.ylabel(r'Flux, [$erg / s / cm^{2} / \AA$]')
plt.title(r'$T_{eff} = %dK, M = %.2f M_{\odot}, log(g) = %.2f$' % (MS_BSS_hdu[0].header['PHXTEFF'],MS_BSS_hdu[0].header['PHXMASS']/1000/M_sun,MS_BSS_hdu[0].header['PHXLOGG']))

# plt.figure(2)
# plt.plot(wave_data,MS_BSS_flux_data,'b-')
# plt.xlabel(r'Wavelength, [$\AA$]')
# plt.xlim(xmin,xmax)
# plt.ylabel(r'Flux, [$erg / s / cm^{2} / \AA$]')
# plt.title(r'$T_{eff} = %dK, M = %.2f M_{\odot}, log(g) = %.2f$' % (MS_BSS_hdu[0].header['PHXTEFF'],MS_BSS_hdu[0].header['PHXMASS']/1000/M_sun,MS_BSS_hdu[0].header['PHXLOGG']))

plt.show()







