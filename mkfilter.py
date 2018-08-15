#!/usr/bin/env python
#
# Marc Pound - mpound@umd.edu 
#
# Script to convolve passband response filters with SEDFitter models
#
# Filter response files (xx.dat)  with 1st line # wav = NNNN
# have had the first line added by me so they can be imported into
# sedfitter.   For value of wav, I have used the mean wavelengths
# for each filter from the Spanish Virtual Observatory, VOSA 
# http://svo2.cab.inta-csic.es/theory/fps/index.php.
# The filter reponse files were also downloaded from VOSA.
# 
# For technical description of PACS filters see
# See http://herschel.esac.esa.int/twiki/pub/Public/PacsCalibrationWeb/cc_report_v1.pdf
# My reasoning with this code is as follows:
# Herschel passbands refer to constant energy spectrum nu*Fnu = constant.
# Therefore, following appendix A of Robitaille etal 
# http://iopscience.iop.org/article/10.1086/512039/fulltext/
# we can use the same normalization as for Spitzer IRAC.
#
# Fnu(quoted) = Int[Fnu(actual)(nu_0/nu)R(nu)dnu] / Int[(nu_0/nu)^2 R(nu)dnu]
#
# So compute:
#
#   denom = Int[(nu_0/nu)^2 R(nu)dnu] and then 
#   f.response = (nu_0/nu)*Rnu/denom
# 
# GAIA and SDSS filters functions are also "energy counting" so follow
# same prescription.


import sys
import os
import glob
from astropy import units as u
import numpy as np

# Customize directories as necessary
#sedfit_base_dir = '/n/lupus3/mpound/filter_convolve'
sedfit_base_dir = '/subaruraid/mpound/sedfittermodels'
model_dir       = sedfit_base_dir+"/models_r17/"
base_dir        = '/n/lupus3/mpound/YSOproject'
# relative directory in this repository
filter_dir      = 'filter_responses/'

sys.path.append(sedfit_base_dir+'/sedfitter')
sys.path.append(base_dir)

from sedfitter.utils import integrate
from sedfitter.filter import Filter
from sedfitter.convolve import convolve_model_dir
import filtermanage as fm
import matplotlib.pyplot as plt

#model_dir = sedfit_base_dir+"/sedfittermodels/models_r06"
models = ["s-pbhmi", "s-pbsmi", "s-p-hmi", "sp--hmi", 
          "s-p-smi", "sp--smi", "spubsmi", "s---s-i", 
          "s---smi", "s-ubsmi", "s-u-hmi" ]
#models = [ "s---s-i"]
         

# Note we use GAIA 2nd Data release filter shapes
filter_info = [
("GAIA-GAIA2r.G.dat",fm.GAIA_G2r),
("GAIA-GAIA2r.Gbp.dat",fm.GAIA_B2r),
("GAIA-GAIA2r.Grp.dat",fm.GAIA_R2r),
("GAIA-GAIA2.G.dat",fm.GAIA_G2),
("GAIA-GAIA2.Gbp.dat",fm.GAIA_B2),
("GAIA-GAIA2.Grp.dat",fm.GAIA_R2),
#("SLOAN-SDSS.u.dat",fm.SDSS_u),
#("SLOAN-SDSS.g.dat",fm.SDSS_g),
#("Generic-Bessell.V.dat",fm.BESSELL_V),
#("WISE-WISE.W1.dat",fm.WISE1),
#("SLOAN-SDSS.r.dat",fm.SDSS_r),
#("SLOAN-SDSS.i.dat",fm.SDSS_i),
#("SLOAN-SDSS.z.dat",fm.SDSS_z),
]

do_plot=False
do_convolve=True
do_normalize=False

filters = []
for filt,name in filter_info:
    # Read in unitless transmission function from file. 
    # This will get stored into f.response, but it is not yet the 
    # "response function" as defined in 
    # http://sedfitter.readthedocs.io/en/stable/convolution.html
    print("Creating %s filter from %s " % (name,filter_dir+filt))
    f = Filter.read(filter_dir+filt)
    f.name = name
    print(f.wav[0],f.nu[0].to(u.GHz))
    f.wav *= 1E-4  # Convert angstroms to microns!
    if f.wav[0] < f.wav[-1]:
      # filters must be in increasing frequency or convolution fails!
      # See https://github.com/astrofrog/sedfitter/issues/59
      f.wav = f.wav[::-1]
      f.response = f.response[::-1]

    f.nu = f.wav.to(u.Hz, equivalencies=u.spectral())
    print(f.wav[0],f.central_wavelength,f.nu[0].to(u.GHz))

    # this does not seem to make much difference because normalize() does
    # something similar
    if do_normalize:
        transmission_curve = np.copy(f.response)
    #    print("nu0,TC[0],nu[1/2],TC[1/2]",f.nu[0],f.response[0],f.nu[len(f.nu)/2],f.response[len(f.response)/2])
        nu_0 = f.central_wavelength.to(u.Hz,equivalencies=u.spectral())
        integrand = f.response*(nu_0/f.nu)**2 
        denom = integrate(f.nu,integrand) 
        #print("denom=",denom)
        response = f.response*(nu_0/f.nu)/denom
        #save_response = np.copy(response)
        f.response = response
        print(f.response)

    if do_plot:
       plt.plot(f.wav,f.response)
       plt.title(f.name)
       plt.show()
    f.normalize()
#    print("TC[0],R[0],max(TC),max(R) ",transmission_curve[0],f.response[0],np.max(transmission_curve),np.max(f.response))
#    print("normalized R[0],max(R)",f.response[0],np.max(f.response))
    # almost no difference betwen before and after normalization
    #print(save_response - f.response) 
    #print(f)
    #print(f.response)
    filters.append(f)

if do_convolve:
    for m in models:
        print("Convolving %s..." % m)
        model_path = model_dir+m
        convolve_model_dir(model_path, filters)
