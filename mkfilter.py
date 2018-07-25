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
from sedfitter.utils import integrate
from sedfitter.filter import Filter
from sedfitter.convolve import convolve_model_dir
import numpy as np
import filtermanage as fm

# Customize directories as necessary
sedfit_base_dir = '/n/subaruraid/mpound'
base_dir        = '/n/lupus3/mpound/YSOproject'
filter_dir      = base_dir+'/filter_convolve'
model_dir       = sedfit_base_dir+"/sedfittermodels/models_r17/"

sys.path.append(sedfit_base_dir+'/sedfitter')
sys.path.append(base_dir)

#model_dir = sedfit_base_dir+"/sedfittermodels/models_r06"
models = ["s-pbhmi", "s-pbsmi", "s-p-hmi", "sp--hmi", 
          "s-p-smi", "sp--smi", "spubsmi", "s---s-i", 
          "s---smi", "s-ubsmi", "s-u-hmi" ]

# Note we use GAIA 2nd Data release filter shapes
filter_info = [
("GAIA-GAIA2r.G.dat",fm.GAIA_G),
("GAIA-GAIA2r.Gbp.dat",fm.GAIA_B),
("GAIA-GAIA2r.Grp.dat",fm.GAIA_R),
("SLOAN-SDSS.u.dat",fm.SDSS_u),
("SLOAN-SDSS.g.dat",fm.SDSS_g),
("SLOAN-SDSS.r.dat",fm.SDSS_r),
("SLOAN-SDSS.i.dat",fm.SDSS_i),
("SLOAN-SDSS.z.dat",fm.SDSS_z),
]

filters = []
for filt,name in filter_info:
    # Read in unitless transmission function from file. 
    # This will get stored into f.response, but it is not yet the 
    # "response function" as defined in 
    # http://sedfitter.readthedocs.io/en/stable/convolution.html
    f = Filter.read(filter_dir+filt)
    f.name = name
    transmission_curve = np.copy(f.response)
    print("Creating %s filter",%name)
#    print("nu0,TC[0],nu[1/2],TC[1/2]",f.nu[0],f.response[0],f.nu[len(f.nu)/2],f.response[len(f.response)/2])
    nu_0 = f.central_wavelength.to(u.Hz,equivalencies=u.spectral())
    integrand = f.response*(nu_0/f.nu)**2 
    denom = -integrate(f.nu,integrand) # take negative because frequncy f.nu is decreasing 
    #print("denom=",denom)
    response = f.response*(nu_0/f.nu)/denom
    #save_response = np.copy(response)
    f.response = response
#    print("TC[0],R[0],max(TC),max(R) ",transmission_curve[0],f.response[0],np.max(transmission_curve),np.max(f.response))
    f.normalize()
#    print("normalized R[0],max(R)",f.response[0],np.max(f.response))
    # almost no difference betwen before and after normalization
    #print(save_response - f.response) 
    filters.append(f)

foreach m in models:
    print("Convolving %s..." % m)
    model_path = model_dir+m
    convolve_model_dir(model_path, filters)
