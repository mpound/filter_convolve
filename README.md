# filter_convolve
Convolution of passband filters with SEDFitter models


Filter response files (dat)  with 1st line # wav = NNNN
have had the first line added by me so they can be imported into
sedfitter ( see mkfilter.py ).   For wav, I have used the mean wavelengths
for each filter from the Spanish Virtual Observatory.
VOSA http://svo2.cab.inta-csic.es/theory/fps/index.php.
The filter reponse files were also downloaded from VOSA.


WISE REFERENCE

http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html#WISEZMA

SLOAN
http://classic.sdss.org/dr7/instruments/imager/index.html

http://classic.sdss.org/dr7/algorithms/fluxcal.html

HERSCHEL

See http://herschel.esac.esa.int/twiki/pub/Public/PacsCalibrationWeb/cc_report_v1.pdf
Passbands refer to constant energy spectrum nu*Fnu = constant.
Therefore, following appendix A of Robitaille etal 
http://iopscience.iop.org/article/10.1086/512039/fulltext/
we can use the same normalization as for Spitzer IRAC

Fnu(quoted) = Int[Fnu(actual)(nu_0/nu)R(nu)dnu] / Int[(nu_0/nu)^2 R(nu)dnu]

So need to compute:
denom = Int[(nu_0/nu)^2 R(nu)dnu] 
f.response = (nu_0/nu)*Rnu/denom

Acknowledgement of SVO Filter Profile Service:

If your research benefits from the use of the SVO Filter Profile Service, we would appreciate if you could include the following acknowledgment in your publication:
This research has made use of the SVO Filter Profile Service (http://svo2.cab.inta-csic.es/theory/fps/) supported from the Spanish MINECO through grant AyA2014-55216
and we would appreciate if you could include the following references in your publication:

The SVO Filter Profile Service. Rodrigo, C., Solano, E., Bayo, A. http://ivoa.net/documents/Notes/SVOFPS/index.html
The Filter Profile Service Access Protocol. Rodrigo, C., Solano, E. http://ivoa.net/documents/Notes/SVOFPSDAL/index.html
