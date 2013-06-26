# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import astropy.io.fits as pyfits
import astropy
import astropy.io.ascii as asciitable
from copy import deepcopy
import glob
import os
from scipy import optimize
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.io import fits


import input_info
reload(input_info)
import photometry_wcs
reload(photometry_wcs)
import photometry_both
reload(photometry_both)

datapath = input_info.resultfolder + '*YSO*/*_wcs.fits' 
imlist = glob.glob(datapath)
imlist.sort()

# this finds the observations which are not actually centered on the cluster.
badlist = photometry_both.make_badlist(imlist)
# this should be empty because only the good nights should have been copied to the resultfolder.
print badlist

# create the filenames for the individual PSF star files
outlist = deepcopy(imlist)
for i in np.arange(0,len(outlist)):
    outlist[i] = outlist[i] + '.pstbyhand'

# get a list of the results of the aperture photometry (.mag.1 files):
magpath = input_info.resultfolder + '*YSO*/*_wcs.fits.magapphot' 
maglist=glob.glob(magpath)
maglist.sort()

# now do the actual conversion of the psf star list to image coordinates for each observation.
for i in np.arange(0, len(imlist)):
    photometry_wcs.make_pst_file_from_ds9reg(input_info.psfstarfile, maglist[i], imlist[i], outlist[i])
    plt.close()




