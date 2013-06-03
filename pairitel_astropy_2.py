# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import pyfits
import astropy
import astropy.io.ascii as asciitable
from copy import deepcopy
import glob
import os
from scipy import optimize
from astropy.table import Table, Column

import input_info
reload(input_info)
import photometry_wcs
reload(photometry_wcs)
import photometry_both
reload(photometry_both)

# take the region file containing the psf photometry stars of the masterimage, and make .coo-like files from that for each image. This means getting the image coordinates for the star centers in each image, and giving them fantasy magnitudes to start with.
fantasymag = 20.

# get all image names for psf photometry (i.e. all but the masterimage):
inpath = input_info.resultfolder + '*YSO*/*_wcs.fits' 
imlist=glob.glob(inpath)
imlist.sort()
#imlist.remove(input_info.masterimage)

# make list of filenames for the output files:
outlist = deepcopy(imlist)
for i in np.arange(0, len(imlist)):
    outlist[i] = outlist[i].replace('_wcs.fits', '_wcs.mastercoo')

# make coordinate files for each individual image which will be used for the psf photometry of each image.
for i in np.arange(0, len(imlist)):
    print i
    photometry_wcs.make_mastercoos_for_images(input_info.masterimage, input_info.masterregfile, imlist[i], outlist[i], fantasymag, ds9=True)



