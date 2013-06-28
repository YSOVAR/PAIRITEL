# -*- coding: utf-8 -*-
# start /usr/local/bin/ipython --colors lightbg # use this python because it knows where pyraf is.

import astropy.io.fits as pyfits
import shutil
from pyraf import iraf
import glob
import os
import sys
import astropy.io.ascii as ascii
import numpy as np
import string
import pickle
from copy import deepcopy

import input_info
reload(input_info)
sys.path.append(input_info.pairitel_scripts_path)
import photometry
reload(photometry)
import photometry_both
reload(photometry_both)

iraf.noao()
iraf.digiphot()
iraf.daophot()


# get a list of the images and copy them into new files on which the coordinate correction will be performed.
datapath_sky = input_info.resultfolder + '*YSO*/*_coadd.fits' 
datalist = glob.glob(datapath_sky)
datalist.sort()
for f in datalist:
    shutil.copy(f,f.replace('_coadd.fits', '_wcs_unnormed.fits'))


# correct coordinate shifts by comparing each image to the 2MASS catalogue and applying necessary coordinate corrections.

wcspath = input_info.resultfolder + '*YSO*/*_wcs_unnormed.fits' 
wcslist = glob.glob(wcspath)
wcslist.sort()
# use different centering parameters for that; set back to original values afterwards.
iraf.centerpars.calgorithm = "centroid"
iraf.centerpars.cbox = 10.
iraf.centerpars.maxshift = 10.
# correct coordinates; do this three times for best results.
photometry.correct_coordinates(wcslist[:], radius=10., twomassmag = 17.)
photometry.correct_coordinates(wcslist[:], radius=10., twomassmag = 17.)
photometry.correct_coordinates(wcslist[:], radius=10., twomassmag = 17.)
# reset centering parameters to original values.
iraf.centerpars.cbox = 9.
iraf.centerpars.maxshift = 1.

