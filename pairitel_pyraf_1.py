# -*- coding: utf-8 -*-
#write this as script for now, turn parts in sensible functions later 
#keep all intermediate steps so I can track down errors
#change this behavious later e.g. use os.tmpfile or maybe tempfile.XXX()

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
import photometry
reload(photometry)
import photometry_both
reload(photometry_both)


# get a list of the coordinate-corrected images.
datapath_extractsources = input_info.resultfolder + '*YSO*/*_wcs.fits' 
extractlist = glob.glob(datapath_extractsources)
extractlist.sort()

# perform aperture photometry on those images; found sources are also saved as a ds9 region file for checking by eye.
photometry.do_aperture_photometry(extractlist, input_info.photfilesuffix)

#check images for best sharpness.
datapath = input_info.resultfolder + '*YSO*/*_wcs.fits.coo.1'
datalist = glob.glob(datapath)
datalist.sort()

(sharpnesses, bestfiles, bestsharps, mostfiles, stars_found) = photometry.sort_by_apertphotresults(datalist, nbest=8)
 
#print 'These are the sharpest images.\n'
#for b in bestfiles: print b[0:-6]
#for s in bestsharps: print s

print 'These are the images in which the most stars were found.\n'
for b in mostfiles: print b[0:-6]
for s in stars_found: print s

print('Choose one as your master image for the source extraction (best look at them by hand in ds9 first). Then define your chosen master image in the input_info.py file.')


