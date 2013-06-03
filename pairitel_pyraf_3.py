# -*- coding: utf-8 -*-
#write this as script for now, turn parts in sensible functions later 
#keep all intermediate steps so I can track down errors
#change this behavious later e.g. use os.tmpfile or maybe tempfile.XXX()

# /usr/local/bin/ipython --colors lightbg # use this python because it knows where pyraf is.

#import matplotlib.pyplot as plt
import pyfits
import shutil
from pyraf import iraf
import glob
import os
import sys
import asciitable
#sys.path.append("/data/hguenther/Dropbox/code/python")
#import atpyextensions
#sys.path.append("/data/hguenther/Dropbox/code/python/atpy")
#import atpy
#from atpyextensions import catalog
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


# get all image names for psf photometry:
inpath = input_info.resultfolder + '*YSO*/*_wcs.fits' 
imlist=glob.glob(inpath)
imlist.sort()
#imlist.remove(input_info.masterimage)

# get list of all psf star files:
pstpath = input_info.resultfolder + '*YSO*/*_wcs.fits.pstbyhand' 
pstlist=glob.glob(pstpath)
pstlist.sort()
#pstlist.remove(input_info.masterimage + '.pstbyhand')

photometry.do_psf_photometry_with_coo(imlist[:], satmag=input_info.satmag, psfstarlist=pstlist[:])

#print iraf.daopars.psfrad


