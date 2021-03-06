# -*- coding: utf-8 -*-
#write this as script for now, turn parts in sensible functions later 
#keep all intermediate steps so I can track down errors
#change this behavious later e.g. use os.tmpfile or maybe tempfile.XXX()

# /usr/local/bin/ipython --colors lightbg # use this python because it knows where pyraf is.

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

iraf.noao()
iraf.digiphot()
iraf.daophot()


# do psf photometry for master image.
masterpsffile = input_info.masterimage + '.pstbyhand'

threshold_secondrun = 10.

photometry.do_psf_photometry([input_info.masterimage], satmag=input_info.satmag, photfilesuffix=input_info.photfilesuffix,psfstarfile=masterpsffile, thresh = threshold_secondrun, psfcleaningradius=input_info.pfscleaningradius)





