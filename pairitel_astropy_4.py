# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import astropy
import astropy.io.fits as pyfits
import astropy.io.ascii as ascii
from copy import deepcopy
import glob
import os
from scipy import optimize
from astropy.table import Table, Column
import sys
import YSOVAR
from YSOVAR import atlas
from YSOVAR.great_circle_dist import dist_radec, dist_radec_fast
from astropy.wcs import WCS
from astropy.io import fits


import urllib
import StringIO

import input_info
reload(input_info)
import photometry_wcs
reload(photometry_wcs)
import photometry_both
reload(photometry_both)
import calibration
reload(calibration)



# find all .als.1 files:
filepath =  input_info.resultfolder + '*YSO*/*.als.1' 
filelist = glob.glob(filepath)
filelist.sort()

# collect all raw magnitudes and their errors into lightcurves, matched onto the positions in the masterreg file. Automatically writes this to disk.
ptl = calibration.make_lcs_from_nights(filelist, r_match=0.1*1./3600., masterregfile=input_info.masterregfile, outname=input_info.resultfolder + 'rawlcs.dat')
