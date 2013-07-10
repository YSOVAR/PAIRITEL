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

iraf.datapars.fwhmpsf = input_info.fwhm_psf
iraf.fitskypars.annulus = input_info.radius_annulus
iraf.fitskypars.dannulus = input_info.width_annulus
iraf.photpars.apertures = input_info.fwhm_apertures
iraf.daopars.function = input_info.psf_function
iraf.daopars.varorder = input_info.psf_varorder
iraf.daopars.fitrad = input_info.psf_fitrad
iraf.daopars.psfrad = input_info.psfrad
iraf.daopars.fitsky = input_info.fitsky
iraf.daopars.sannulus = input_info.skyannulus
iraf.daopars.wsannulus = input_info.width_skyannulus
iraf.daopars.groupsky = input_info.groupsky
iraf.fitskypars.skyvalue = input_info.skyvalue

