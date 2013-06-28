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


# now test visually in ds9 for images that have not been properly aligned.
# put info for those which have not been properly aligned into the manual wcs correction file.
if input_info.badwcs == "yes":
    corrbyhand = ascii.read(input_info.wcs_corrbyhand_file)
    for i in np.arange(0, len(corrbyhand)):
        # remove file if beyond saving
        if corrbyhand[i]["beyond_saving"] == "yes":
            print("removing files which are beyond saving...")
            if os.path.isfile(corrbyhand[i]["filename_incl_path"]): os.remove(corrbyhand[i]["filename_incl_path"])
        else:
            # if saveable, update header with manual coordinates.
            print("updating headers with manual coordinates...")
	    iraf.cd(os.path.dirname(corrbyhand[i]["filename_incl_path"]))
	    f = os.path.basename(corrbyhand[i]["filename_incl_path"])
	    iraf.hedit(f, "CRVAL1", corrbyhand[i]["ra_good"], verify="no")
	    iraf.hedit(f, "CRPIX1", corrbyhand[i]["x_pixel"], verify="no")
	    iraf.hedit(f, "CRVAL2", corrbyhand[i]["dec_good"], verify="no")
	    iraf.hedit(f, "CRPIX2", corrbyhand[i]["y_pixel"], verify="no")
            # run the automatic coordinate correction again:
            print("running refined coordinate matching for updated files...")
            iraf.centerpars.cbox = 10.
	    iraf.centerpars.maxshift = 10.
            photometry.correct_coordinates([corrbyhand[i]["filename_incl_path"]], radius=10., twomassmag = 17.)
            photometry.correct_coordinates([corrbyhand[i]["filename_incl_path"]], radius=10., twomassmag = 17.)
            photometry.correct_coordinates([corrbyhand[i]["filename_incl_path"]], radius=10., twomassmag = 17.)
            photometry.correct_coordinates([corrbyhand[i]["filename_incl_path"]], radius=10., twomassmag = 17.)
            iraf.centerpars.cbox = 9.
	    iraf.centerpars.maxshift = 1.


