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

if input_info.bad_files_exist == "yes":
    for b in input_info.bad_exposures: badlist.append(b)

imlist = photometry_both.sanitycheck(imlist, badlist)

# create the filenames for the individual PSF star files
outlist = deepcopy(imlist)
for i in np.arange(0,len(outlist)):
    outlist[i] = outlist[i] + '.pstbyhand'

# get a list of the results of the aperture photometry (.mag.1 files):
maglist = deepcopy(imlist)
for i in np.arange(0, len(maglist)): maglist[i] = imlist[i].replace('_wcs.fits', '_wcs.fits'+input_info.photfilesuffix)

n_matched = np.ones(len(imlist),int)*-99999.
# now do the actual conversion of the psf star list to image coordinates for each observation.
for i in np.arange(0, len(imlist)):
#for i in np.arange(0, 1):
    n_matched[i] = photometry_wcs.make_pst_file_from_ds9reg(input_info.psfstarfile, maglist[i], imlist[i], outlist[i])
    plt.close()


print('Include these files in your "bad_exposures": ')
for n in np.where(n_matched == 0)[0]:
    print(str(imlist[n]))



