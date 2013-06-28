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

# apply the coordinate transformations done for the individual images also to their respective masks. This is not strictly necessary, because the mask trimming and dividing is done on a pixel-by-pixel basis, but it's nicer this way when you want to look at the mask in ds9 and compare things by eye.
trimpath = input_info.resultfolder + '*YSO*/*_wcs_unnormed.fits'
trimlist = glob.glob(trimpath)
trimlist.sort()
masklist = deepcopy(trimlist)
for i in np.arange(0, len(masklist)): masklist[i] = masklist[i].replace('_wcs_unnormed.fits', '_coadd.weight.fits')

for f in masklist:
    image = f.replace('_coadd.weight.fits','_wcs_unnormed.fits')
    for headerkeyword in ['RADECSYS', 'CTYPE1', 'CRVAL1', 'CRPIX1', 'CD1_1', 'CD1_2', 'CTYPE2', 'CRVAL2', 'CRPIX2', 'CD2_1', 'CD2_2']:
        iraf.imgets(image, headerkeyword)
        iraf.hedit(images=f, fields=headerkeyword, value=iraf.imgets.value, verify="no")



# before the trimming and normalizing, delete the respective files which may have been produced in earlier runs. Otherwise IRAF will refuse to overwrite the files.
for fl in glob.glob(input_info.resultfolder + '*_trimmed.fits*'):
    os.remove(fl)

for fl in glob.glob(input_info.resultfolder + '*_normed.fits*'):
    os.remove(fl)

for fl in glob.glob(input_info.resultfolder + '*_wcs.fits*'):
    os.remove(fl)


# now trim the images so that only parts with good exposure are retained:
print "trimming..."
for i in np.arange(0, len(masklist)):
#for i in np.arange(0, 1):
    photometry.apply_rectangle_mask(trimlist[i], masklist[i], input_info.threshold, input_info.min_width, input_info.min_height)



print "normalizing mask..."
for i in np.arange(0, len(masklist)):
#for i in np.arange(0, 1):
    photometry.mask_norm(masklist[i].replace('_coadd.weight.fits', '_coadd.weight_trimmed.fits'))


print "applying mask..."
for i in np.arange(0, len(masklist)):
    photometry.divide_by_mask(trimlist[i].replace('.fits', '_trimmed.fits'), masklist[i].replace('.fits', '_normed.fits'))


