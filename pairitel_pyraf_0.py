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


# get list of all long exposure files, both actual target observations and the weight files.
datapath = input_info.rawdatafolder + '*YSO*/*_long_*coadd*fits'
imlist = glob.glob(datapath)
imlist.sort()

# this finds the observations which are not actually centered on the cluster.
badlist = photometry_both.make_badlist(imlist)

# sanitycheck excludes observations which do not include the cluster from the analysis file list.
imlist = photometry_both.sanitycheck(imlist, badlist)


# this copies all data which will be used in the analysis to the resultfolder (so that you have a copy of the unaltered data in the original folder), and then, in the resultfolder, trims the images and normalizes them by the exposure time mask.
# also adds a header keyword for the readout noise (is missing in the original files).
readoutnoise = 10.
photometry.prepare_files(imlist, input_info.resultfolder, readoutnoise)

# set various parameters for the Pairitel observations
photometry.set_Pairitel_params()

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




# now trim the images so that only parts with good exposure are retained:
trimpath = input_info.resultfolder + '*YSO*/*_wcs_unnormed.fits'
trimlist = glob.glob(trimpath)
trimlist.sort()
masklist = deepcopy(trimlist)
for i in np.arange(0, len(masklist)): masklist[i] = masklist[i].replace('_wcs_unnormed.fits', '_coadd.weight.fits')

print "trimming..."
for i in np.arange(0, len(masklist)):
#for i in np.arange(0, 1):
    photometry.apply_rectangle_mask(trimlist[i], masklist[i], input_info.threshold, input_info.min_width, input_info.min_height)



print "normalizing mask..."
for i in np.arange(0, len(masklist)):
#for i in np.arange(0, 1):
    photometry.mask_norm(masklist[i].replace('_coadd.weight.fits', '_coadd.weight_trimmed.fits'))


#iraf.cd('/swiper.real/kpoppen/IR/GD1215/output/YSO.47.1')
#os.remove('test.fits')
#iraf.imcopy('h_long_YSO.47.1_coadd.weight.fits'+'[100:200,100:200]', 'test.fits')
#iraf.imarith('test.fits','/', 2., 'test2.fits')


print "applying mask..."
for i in np.arange(0, len(masklist)):
    photometry.divide_by_mask(trimlist[i].replace('.fits', '_trimmed.fits'), masklist[i].replace('.fits', '_normed.fits'))


#datapath_sky = input_info.resultfolder + '*YSO*/*_wcs.fits' 
#datalist = glob.glob(datapath_sky)
#datalist.sort()

#iraf.centerpars.calgorithm = "centroid"
#iraf.centerpars.cbox = 30.
#iraf.centerpars.maxshift = 20.

#photometry.correct_coordinates(datalist[:], radius=10., twomassmag = 16.)

#matchout=iraf.msccmatch('h_long_YSO.47.1_wcs.fits','2mass_YSO.47.1_h_16.cat',interactive=False, Stdout=1,cfrac=0.95,rms=5, maxshif=1000., nsearch=100, nfit=3, update='yes')

#iraf.cd('/swiper.real/kpoppen/IR/GD1215/output/YSO.48.10')

#matchout=iraf.msccmatch('h_long_YSO.48.10_wcs.fits','2mass_YSO.48.10_h_16.cat',interactive=False, Stdout=1,cfrac=0.95,rms=10, maxshif=1000., nsearch=100, nfit=3, update='yes')


