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


# get all image names for psf photometry:
inpath = input_info.resultfolder + '*YSO*/*_wcs.fits' 
imlist=glob.glob(inpath)
imlist.sort()
badlist = photometry_both.make_badlist(imlist)
if input_info.bad_files_exist == "yes":
    for b in input_info.bad_exposures: badlist.append(b)

imlist = photometry_both.sanitycheck(imlist, badlist)

# get list of all psf star files:
pstlist = deepcopy(imlist)
for i in np.arange(0, len(imlist)):
    pstlist[i] = imlist[i].replace('_wcs.fits', '_wcs.fits.pstbyhand')


photometry.do_psf_photometry_with_coo(imlist[:], satmag=input_info.satmag,photfilesuffix=input_info.photfilesuffix, psfstarlist=pstlist[:], psfcleaningradius=input_info.pfscleaningradius)









#### debugging remainders:

#imagebi = imlist[3]
#psfcleaningradius=input_info.pfscleaningradius
#photfilesuffix=input_info.photfilesuffix
#satmag=20.
#psfstarfile=pstlist[3]

#photometry.remove_previous_psfphot_files(imagebi)

##photometry.psf_star_fitting_2runs(imagebi, psfstarfile,photfilesuffix,satmag, psfcleaningradius=10.)

#sky,skydev=photometry.get_sky(imagebi)
#iraf.datapars.datamin='INDEF'
#iraf.datapars.sigma=skydev
## get inital list of psf stars (either selected by hand or let iraf select) and make a first-round psf model:
#original_psfrad=iraf.daopars.psfrad
#iraf.daopars.psfrad = 15.

#iraf.daophot.psf(imagebi,photfile=imagebi+photfilesuffix,pstfile=imagebi+'.pstbyhand',psfimage='default',opstfile='default',groupfile='default', interactive=False,verify=False,nclean=10)

#photometry.psf_star_checking(photometry.get_last_iraf(imagebi+'.pst'), satmag)


##use a smaller radius for nstar and substar to substract the "core only" of close neighbours
## find neighbors of psf stars which are in the fitting radius.
#iraf.nstar(imagebi,photometry.get_last_iraf(imagebi+'.psg'),'default','default','default', verbose = False)
## subtract the first-round psf model of those neighboring stars from the psf stars (so that wings of psf stars are star-free).
#iraf.substar(imagebi,photometry.get_last_iraf(imagebi+'.nst'), photometry.get_last_iraf(imagebi+'.pst'), photometry.get_last_iraf(imagebi+'.psf',suffix='fits'), 'default', verbose=True, Stdout=1)
##then reset psfrad to original value
#iraf.daopars.psfrad = original_psfrad
## now get the second-round psf model by fitting the neighbor-subtracted psf stars.
#iraf.daophot.psf(photometry.get_last_iraf(imagebi+'.sub',suffix='fits'), imagebi+photfilesuffix,photometry.get_last_iraf(imagebi+'.pst'), photometry.get_next_iraf(imagebi+'.psf',suffix='fits'), photometry.get_next_iraf(imagebi+'.pst'), photometry.get_next_iraf(imagebi+'.psg'), interactive=False,verify=False,nclean=10)


#print iraf.daopars.psfrad






