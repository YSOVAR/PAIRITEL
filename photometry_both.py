# -*- coding: utf-8 -*-
# this containes routines that are used by both the pyraf and the astropy scripts.

import glob
import pyfits
import numpy as np
import input_info

def make_badlist(imlist):
    badlist = []
    for n in imlist:
	hdulist = pyfits.open(n)
	ra = hdulist[0].header['CRVAL1']
	dec = hdulist[0].header['CRVAL2']
	if ( (np.abs(float(ra) - input_info.ra_cluster) > input_info.width_images) | (np.abs(float(dec) - input_info.dec_cluster) > input_info.width_images) ):
	    # check if it is already included in badlist:
	    if not(n.replace('_coadd.fits', '').replace('_coadd.weight.fits','') in badlist):
	        # save only beginning of filename in badlist (suffix may change for other file types.)
	        badlist.append(n.replace('_coadd.fits', '').replace('_coadd.weight.fits',''))
    
    return badlist

def sanitycheck(filelist, badlist):
    newlist = []
    for f in filelist:
        test=0
        for b in badlist:
            if (b in f):
                test=1
        if test==0:
            newlist.append(f)
    return newlist

