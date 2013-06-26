# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import astropy.io.fits as pyfits
import astropy
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


# now take only the sources whose light curves are pretty constant and redo the calibration to 2MASS.
threshold = input_info.threshold_lc
ind_const = mycloud[mycloud['stddev_Hcal1']<threshold]['PID']
i_const = np.array([])
for i in ind_const:
    ind = np.where(ptl['PID'] == i)[0]
    i_const = np.append(i_const, ind)

i_const_2mass = np.intersect1d( i_const, np.where( (ptl['2H']>0) & (ptl['2J']>0) & (ptl['2K']>0) )[0])

(ptl, chisq2, fiterror2) = calibration.do_calibration_for_subset(ptl, times, i_const_2mass, 'mag_raw', 'mag_raw_err', 'mag_cal2', 'mag_cal2_err')

# plot the old and the new scatter introduced by the fit.
plt.clf()
plt.hist(fiterror, bins=20, alpha=0.5)
plt.hist(fiterror2, bins=20, alpha=0.5)
plt.xlabel('mag')
plt.legend(['old fit scatter','new fit scatter'])

# save calibrated light curves in YSOVAR2-like format:
print('Writing results to disk...')
ascii.write(ptl, input_info.resultfolder + 'cal2lcs.dat')
for b in ['j', 'h', 'k']:
    ind = np.where(ptl['filter'] == b)[0]
    ascii.write(ptl[ind], input_info.resultfolder + b + '_cal2lcs.dat')

# add calibrated light curves to mycloud.
for b in ['j', 'h', 'k']:
    lc = ascii.read(input_info.resultfolder + b + '_cal2lcs.dat')
    # find out which columns to use
    cross_ids = atlas.makecrossids_all(mycloud, lc, 1./3600., ra1 = 'ra', dec1 = 'dec', ra2 = 'ra', dec2 = 'de')
    mycloud.add_mags(lc, cross_ids, ['mag_cal2','mag_cal2_err','t'], b.upper() + 'cal2')


# sort all final calibrated lcs by time (right now the time coordinate is unsorted):
print('sorting lcs by time...')
for i in np.arange(0, len(mycloud)):
    for b in ['J', 'H', 'K']:
        if len(mycloud[i].lclist[0]['t'+b+'cal2']) > 0:
	    stuff = zip(mycloud[i].lclist[0]['t'+b+'cal2'], mycloud[i].lclist[0]['m'+b+'cal2'], mycloud[i].lclist[0]['m'+b+'cal2_error']) # sort by first argument in zip.
	    stuff.sort()
	    (mycloud[i].lclist[0]['t'+b+'cal2'], mycloud[i].lclist[0]['m'+b+'cal2'], mycloud[i].lclist[0]['m'+b+'cal2_error']) = zip(*stuff)



# plot a source for illustrative purposes:

ind = np.where( (mycloud['mean_Hcal2'] < 15))[0]

if len(ind) > 0:
    ind = random.sample(ind, 1)[0]
    plt.clf()
    plt.errorbar(mycloud[ind].lclist[0]['tJcal2'], mycloud[ind].lclist[0]['mJcal2'], yerr = mycloud[ind].lclist[0]['mJcal2_error'], fmt='.-', capsize=0 )
    plt.errorbar(mycloud[ind].lclist[0]['tHcal2'], mycloud[ind].lclist[0]['mHcal2'], yerr = mycloud[ind].lclist[0]['mHcal2_error'], fmt='.-', capsize=0  )
    plt.errorbar(mycloud[ind].lclist[0]['tKcal2'], mycloud[ind].lclist[0]['mKcal2'], yerr = mycloud[ind].lclist[0]['mKcal2_error'], fmt='.-', capsize=0  )
    plt.title('a sample source from your cluster (' + str(ind) + ')')
    plt.xlabel('t (MJD)')
    plt.ylabel('mag')
    plt.legend(['J', 'H', 'K'], 'best')








