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




# Read the uncalibrated light curves from disk:
ptl = ascii.read(input_info.resultfolder + 'rawlcs.dat')

# now transform this into a ysovar atlas object.
stars = atlas.dict_from_csv(csvfile=input_info.resultfolder + 'rawlcs.dat',  match_dist = 0., min_number_of_times = 5, channels = {'j': 'J', 'h': 'H', 'k': 'K'}, data = [], floor_error = {'j': 0., 'h': 0., 'k': 0.}, mag = 'mag_raw', emag = 'mag_raw_err', time = 't', bg = None, source_name = 'PID',  verbose = True, readra='ra', readdec='de', sourceid='PID', channelcolumn='filter')

atlas.check_dataset(stars)

mycloud = atlas.YSOVAR_atlas(lclist = stars)

# add the pairitel id as float to the atlas:
mycloud.add_catalog_data(ptl, names = ['PID'], ra1='ra', dec1='dec', ra2='ra', dec2='de')

# now find matching 2MASS sources for that.
# load 2mass catalogue from disk:
data2mass = ascii.read(input_info.resultfolder + '2mass_for_cal.dat')
data2mass.rename_column('ra', 'ra2mass')
data2mass.rename_column('dec', 'dec2mass')
#twomass_names = ['2massid', 'ra2mass', 'dec2mass', '2J', '2H', '2K']
twomass_names = ['2J', '2H', '2K']

# add 2MASS magnitudes to the atlas object
mycloud.add_catalog_data(data2mass, names = twomass_names, ra1='ra', dec1='dec', ra2='ra2mass', dec2='dec2mass')

# also add 2mass magnitudes to the ptl array - this makes the light curve calibration easier than doing the calibration directly with the atlas object.
ptl.add_column(Column(data = np.ones(len(ptl['PID']))*-99999., name='2J'))
ptl.add_column(Column(data = np.ones(len(ptl['PID']))*-99999., name='2H'))
ptl.add_column(Column(data = np.ones(len(ptl['PID']))*-99999., name='2K'))
for i in np.arange(0, len(ptl)):
    ind = np.where( mycloud['PID'] == ptl[i]['PID'] )[0]
    if mycloud['2J'][ind]> 0: ptl[i]['2J'] = mycloud['2J'][ind]
    if mycloud['2H'][ind]> 0: ptl[i]['2H'] = mycloud['2H'][ind]
    if mycloud['2K'][ind]> 0: ptl[i]['2K'] = mycloud['2K'][ind]


## read additional data on YSO classes from catalog if catalog supplied:
#if input_info.catalogfile != '':
    #cat = ascii.read(input_info.catalogfile)
    ## take only the data part from the masked array (this is a workaround because there's a bug with renaming columns for masked arrays in the current astropy version)
    #mycloud.add_catalog_data(cat, names = ['col7'], ra1='ra', dec1='dec', ra2='col2', dec2='col3')
    #print np.where((mycloud['col7'] < 1) & (mycloud['2H']>0) )


# fit raw Pairitel magnitudes to 2mass magnitudes using a subset of stars. This set can be:
# A) all stars with 2mass detections (quite large scatter; systematic error by fit will overestimate real error because there'll be real variables in the fit)
# B) 2mass detected stars which have quite constant light curves and thus serve as "empiric non-YSOs". For this to work, one first needs to do a rough calibration using all stars in order to get calibrated light curves (albeit with too large errors).

times = np.array(list(set(ptl['t']))) # detour over the list, because sets cannot be transformed directly into np.arrays

# take all sources to fit (in this case all with 2mass data)
i_ptlsubset = np.where( (ptl['2H']>0) & (ptl['2J']>0) & (ptl['2K']>0) )[0]


(ptl, chisq, fiterror) = calibration.do_calibration_for_subset(ptl, times, i_ptlsubset, 'mag_raw', 'mag_raw_err', 'mag_cal1', 'mag_cal1_err')

print('Writing results to disk...')
ascii.write(ptl, input_info.resultfolder + 'cal1lcs.dat')
for b in ['j', 'h', 'k']:
    ind = np.where(ptl['filter'] == b)[0]
    ascii.write(ptl[ind], input_info.resultfolder + b + '_cal1lcs.dat')

for b in ['j', 'h', 'k']:
    lc = ascii.read(input_info.resultfolder + b + '_cal1lcs.dat')
    # find out which columns to use
    cross_ids = atlas.makecrossids_all(mycloud, lc, 1./3600., ra1 = 'ra', dec1 = 'dec', ra2 = 'ra', dec2 = 'de')
    mycloud.add_mags(lc, cross_ids, ['mag_cal1','mag_cal1_err','t'], b.upper() + 'cal1')


# calculate mean of mycloud:
print('calculating standard deviation in preliminary light curves...')
print mycloud['stddev_Hcal1']

plt.clf()
plt.hist(mycloud['stddev_Hcal1'], bins=50)
plt.xlabel('standard deviation of H band light curves')


