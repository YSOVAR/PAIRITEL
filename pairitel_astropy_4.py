# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import pyfits
import astropy
import astropy.io.ascii as asciitable
from copy import deepcopy
import glob
import os
from scipy import optimize
from astropy.table import Table, Column

import urllib
import StringIO

import input_info
reload(input_info)
import photometry_wcs
reload(photometry_wcs)
import photometry_both
reload(photometry_both)


# load 2mass catalogue from disk:
data2mass = asciitable.read(input_info.resultfolder + '2mass_for_cal.dat')

# find all .als.* files - this should only be one file per image, except for the masterimage which has als.1, als.2, and als.3:
filepath =  input_info.resultfolder + '*YSO*/*.als.1' 
filelist = glob.glob(filepath)
filelist.sort()
# get rid of the two old als files of the masterimage:
#filelist.remove(input_info.masterimage + '.als.1')
#filelist.remove(input_info.masterimage + '.als.2')

f = plt.figure()
chisq = np.ones(len(filelist))*-99999.


# now transform found pairitel sources into wcs, find matching 2mass sources, and calibrate pairitel against 2mass. I use nominal errors of 0.1 mag for all pairitel sources instead of the output errors from the psf photometry. This is because sometimes there is large variability in a bright source (which has a small photometric error), which then would mess up the whole line fitting because that data point would have so much weight.
for filename in filelist[:]:
    
    imagefile = filename[0:-6]
    
    pairitel = asciitable.read(filename, Reader=asciitable.Daophot)
    checkband = os.path.basename(filename)[0]
    if checkband == 'h':
	(band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('2H', 'PH', 'PH_err', 'PH_cal', 'PH_cal_err')
    else:
	if checkband == 'k':
	    (band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('2K', 'PK', 'PK_err', 'PK_cal', 'PK_cal_err')
	else:
	    if checkband == 'j':
		(band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('2J', 'PJ', 'PJ_err', 'PJ_cal', 'PJ_cal_err')
    
    print checkband
    
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='RA'))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='DEC'))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mp))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpe))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2m))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpc))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpce))
    
    pairitel[band2mp] = pairitel['MAG']
    pairitel[band2mpe] = pairitel['MERR']
    
    # translate X,Y image coordinates into wcs:
    gc = aplpy.FITSFigure(imagefile)
    for i in np.arange(0,len(pairitel)):
	(pairitel[i]['RA'], pairitel[i]['DEC']) = gc.pixel2world(pairitel[i]['XCENTER'], pairitel[i]['YCENTER'])
    
    plt.close()
    hdulist = pyfits.open(imagefile)
    t = np.float(hdulist[0].header['HJULDATE'])
    hdulist.close()
    
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*t, name='t'))
    
    
    # find best-matching 2mass source for each pairitel source; if distance small enough, add the 2mass magnitude to the pairitel table.
    
    matchradius = 1./3600.
    
    for i in np.arange(0, len(pairitel)):
	dist = np.sqrt((data2mass['ra'] - pairitel[i]['RA'])**2 + (data2mass['dec'] - pairitel[i]['DEC'])**2)
	mindist = dist.min()
	if mindist < matchradius:
	    ind = np.where(dist == mindist)[0]
	    pairitel[i][band2m] = data2mass[ind][band2m]
    
    print np.where(pairitel[band2m]>-99999.)
    
    good = np.where(pairitel[band2m]>-99999.)[0]
    
    f.clf()
    ax = f.add_subplot(111)
    ax.plot(pairitel[good][band2m], pairitel[good][band2mp], '.')
    #plt.show()
    
    # fit straight line to that.
    # define our (line) fitting function:
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    
    #error = pairitel[good][band2mpe]
    error = 0.1 * np.ones(len(good))
    
    pinit = [1.0, -1.0]
    out = optimize.leastsq(errfunc, pinit, args=(pairitel[good][band2m], pairitel[good][band2mp], error), full_output=1)
    pfinal = out[0]
    x = np.linspace(np.min(pairitel[good][band2m]), np.max(pairitel[good][band2m]), 100)
    ax.plot(x, pfinal[0] + pfinal[1]*x)
    
    # calculate new pairitel magnitudes.
    pairitel[band2mpc] = (pairitel[band2mp] - pfinal[0])/pfinal[1]
    pairitel[band2mpce] = (pairitel[band2mpe])/pfinal[1]
    
    print pairitel[good][band2mpc]
    print pairitel[good][band2mp]
    print pairitel[good][band2m]
    
    n = np.where(filename == np.array(filelist))[0][0]
    #chisq[n] = np.sum(((pairitel[good][band2mp] - (pfinal[0] + pfinal[1]*pairitel[good][band2m]))/(pairitel[good][band2mpe]))**2.)/len(pairitel[good])
    chisq[n] = np.sum(((pairitel[good][band2mp] - (pfinal[0] + pfinal[1]*pairitel[good][band2m]))/(error))**2.)/len(pairitel[good])
    
    ax.text(0.3, 0.7, 'red. chi**2 = ' + str(np.round(chisq[n],1)), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    plt.show()
    
    # write results to file.
    asciitable.write(pairitel, os.path.dirname(filename) + '/calibratedmags_' + checkband + '.dat')

# if necessary, do some manual checking for high chisquared values.
print np.where(chisq > 5)[0]
plt.figure()
plt.hist(chisq, bins=100)
plt.title('Distribution of reduced chi**2 values')


