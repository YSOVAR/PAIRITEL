# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import aplpy
import astropy.io.fits as pyfits
import astropy
from astropy.io import ascii
from astropy.io import fits
from astropy.wcs import WCS
from copy import deepcopy
import glob
import os
from scipy import optimize
from astropy.table import Table, Column
import YSOVAR
from YSOVAR import atlas
from YSOVAR.great_circle_dist import dist_radec, dist_radec_fast

import urllib
import StringIO

import input_info
reload(input_info)
import photometry_wcs
reload(photometry_wcs)
import photometry_both
reload(photometry_both)



def check_band(filename):
    checkband = os.path.basename(filename)[0].lower()
    return checkband



def make_lcs_from_nights(filenames, r_match=0.1*1./3600., masterregfile=input_info.masterregfile, outname=input_info.resultfolder + 'rawlcs.dat'):
    # match sources into light curves. matching radius of 0.1 arcsec is good for this.
    # first get coordinates of sources in masterregfile:
    (ra, dec) = photometry_wcs.coords_from_ds9(input_info.masterregfile)
    p_id = np.arange(0, len(ra))
    # names for the columns of the output table:
    names = ['PID', 'ra', 'de', 't', 'filter', 'mag_raw', 'mag_raw_err', 'mag_cal1', 'mag_cal1_err', 'mag_cal2', 'mag_cal2_err']
    dtypes = ['<i8', '<f8', '<f8', '<f8', '|S1', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8']
    dummy = [[]]
    for i in np.arange(0, len(names)-1):
	dummy.append([])
    
    tab = Table(dummy, names=names, dtypes=dtypes)
    
    for filename in filenames:
	# get corresponding imagefile to .als.1 file:
	imagefile = filename[0:-6]
	
	pairitel = ascii.read(filename, Reader=ascii.Daophot)
	pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='RA'))
	pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='DEC'))

	checkband = check_band(filename)
	# translate X,Y image coordinates into wcs:
	hdus = fits.open(imagefile)
	wcs = WCS(hdus[0].header)
	pairitel['RA'], pairitel['DEC'] = wcs.wcs_pix2world(pairitel['XCENTER'], pairitel['YCENTER'], 1)
	t = np.float(hdus[0].header['HJULDATE'])
	hdus.close()
	pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*t, name='t'))
	
	# match nighlty source detections to positions in masterregfile.
	
	for i in np.arange(0, len(p_id)):
	    # matching with masterregfile:
	    dist = dist_radec_fast(ra[i], dec[i], pairitel['RA'], pairitel['DEC'], scale = r_match, unit ='deg')
	    if dist.min() < r_match:
		ind = np.where(dist == dist.min())[0]
		tab.add_row(np.ones(len(names))*-99999.)
		tab[-1]['PID'] = i
		tab[-1]['ra'] = ra[i]
		tab[-1]['de'] = dec[i]
		tab[-1]['t'] = pairitel[ind]['t'] - 2400000.5
		tab[-1]['mag_raw'] = pairitel[ind]['MAG']
		tab[-1]['mag_raw_err'] = pairitel[ind]['MERR']
		tab[-1]['filter'] = str(checkband)
	
    ascii.write(tab, outname, delimiter=',')
    return tab


def calc_rchisq(values, fitvalues, errors):
    # calculated the reduced chi squared value for data to a fit.
    rchisq = np.sum( ((values - fitvalues)/(errors))**2. )/len(values)
    return rchisq





def fit_2mass(pmags, perrs, twomassmags, fig, withplot=True):
    # now do the fit:
    # define our (line) fitting function (with slope fixed at 1):
    fitfunc = lambda p, x: p[0] + x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    pinit = [0.] # initital guess for intercept

    p = optimize.leastsq(errfunc, pinit, args=(twomassmags, pmags, perrs), full_output=1)[0]

    chisq = calc_rchisq(pmags - p[0], twomassmags, perrs)
    fiterror = np.std((pmags - p[0]) - twomassmags)

    if withplot == True:
	fig.clf()
        ax = fig.add_subplot(111)
	ax.plot(twomassmags, pmags, '.')
	x = np.linspace(np.min(twomassmags), np.max(twomassmags), 100)
	ax.plot(x, p[0] + x, 'b-')
	ax.text(0.3, 0.7, 'red. chi**2 = ' + str(np.round(chisq,1)), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
	plt.show()

    return (p[0], fiterror, chisq)


def do_calibration_for_subset(ptl, t, ind_subset, mag_in_name, err_in_name, mag_out_name, err_out_name):
    
    chisqs = np.array([])
    fiterrors = np.array([])
    f1 = plt.figure()
    
    for i in np.arange(0, len(t)):
	# find entries in ptl which are from a given time stamp: (filters have different time stamps in the same night)
	i_night = np.where( ptl['t'] == t[i])[0]
	#print i_night
	filt = list(set(ptl[i_night]['filter']))
	print filt
	#if len(filt) > 1:
	    #print "There's something wrong with your timestamps."
	
	for f in filt:
	    twomassmag = '2' + f.upper()
	    
	    # take only the constant sources of that night for the fit:
	    i_fit = np.intersect1d(np.array(ind_subset, dtype=int), i_night)
	    i_fit = i_fit[ptl[i_fit]['filter'] == f]
	    
	    # compute the fit
	    (par, fiterror, chisq) = fit_2mass(ptl[i_fit][mag_in_name], ptl[i_fit][err_in_name], ptl[i_fit][twomassmag], fig=f1)
	    
	    # save corrected magnitudes into ptl (for all sources of that night):
	    i_all = i_night[ptl[i_night]['filter'] == f]
	    for i in i_all:
	        ptl[i][mag_out_name] = ptl[i][mag_in_name] - par
	        ptl[i][err_out_name] = np.sqrt(ptl[i][err_in_name]**2 + fiterror**2)
	    
	    chisqs = np.append(chisqs, chisq)
	    fiterrors = np.append(fiterrors, fiterror)
	
    return (ptl, chisqs, fiterrors)



