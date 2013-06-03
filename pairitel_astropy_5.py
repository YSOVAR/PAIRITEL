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
# YSOVAR python scripts need to be in the PYTHONPATH for the following line:
from great_circle_dist import dist_radec, dist_radec_fast

import urllib
import StringIO

import input_info
reload(input_info)
import photometry_wcs
reload(photometry_wcs)
import photometry_both
reload(photometry_both)


# uncomment this to visually inspect the positions of the sources found for every nights.

## find files with the calibrated magnitudes:
#filepath =  input_info.resultfolder + '*YSO*/calibratedmags*.dat' 
#filelist = glob.glob(filepath)
#filelist.sort()

#gc = aplpy.FITSFigure(input_info.masterimage)
#gc.show_grayscale()

#print 'Overplotting positions of found sources of all nights and bands. This may take a while.'
#for i in np.arange(0, len(filelist)):
    #a = asciitable.read(filelist[i])
    #gc.show_markers(np.array(a['RA']), np.array(a['DEC']), marker = '+', edgecolor='r', facecolor='r')


#master = os.path.dirname(input_info.masterimage) + '/calibratedmags_' + os.path.basename(input_info.masterimage)[0] + '.dat'
#masterals = input_info.masterimage + '.als.1'
#b = asciitable.read(master)
#gc.show_markers(np.array(b['RA']), np.array(b['DEC']), marker = '+', edgecolor='b', facecolor='b')


# match sources into light curves. matching radius of 0.1 arcsec is good for this.
# for visual check, uncomment the stuff above - it overplots the positions of the detected sources in each night.
r_match = 0.1 * 1./3600.

bands = ['H', 'K', 'J']

(ra, dec) = photometry_wcs.coords_from_ds9(input_info.masterregfile)
p_id = np.arange(0, len(ra))

def make_lcs_from_disk():
    for b in bands[:]:
	a1=[]
	a2=[]
	a3=[]
	a4=[]
	a5=[]
	a6=[]
	a7=[]
	a8=[]
	a9=[]
	a10=[]
	names = ['Pair_ID', 'ra', 'dec', 't', 'H', 'H_err', 'J', 'J_err', 'K', 'K_err']
	tab = Table([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10], names=names)
	
	print b
	filepath =  input_info.resultfolder + '*YSO*/calibratedmags_' + b.lower() + '.dat' 
	filelist = glob.glob(filepath)
	filelist.sort()
	for i in p_id:
	    print i
	    for j in np.arange(0, len(filelist[:])):
		# find source in calibratedmags file of a given night.
		a = asciitable.read(filelist[j])
		# look only at the ones which are near in dec:
		ind_near = np.where(np.abs(a['DEC'] - dec[i]) <= r_match)[0]
		# do full distance calculation for those:
		if len(ind_near) > 0:
		    distance = np.sqrt((ra[i] - a[ind_near]['RA'])**2 + (dec[i] - a[ind_near]['DEC'])**2)
		    ind = np.where(distance == distance.min())[0][0]
		    tab.add_row(np.ones(len(names))*np.nan)
		    tab[-1]['Pair_ID'] = i
		    tab[-1]['ra'] = ra[i]
		    tab[-1]['dec'] = dec[i]
		    tab[-1]['t'] = a[0]['t'] - 2400000.5 # all lines in a given nightly file have the same time stamp.
		    tab[-1][b] = a[ind_near[ind]]['P' + b + '_cal']
		    tab[-1][b + '_err'] = a[ind_near[ind]]['P' + b + '_cal_err_corr']
	
	asciitable.write(tab, input_info.resultfolder + b + '_pairitel.csv' )


def make_lcs():
    for b in bands[:]:
	a1=[]
	a2=[]
	a3=[]
	a4=[]
	a5=[]
	a6=[]
	a7=[]
	a8=[]
	a9=[]
	a10=[]
	names = ['Pair_ID', 'ra', 'dec', 't', 'H', 'H_err', 'J', 'J_err', 'K', 'K_err']
	tab = Table([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10], names=names)
	print b
	for i in p_id:
	    print i
	
	    distance = dist_radec_fast(ra[i], dec[i], datadict['data' + b]['RA'], datadict['data' + b]['DEC'], scale = r_match, unit ='deg')
	    ind = np.where(distance < r_match)[0]
	    
	    for i1 in ind:
	        
	        tab.add_row(np.ones(len(names))*np.nan)
		tab[-1]['Pair_ID'] = i
		tab[-1]['ra'] = ra[i]
		tab[-1]['dec'] = dec[i]
		tab[-1]['t'] = datadict['data' + b]['t'][i1] - 2400000.5 # all lines in a given nightly file have the same time stamp.
		tab[-1][b] = datadict['data' + b]['P' + b + '_cal'][i1]
		tab[-1][b + '_err'] = datadict['data' + b]['P' + b + '_cal_err_corr'][i1]
	    
	asciitable.write(tab, input_info.resultfolder + b + '_pairitel.csv' )

# if pairitel_astropy_5.oy has been run before, use this:
# (it collects the calibrated magnitudes of each night into light curves for each source)
make_lcs()

# if something has chrashed and the calibrated data needs to be read from hard drive, run this:
# (it does the same, but takes longer)
#make_lcs_from_disk()
