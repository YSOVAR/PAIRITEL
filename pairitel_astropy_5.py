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

# collect values into light curves:

# first get master coordinates from masterregfile.


#(ra, dec) = photometry_wcs.coords_from_ds9(input_info.masterregfile)
#p_id = np.arange(0, len(ra))

## make one LC:
#filepath =  input_info.resultfolder + '*YSO*/calibratedmags_h.dat' 
#filelist = glob.glob(filepath)
#filelist.sort()
#imagepath =  input_info.resultfolder + '*YSO*/h_long*_wcs.fits' 
#imagelist = glob.glob(imagepath)
#imagelist.sort()


#b = asciitable.read(os.path.dirname(input_info.masterimage) + '/calibratedmags_h.dat')
#ind = np.where(b['PH_cal'] < 11.)[0][4]
#distance = np.sqrt((ra - b[ind]['RA'])**2 + (dec - b[ind]['DEC'])**2)
#ind_file = np.where(distance == distance.min())[0][0]

#lc = []
#t = []
#e = []
#n = ind_file
##n=550

#for i in np.arange(0, len(filelist)):
    #print i
    #a = asciitable.read(filelist[i])
    #distance = np.sqrt((ra[n] - a['RA'])**2 + (dec[n] - a['DEC'])**2)
    #if distance.min() < 0.1*1./3600.:
        #ind = np.where(distance == distance.min())[0][0]
        #lc.append(a[ind]['PH_cal'])
        #e.append(a[ind]['PH_cal_err'])
        #hdulist = pyfits.open(imagelist[i])
        #t.append(np.float(hdulist[0].header['HJULDATE']))
        #hdulist.close()


#errorbar(t, lc, e, fmt = '.',  )
#plt.gca().invert_yaxis()

# ok, needed for reading into ysovar atlas: a file which has the entries:
# source_id mag   mag_err time
# 0         15.0  0.1     234...
# 0         14.9  0.11    234...
# etc. for every source id.


r_match = 0.1 * 1./3600.

bands = ['H', 'K', 'J']

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
#tab.add_row([0,0,0,0,0,0,0,0,0,0])
tab.add_row(np.zeros(10))

(ra, dec) = photometry_wcs.coords_from_ds9(input_info.masterregfile)
p_id = np.arange(0, len(ra))


for b in bands[:]:
    print b
    filepath =  input_info.resultfolder + '*YSO*/calibratedmags_' + b.lower() + '.dat' 
    filelist = glob.glob(filepath)
    filelist.sort()
    imagepath =  input_info.resultfolder + '*YSO*/' + b.lower() + '_long*_wcs.fits' 
    imagelist = glob.glob(imagepath)
    imagelist.sort()
    for i in p_id:
        print i
	for j in np.arange(0, len(filelist[:])):
	    #print filelist[j]
	    # get time stamp of observation.
	    # Ah, this would be faster if I had collected the HJDs already when I made the calibratedmags files. Bummer.
	    hdulist = pyfits.open(imagelist[j])
	    t = np.float(hdulist[0].header['HJULDATE'])
	    hdulist.close()
	    # find source in calibratedmags file of a given night.
	    a = asciitable.read(filelist[j])
	    distance = np.sqrt((ra[i] - a['RA'])**2 + (dec[i] - a['DEC'])**2)
	    if distance.min() < r_match:
		ind = np.where(distance == distance.min())[0][0]
		m = a[ind]['P' + b + '_cal']
		merr = a[ind]['P' + b + '_cal_err']
		tab.add_row(np.ones(len(names))*np.nan)
		tab[-1]['Pair_ID'] = i
		tab[-1]['ra'] = ra[i]
		tab[-1]['dec'] = dec[i]
		tab[-1]['t'] = t
		tab[-1][b] = m
		tab[-1][b + '_err'] = merr
	    
	


asciitable.write(tab, input_info.resultfolder + 'JHK_pairitel.csv' )



