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
import sys
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



# now transform found pairitel sources into wcs, find matching 2mass sources, and calibrate pairitel against 2mass. I use nominal errors of 0.1 mag for all pairitel sources instead of the output errors from the psf photometry. This is because sometimes there is large variability in a bright source (which has a small photometric error), which then would mess up the whole line fitting because that data point would have so much weight.

def check_band(filename):
    checkband = os.path.basename(filename)[0]
    (band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('', '', '', '', '')
    if checkband.lower() == 'h':
	(band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('2H', 'PH', 'PH_err', 'PH_cal', 'PH_cal_err')
    else:
	if checkband.lower() == 'k':
	    (band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('2K', 'PK', 'PK_err', 'PK_cal', 'PK_cal_err')
	else:
	    if checkband.lower() == 'j':
		(band2m, band2mp, band2mpe, band2mpc, band2mpce) = ('2J', 'PJ', 'PJ_err', 'PJ_cal', 'PJ_cal_err')
    
    return(checkband, band2m, band2mp, band2mpe, band2mpc, band2mpce)

def match_pairitel_2mass_g12(filename, data2mass, g12filename, matchradius, writetofile=False):
    imagefile = filename[0:-6]
    
    pairitel = asciitable.read(filename, Reader=asciitable.Daophot)
    (checkband, band2m, band2mp, band2mpe, band2mpc, band2mpce) = check_band(filename)
    
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='RA'))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='DEC'))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mp))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpe))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2m))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpc))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpce))
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='ysoflag'))
    
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
    # and check which 2mass-pairitel sources are not YSOs.
    g12 = asciitable.read(g12filename)
    g12 = g12._data.data
    g12 = Table(g12)
    
    
    for i in np.arange(0, len(pairitel)):
        # matching with 2mass:
	dist = dist_radec_fast(pairitel[i]['RA'], pairitel[i]['DEC'], data2mass['ra'], data2mass['dec'], scale = matchradius, unit ='deg')
	if dist.min() < matchradius:
	    ind = np.where(dist == dist.min())[0]
	    pairitel[i][band2m] = data2mass[ind][band2m]
	
	# matching with Guenther+2012 and flagging stuff that is a YSO:
	#ind_near = np.where(np.abs(g12['DEdeg'] - pairitel[i]['DEC']) <= r_match)[0]
	#if len(ind_near) > 0:
	distg12 = dist_radec_fast(pairitel[i]['RA'], pairitel[i]['DEC'], g12['RAdeg'], g12['DEdeg'], scale = matchradius, unit ='deg')
	ig12 = np.where(distg12 == distg12.min())[0]
	if ((distg12.min() < matchradius) & ((g12[ig12]['Class'] == 'star') | (g12[ig12]['Class'] == ''))):
	    pairitel[i]['ysoflag'] = 0.
    
    print np.where(pairitel[band2m]>-99999.)
    
    if writetofile == True:
        asciitable.write(pairitel, os.path.dirname(filename) + '/pairitel_matched_' + checkband + '.dat')
    
    return pairitel




def fit_to_2mass(filename, pairitel, chisq, chisq_const, chisq_real, chisq_real_const, writetofile=False):
    
    (checkband, band2m, band2mp, band2mpe, band2mpc, band2mpce) = check_band(filename)
    
    pairitel = asciitable.read(os.path.dirname(filename)+'/pairitel_matched_'+os.path.basename(filename)[0]+'.dat')
    good = np.where(pairitel[band2m]>-99999.)[0]
    good_const = np.where((pairitel[band2m]>-99999.) & (pairitel['ysoflag'] > -99999.))[0]
    print good
    print good_const
    
    f.clf()
    ax = f.add_subplot(111)
    ax.plot(pairitel[good][band2m], pairitel[good][band2mp], '.')
    ax.plot(pairitel[good_const][band2m], pairitel[good_const][band2mp], 'r.')
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
    pfinal_const = optimize.leastsq(errfunc, pinit, args=(pairitel[good_const][band2m], pairitel[good_const][band2mp], error[0:len(good_const)]), full_output=1)[0]
    
    x = np.linspace(np.min(pairitel[good][band2m]), np.max(pairitel[good][band2m]), 100)
    ax.plot(x, pfinal[0] + pfinal[1]*x, 'b-')
    ax.plot(x, pfinal_const[0] + pfinal_const[1]*x, 'r-')
    
    # calculate new pairitel magnitudes, using the 2mass fit to the non-YSO sources
    pairitel[band2mpc] = (pairitel[band2mp] - pfinal_const[0])/pfinal_const[1]
    pairitel[band2mpce] = (pairitel[band2mpe])/pfinal_const[1]
    
    #print pairitel[good][band2mpc]
    #print pairitel[good][band2mp]
    #print pairitel[good][band2m]
    
    # calculate several chi-squared values.
    # chisq = all 2mass matches, using constant errors
    n = np.where(filename == np.array(filelist))[0][0]
    chisq[n] = np.sum(((pairitel[good][band2mp] - (pfinal[0] + pfinal[1]*pairitel[good][band2m]))/(error))**2.)/len(pairitel[good])
    # chisq_const = 2mass matches which are non-YSOs, using constant errors 
    chisq_const[n] = np.sum(((pairitel[good_const][band2mp] - (pfinal[0] + pfinal[1]*pairitel[good_const][band2m]))/(error[0:len(good_const)]))**2.)/len(pairitel[good_const])
    # chisq_real =  all 2mass matches, using daophot errors
    chisq_real[n] = np.sum(((pairitel[good][band2mp] - (pfinal[0] + pfinal[1]*pairitel[good][band2m]))/(pairitel[good][band2mpe]))**2.)/len(pairitel[good])
    # chisq_real_const = 2mass matches which are non-YSOs, using daophot errors
    chisq_real_const[n] = np.sum(((pairitel[good_const][band2mp] - (pfinal[0] + pfinal[1]*pairitel[good_const][band2m]))/(pairitel[good_const][band2mpe]))**2.)/len(pairitel[good_const])
    
    # collect the difference of the calibrated Pairitel magnitudes to the fit derived from the non-YSO sources.
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name='deltafit'))
    for ind in good_const:
        pairitel[ind]['deltafit'] = (pairitel[ind][band2mpc] - pairitel[ind][band2m])
    
    ax.text(0.3, 0.7, 'red. chi**2 = ' + str(np.round(chisq_real[n],1)), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    ax.text(0.3, 0.8, 'red. chi**2, non-YSO = ' + str(np.round(chisq_real_const[n],1)), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    plt.show()
    
    # correct errors so that reduced chisquared = 1 for non-YSO sources.
    pairitel.add_column(Column(data = np.ones(len(pairitel['ID']))*-99999., name=band2mpce + '_corr'))
    if chisq_real_const[n] > 1.:
        pairitel[band2mpce + '_corr'] = pairitel[band2mpce] * np.sqrt(chisq_real_const[n])
    else:
        pairitel[band2mpce + '_corr'] = pairitel[band2mpce]
    
    # write results to file.
    if writetofile == True:
        asciitable.write(pairitel, os.path.dirname(filename) + '/calibratedmags_' + checkband + '.dat')
    
    return (chisq, chisq_const, chisq_real, chisq_real_const, pfinal, pfinal_const, pairitel)

def get_all_caldata():
    filepath =  input_info.resultfolder + '*YSO*/calibratedmags_j.dat'
    filelist = glob.glob(filepath)
    filelist.sort()
    
    mydict = dict({})
    
    for b in ['j', 'h', 'k']:
        filepath =  input_info.resultfolder + '*YSO*/calibratedmags_'+ b +'.dat'
        filelist = glob.glob(filepath)
        filelist.sort()
            
        data = asciitable.read(filelist[0])
	for filename in filelist[:]:
	    dat_tmp = asciitable.read(filename)
	    for i in np.arange(0, len(dat_tmp)):
		data.add_row(dat_tmp[i])
        
        mydict['data' + b.upper()] = data
    
    return mydict


def check_calmags(data, band):
    (checkband, band2m, band2mp, band2mpe, band2mpc, band2mpce) = check_band(band)
    good = np.where(data['deltafit']>-99999.)[0]
    #plt.figure()
    #plt.plot(data[good][band2mpc], data[good]['deltafit'], 'bo')
    #plt.plot(data[good][band2mpc], data[good]['deltafit']/data[good][band2mpce]/50., 'ro')
    
    #good_bright = np.where((data['deltafit']>-99999.) & (data[band2mpc] < 13.))[0]
    #good_dim = np.where((data['deltafit']>-99999.) & (data[band2mpc] >= 13.))[0]
    
    #plt.figure()
    #plt.hist(data[good_bright]['deltafit'], bins=50, alpha=0.5)
    #plt.hist(data[good_dim]['deltafit'], bins=50, alpha=0.5)
    
    #plt.figure()
    #plt.hist(data[good_bright]['deltafit']/data[good_bright][band2mpce], bins=50, alpha=0.5)
    #plt.hist(data[good_dim]['deltafit']/data[good_dim][band2mpce], bins=50, alpha=0.5)
    
    plt.clf()
    errs=data[good]['deltafit']/data[good][band2mpce + '_corr']
    plt.hist(errs, bins=50, alpha=0.5)
    errs.sort()
    i_low = int(0.16*len(errs))
    i_high = int(0.84*len(errs))
    #print errs[i_low], errs[i_high]
    
    errs_new=data[good]['deltafit']/data[good][band2mpce]/np.mean(np.array([np.abs(errs[i_low]), np.abs(errs[i_high])]))
    #plt.figure()
    #plt.hist(errs_new, bins=50, alpha=0.5)
    errs_new.sort()
    i_low = int(0.16*len(errs_new))
    i_high = int(0.84*len(errs_new))
    #print errs_new[i_low], errs_new[i_high]
    
    return np.mean(np.array([np.abs(errs[i_low]), np.abs(errs[i_high])]))



# load 2mass catalogue from disk:
data2mass = asciitable.read(input_info.resultfolder + '2mass_for_cal.dat')

# find all .als.1 files:
filepath =  input_info.resultfolder + '*YSO*/*.als.1' 
filelist = glob.glob(filepath)
filelist.sort()

# find matching sources in 2mass and guenther+2012 for each source in each night.
for filename in filelist[:]:
    pairitel = match_pairitel_2mass_g12(filename, data2mass, g12filename='/swiper.real/kpoppen/IR/YSOVAR/iras20050/Guenther2012_table.txt', matchradius=1./3600., writetofile=True)


# initialize output arrays
chisq = np.ones(len(filelist))*-99999.
chisq_const = np.ones(len(filelist))*-99999.
chisq_real = np.ones(len(filelist))*-99999.
chisq_real_const = np.ones(len(filelist))*-99999.

# initialize plot window
f = plt.figure()

# fit pairitel data to 2mass data, using non-YSO sources as identified by guenther+2012
for filename in filelist[:]:
    (chisq, chisq_const, chisq_real, chisq_real_const, pfinal, pfinal_const, pairitel) = fit_to_2mass(filename, pairitel, chisq, chisq_const, chisq_real, chisq_real_const, writetofile=True)

# collect all data from all night in a dictionary
datadict = get_all_caldata()



#corrfactor = np.zeros(3)
#for b in np.array(['J', 'H', 'K']):
    #corrfactor[np.where(b == np.array(['J', 'H', 'K']))[0]] = check_calmags(datadict['data'+b][:], b)

#print datadict['dataJ']['PJ_cal_err_corr']/datadict['dataJ']['PJ_cal_err']

#ratio = datadict['dataJ']['PJ_cal_err_corr']/datadict['dataJ']['PJ_cal_err']
#good = np.where(np.isnan(ratio) == False) 
#print ratio[good].min()
#print ratio[good].max()

#plt.hist(ratio[good], bins=100)




