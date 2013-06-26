# -*- coding: utf-8 -*-

# these are functions used by script_pairitel_wcs.py.

import numpy as np
import matplotlib.pyplot as plt
import aplpy
import astropy.io.fits as pyfits
import astropy
import astropy.io.ascii as asciitable
from copy import deepcopy
import glob
from scipy import optimize
from astropy.table import Table, Column
from astropy.wcs import WCS
from astropy.io import fits



# this is a dummy iraf header for the output .pst files.
headerpst = [
'#K IRAF       = NOAO/IRAFV2.14.1        version    %-23s',
'#K USER       = kpoppen                 name	   %-23s',
'#K HOST       = swiper                  computer   %-23s',
'#K DATE       = 0000-00-00              yyyy-mm-dd %-23s',
'#K TIME       = 00:00:00                hh:mm:ss   %-23s',
'#K PACKAGE    = apphot                  name	   %-23s',
'#K TASK       = phot                    name	   %-23s',
'#',
'#K SCALE      = 1.                      units	   %-23.7g',
'#K FWHMPSF    = 4.5                     scaleunit  %-23.7g',
'#K EMISSION   = yes                     switch     %-23b',
'#K DATAMIN    = -7.9655                 counts     %-23.7g',
'#K DATAMAX    = 15000.                  counts     %-23.7g',
'#K EXPOSURE   = EXPTIME                 keyword    %-23s',
'#K AIRMASS    = AIRMASS                 keyword    %-23s',
'#K FILTER     = FILTER                  keyword    %-23s',
'#K OBSTIME    = DATE-OBS                keyword    %-23s',
'#',
'#K NOISE      = poisson                 model	   %-23s',
'#K SIGMA      = 1.534                   counts     %-23.7g',
'#K GAIN       = GAIN                    keyword    %-23s',
'#K EPADU      = 80.21409                e-/adu     %-23.7g',
'#K CCDREAD    = 10.0                    keyword    %-23s',
'#K READNOISE  = 0.                      e-         %-23.7g',
'#',
'#K CALGORITHM = none                    algorithm  %-23s',
'#K CBOXWIDTH  = 9.                      scaleunit  %-23.7g',
'#K CTHRESHOLD = 0.                      sigma	   %-23.7g',
'#K MINSNRATIO = 1.                      number     %-23.7g',
'#K CMAXITER   = 10                      number     %-23d',
'#K MAXSHIFT   = 1.                      scaleunit  %-23.7g',
'#K CLEAN      = no                      switch     %-23b',
'#K RCLEAN     = 1.                      scaleunit  %-23.7g',
'#K RCLIP      = 2.                      scaleunit  %-23.7g',
'#K KCLEAN     = 3.                      sigma	   %-23.7g',
'#',
'#K SALGORITHM = centroid                algorithm  %-23s',
'#K ANNULUS    = 100.                    scaleunit  %-23.7g',
'#K DANNULUS   = 10.                     scaleunit  %-23.7g',
'#K SKYVALUE   = 0.                      counts     %-23.7g',
'#K KHIST      = 3.                      sigma	   %-23.7g',
'#K BINSIZE    = 0.1                     sigma	   %-23.7g',
'#K SMOOTH     = no                      switch     %-23b',
'#K SMAXITER   = 10                      number     %-23d',
'#K SLOCLIP    = 0.                      percent    %-23.7g',
'#K SHICLIP    = 0.                      percent    %-23.7g',
'#K SNREJECT   = 50                      number     %-23d',
'#K SLOREJECT  = 3.                      sigma	   %-23.7g',
'#K SHIREJECT  = 3.                      sigma	   %-23.7g',
'#K RGROW      = 0.                      scaleunit  %-23.7g',
'#',
'#K WEIGHTING  = constant                model	   %-23s',
'#K APERTURES  = 4.0                     scaleunit  %-23s',
'#K ZMAG       = 25.                     zeropoint  %-23.7g',
'#',
'#K IMAGE      = xxxx                    imagename  %-23s',
'#K MAXNPSF    = 10                      number     %-23d',
'#K NEWSCALE   = 1.                      units	   %-23.7g',
'#K PSFRAD     = 11.                     scaleunit  %-23.7g',
'#K FITRAD     = 3.                      scaleunit  %-23.7g',
'#',
'#N ID    XCENTER   YCENTER   MAG         MSKY',
'#U ##    pixels    pixels    magnitudes  counts',
'#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g',
'#'
]

headercoo = [
'#K IRAF       = NOAO/IRAFV2.14.1        version    %-23s',
'#K USER       = kpoppen                 name	   %-23s',
'#K HOST       = swiper                  computer   %-23s',
'#K DATE       = 0000-00-00              yyyy-mm-dd %-23s',
'#K TIME       = 00:00:00                hh:mm:ss   %-23s',
'#K PACKAGE    = apphot                  name	   %-23s',
'#K TASK       = daofind                 name	   %-23s',
'#',
'#K SCALE      = 1.                      units	   %-23.7g',
'#K FWHMPSF    = 4.5                     scaleunit  %-23.7g',
'#K EMISSION   = yes                     switch     %-23b',
'#K DATAMIN    = -19.4022                counts     %-23.7g',
'#K DATAMAX    = 15000.                  counts     %-23.7g',
'#K EXPOSURE   = EXPTIME                 keyword    %-23s',
'#K AIRMASS    = AIRMASS                 keyword    %-23s',
'#K FILTER     = FILTER                  keyword    %-23s',
'#K OBSTIME    = DATE-OBS                keyword    %-23s',
'#',
'#K NOISE      = poisson                 model	   %-23s',
'#K SIGMA      = 3.857                   counts     %-23.7g',
'#K GAIN       = GAIN                    keyword    %-23s',
'#K EPADU      = 70.18733                e-/adu     %-23.7g',
'#K CCDREAD    = 10.0                    keyword    %-23s',
'#K READNOISE  = 0.                      e-         %-23.7g',
'#',
'#K IMAGE      = xxxx                    imagename  %-23s',
'#K FWHMPSF    = 4.5                     scaleunit  %-23.7g',
'#K THRESHOLD  = 4.                      sigma	   %-23.7g',
'#K NSIGMA     = 1.5                     sigma	   %-23.7g',
'#K RATIO      = 1.                      number     %-23.7g',
'#K THETA      = 0.                      degrees    %-23.7g',
'#',
'#K SHARPLO    = 0.2                     number     %-23.7g',
'#K SHARPHI    = 1.                      number     %-23.7g',
'#K ROUNDLO    = -1.                     number     %-23.7g',
'#K ROUNDHI    = 1.                      number     %-23.7g',
'#',
'#N XCENTER   YCENTER   MAG      SHARPNESS   SROUND      GROUND      ID         \ ',
'#U pixels    pixels    #        #           #           #           #          \ ',
'#F %-13.3f   %-10.3f   %-9.3f   %-12.3f     %-12.3f     %-12.3f     %-6d       \ ',
'#'
]


def read_master_sourcefile(filename):
    data = asciitable.read(filename, Reader=asciitable.Daophot)


def read_user_psffile(psfstarfile):
    # reads a list of psf fitting stars provided by the user as a DS9 region file, saved in decimal wcs units.
    f = open(psfstarfile, 'r')
    psfstars = f.read()
    f.close()
    psfstars = psfstars.split('\n')
    
    ra_psf = np.array([])
    dec_psf = np.array([])
    
    for i in np.arange(0, len(psfstars)):
	if psfstars[i][0:6] == 'circle':
	    line = psfstars[i].split('(')
	    line = line[1].split(',')
	    ra_psf = np.append(ra_psf, float(line[0]))
	    dec_psf = np.append(dec_psf, float(line[1]))
    
    return (ra_psf, dec_psf)

def make_pst_file_from_ds9reg_old(psfstarfile, magfile, imagefile, outfile):
    # this writes the data from a user-supplied ds9 region file with stars for PSF fitting in such a format that it looks like an iraf .pst file; this includes the transformation to image coordinates.
    (ra_psf, dec_psf) = read_user_psffile(psfstarfile)
    gc = aplpy.FITSFigure(imagefile)
    print imagefile
    x_psf = deepcopy(ra_psf)
    y_psf = deepcopy(dec_psf)
    for i in np.arange(0,len(ra_psf)):
        (x_psf[i], y_psf[i]) = gc.world2pixel(ra_psf[i], dec_psf[i])
    
    # test if any of those psf stars are outside the actual image (i.e. have x or y coordinates <0):
    if (x_psf<0).any():
        bad = np.where(x_psf < 0)[0]
        x_psf = np.delete(x_psf, bad)
        y_psf = np.delete(y_psf, bad)
    if (y_psf<0).any():
        bad = np.where(y_psf < 0)[0]
        x_psf = np.delete(x_psf, bad)
        y_psf = np.delete(y_psf, bad)
    
    id_psf = np.zeros(len(x_psf), int)
    mag_psf = np.zeros(len(x_psf))
    msky_psf = np.zeros(len(x_psf))
    apdat = asciitable.read(magfile, Reader=asciitable.Daophot)
    apdat = apdat._data.data
    # match the translated coordinates from the ds9 reg file to magnitudes and other values in the .mag file from the master image. This should be easy, because all reasonable PSF fitting stars should have been detected in the aperture photometry.
    for i in np.arange(0, len(x_psf)):
	distance = np.sqrt((apdat['XCENTER']-x_psf[i])**2 + (apdat['YCENTER']-y_psf[i])**2)
	match = np.where(distance == np.min(distance))[0]
	if distance[match] < 1:
	    print 'matched successfully.'
	    id_psf[i] = apdat[match]['ID']
	    mag_psf[i] = apdat[match]['MAG']
	    msky_psf[i] = apdat[match]['MSKY']
    
    allstuff = zip(mag_psf, id_psf, x_psf, y_psf, msky_psf) # sort by first argument in zip.
    allstuff.sort()
    (mag_psf, id_psf, x_psf, y_psf, msky_psf) = zip(*allstuff)
    
    # construct output lines which look like in a .pst file...
    line = []
    for i in np.arange(0, len(x_psf)):
	line.append(str(id_psf[i]).ljust(9) + str(x_psf[i])[0:7].ljust(10) + str(y_psf[i])[0:7].ljust(10) + str(mag_psf[i])[0:6].ljust(12) + str(msky_psf[i])[0:10].ljust(10))
    
    # ...and write the stuff into the file which will be used for fitting the psf stars.
    f = open(outfile, 'w')
    for i in np.arange(0, len(headerpst)):
        f.write(headerpst[i] + '\n')
    
    for i in np.arange(0, len(line)):
        f.write(line[i] + '\n')
    
    f.close()
    



def make_pst_file_from_ds9reg(psfstarfile, magfile, imagefile, outfile):
    # this writes the data from a user-supplied ds9 region file with stars for PSF fitting in such a format that it looks like an iraf .pst file; this includes the transformation to image coordinates.
    (ra_psf, dec_psf) = read_user_psffile(psfstarfile)
    hdus = fits.open(imagefile)
    wcs = WCS(hdus[0].header)
    print imagefile
    x_psf = deepcopy(ra_psf)
    y_psf = deepcopy(dec_psf)
    x_psf, y_psf = wcs.wcs_world2pix(ra_psf, dec_psf, 1)
    hdus.close()
    #for i in np.arange(0,len(ra_psf)):
        #(x_psf[i], y_psf[i]) = gc.world2pixel(ra_psf[i], dec_psf[i])
    
    # test if any of those psf stars are outside the actual image (i.e. have x or y coordinates <0):
    if (x_psf<0).any():
        bad = np.where(x_psf < 0)[0]
        x_psf = np.delete(x_psf, bad)
        y_psf = np.delete(y_psf, bad)
    if (y_psf<0).any():
        bad = np.where(y_psf < 0)[0]
        x_psf = np.delete(x_psf, bad)
        y_psf = np.delete(y_psf, bad)
    
    id_psf = np.zeros(len(x_psf), int)
    mag_psf = np.zeros(len(x_psf))
    msky_psf = np.zeros(len(x_psf))
    apdat = asciitable.read(magfile, Reader=asciitable.Daophot)
    apdat = apdat._data.data
    # match the translated coordinates from the ds9 reg file to magnitudes and other values in the .mag file from the master image. This should be easy, because all reasonable PSF fitting stars should have been detected in the aperture photometry.
    for i in np.arange(0, len(x_psf)):
	distance = np.sqrt((apdat['XCENTER']-x_psf[i])**2 + (apdat['YCENTER']-y_psf[i])**2)
	match = np.where(distance == np.min(distance))[0]
	if distance[match] < 1:
	    print 'matched successfully.'
	    id_psf[i] = apdat[match]['ID']
	    mag_psf[i] = apdat[match]['MAG']
	    msky_psf[i] = apdat[match]['MSKY']
    
    allstuff = zip(mag_psf, id_psf, x_psf, y_psf, msky_psf) # sort by first argument in zip.
    allstuff.sort()
    (mag_psf, id_psf, x_psf, y_psf, msky_psf) = zip(*allstuff)
    
    # construct output lines which look like in a .pst file...
    line = []
    for i in np.arange(0, len(x_psf)):
	line.append(str(id_psf[i]).ljust(9) + str(x_psf[i])[0:7].ljust(10) + str(y_psf[i])[0:7].ljust(10) + str(mag_psf[i])[0:6].ljust(12) + str(msky_psf[i])[0:10].ljust(10))
    
    # ...and write the stuff into the file which will be used for fitting the psf stars.
    f = open(outfile, 'w')
    for i in np.arange(0, len(headerpst)):
        f.write(headerpst[i] + '\n')
    
    for i in np.arange(0, len(line)):
        f.write(line[i] + '\n')
    
    f.close()
    







def dummycheck_psf():
    # this is just for checking that the psf really has there ringlike deviations from a gaussian. Spoiler: yes, they're real.
    image = '/swiper.real/kpoppen/IR/I20050/output/YSO.38.11/h_long_YSO.38.11_wcs.fits'
    gc = aplpy.FITSFigure(image)
    gc.show_colorscale(vmin=-10, vmax=100)
    
    hdulist = pyfits.open(image)
    data = hdulist[0].data
    hdulist.close
    
    plt.imshow(data)
    
    a1=270
    a2=280
    a3=110
    a4=121
    dat_crop = data[a1:a2]
    cutout = np.zeros([a2-a1,a4-a3])
    
    for i in np.arange(0,a2-a1):
	cutout[i][:] = dat_crop[i][a3:a4]
    
    flat = cutout.flatten()
    flat.sort()
    baselevel = median(flat[0:12])

    cutout = cutout - baselevel
    plt.figure()
    plt.imshow(cutout)
    plt.colorbar()
    p = fitgaussian(cutout)
    fit = gaussian(p[0], p[1], p[2], p[3], p[4])
    i_array = deepcopy(cutout)
    i_array[:][:] = 0
    for i in np.arange(0, cutout.shape[0]):
	for j in np.arange(0, cutout.shape[1]):
	    i_array[i][j] = fit(i,j)
    
    plt.figure()
    plt.imshow(i_array)
    plt.colorbar()
    resi = cutout - i_array
    plt.figure()
    plt.imshow(resi)
    plt.colorbar()

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


#def make_mastercoo_for_masterimage(filename, outname, ds9=False):
    ## read the final .als file for the masterimage.
    #data = asciitable.read(filename, Reader=asciitable.Daophot)
    ## needed for .coo: XCENTER, YCENTER, MAG, SHARPNESS, SROUND, GROUND, ID
    ## SROUND and GROUND are not in the .als file; just write zeroes.
    ## subtract 26 from the magnitudes, because .coo files have wrong values. I could probably also put zeroes there.
    #line = []
    #for i in np.arange(0, len(data)):
	#line.append('   ' + str(data[i]['XCENTER']).ljust(10) + str(data[i]['YCENTER']).ljust(10) + str(data[i]['MAG']-26).ljust(9) + str(data[i]['SHARPNESS']).ljust(12) + '0.000'.ljust(12) + '0.000'.ljust(12) + str(data[i]['ID']).ljust(6))
    
    ## now write the mastercoo file.
    #f = open(outname, 'w')
    #for i in np.arange(0, len(headercoo)):
	#f.write(headercoo[i] + '\n')
    
    #for i in np.arange(0, len(line)):
	#f.write(line[i] + '\n')
    
    #f.close()
    
    ## for testing purposes: also make ds9 region file for the newly created coo list. Spoiler: yes, it works. Whooeee!
    #if ds9 == True:
	#f = open(outname + '.reg', 'w')
	#f.write('# Region file format: DS9 version 4.1' + '\n')
	#f.write('# Filename bla' + '\n')
	#f.write('global color=blue dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1' + '\n')
	#f.write('image' + '\n')
	#for i in np.arange(0,len(data)):
	    #f.write('circle(' + str(np.round(data[i]['XCENTER'],3)) + ',' + str(np.round(data[i]['YCENTER'],3)) + ',1.0")' + '\n')
        
	#f.close()



def make_mastercoos_for_images(masterimage, masterreg, imagefile, outfile, fantasymag, ds9=False):
    # this reads the masterregfile derived from the psf photometry of the masterimage, and translates in into image coordinates for each image. It then makes a .coo-like file for each image. This will later be used for the psf photometry of all images.
    
    #data = asciitable.read(mastercoo, Reader=asciitable.Daophot)
    f = open(masterreg, 'r')
    data = f.read()
    f.close()
    data = data.split('\n')
    data = data[3:-1] # first three lines are header, last line is empty.
    ra = np.zeros(len(data))
    dec = np.zeros(len(data))
    for i in np.arange(0, len(data)):
        ra[i] = np.float(data[i].split('(')[1].split(',')[0])
        dec[i] = np.float(data[i].split('(')[1].split(',')[1])
        
    #gc = aplpy.FITSFigure(masterimage)
    x = deepcopy(ra)
    y = deepcopy(dec)
    hdus = fits.open(imagefile)
    wcs = WCS(hdus[0].header)
    print imagefile
    x, y = wcs.wcs_world2pix(ra, dec, 1)
    hdus.close()
    
    ##ra = deepcopy(data['XCENTER'])
    ##dec = deepcopy(data['YCENTER'])
    ### first translate xcenter/ycenter coordinates of mastercoo to wcs:
    ##for i in np.arange(0,len(x)):
	##(ra[i], dec[i]) = gc.pixel2world(data[i]['XCENTER'], data[i]['YCENTER'])
    
    ##plt.close()
    #gc = aplpy.FITSFigure(imagefile)
    ## then translate ra/dec into x/y coordinates in the imagefile.
    #for i in np.arange(0,len(x)):
	#(x[i], y[i]) = gc.world2pixel(ra[i], dec[i])
    
    #plt.close()
    # now put it into a coo-like line structure:
    line = []
    for i in np.arange(0, len(data)):
	line.append('   ' + str(round(x[i],3)).ljust(10) + str(round(y[i],3)).ljust(10) + str(fantasymag).ljust(9) + '0.000'.ljust(12) + '0.000'.ljust(12) + '0.000'.ljust(12) + str(i).ljust(6))
    
    f = open(outfile, 'w')
    for i in np.arange(0, len(headercoo)):
	f.write(headercoo[i] + '\n')

    for i in np.arange(0, len(line)):
	f.write(line[i] + '\n')
    
    f.close()
    
    # for testing purposes: also make ds9 region file for the newly created coo list. Spoiler: yes, it works. Whooeee!
    if ds9 == True:
	f = open(outfile + '.reg', 'w')
	f.write('# Region file format: DS9 version 4.1' + '\n')
	f.write('# Filename bla' + '\n')
	f.write('global color=blue dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1' + '\n')
	f.write('image' + '\n')
	for i in np.arange(0,len(x)):
	    f.write('circle(' + str(np.round(x[i],3)) + ',' + str(np.round(y[i],3)) + ',1.0")' + '\n')
        
	f.close()

def coords_from_ds9(filename):
    f = open(filename, 'r')
    data = f.read()
    f.close()
    data = data.split('\n')
    data = data[3:-1] # first 3 lines header, last line empty
    
    ra = np.zeros(len(data))
    dec = np.zeros(len(data))
    for i in np.arange(0, len(data)):
        ra[i] = np.float(data[i].split('(')[1].split(',')[0])    
        dec[i] = np.float(data[i].split('(')[1].split(',')[1])
    
    return (ra, dec)

