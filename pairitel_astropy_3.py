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

# get absolute magnitudes by fitting to 2mass magnitudes.
# get a 2mass catalogue for the cluster:


ra = input_info.ra_cluster
dec = input_info.dec_cluster
radius = 30.
minmag = 15.

data=urllib.urlopen('http://archive.tuc.noao.edu/cgi-bin/scat?catalog=tmc&ra='+str(ra)+'&dec='+str(dec)+'&sys=J2000&mrad='+str(radius)+'&nstar=-1&mag='+str(minmag)).read()

data = data.split('\n')
data = data[10:len(data)-1]

ind = np.array([], int)
for i in np.arange(0, len(data)):
    if ('L' in data[i]):
        ind = np.append(ind, i)

if len(ind) > 0:
    for i in ind:
        data[i] = '!'


for i in np.arange(0, len(ind)):
    data.remove('!')


id2m = np.zeros([len(data)])
ra2m = np.zeros([len(data)])
dec2m = np.zeros([len(data)])
j2m = np.zeros([len(data)])
h2m = np.zeros([len(data)])
k2m = np.zeros([len(data)])



for i in np.arange(0, len(data)):
    line = data[i].split('\t')
    id2m[i] = line[0]
    ra2m[i] = float(line[1].split(':')[0])*15 + float(line[1].split(':')[1])*15/60. + float(line[1].split(':')[2])*15/3600.
    dec2m[i] = float(line[2].split(':')[0]) + float(line[2].split(':')[1])/60. + float(line[2].split(':')[2])/3600.
    j2m[i] = float(line[3])
    h2m[i] = float(line[4])
    k2m[i] = float(line[5])

names = ('2massid', 'ra', 'dec', '2J', '2H', '2K')

data2mass = Table([id2m, ra2m, dec2m, j2m, h2m, k2m], names=names)

#data2mass.add_column(Column(data = np.ones(len(data2mass['ra']))*-99999., name='PJ'))
#data2mass.add_column(Column(data = np.ones(len(data2mass['ra']))*-99999., name='PH'))
#data2mass.add_column(Column(data = np.ones(len(data2mass['ra']))*-99999., name='PK'))
#data2mass.add_column(Column(data = np.ones(len(data2mass['ra']))*-99999., name='PJ_cal'))
#data2mass.add_column(Column(data = np.ones(len(data2mass['ra']))*-99999., name='PH_cal'))
#data2mass.add_column(Column(data = np.ones(len(data2mass['ra']))*-99999., name='PK_cal'))

asciitable.write(data2mass, input_info.resultfolder + '2mass_for_cal.dat')
