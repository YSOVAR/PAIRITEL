# -*- coding: utf-8 -*-
#write this as script for now, turn parts in sensible functions later 
#keep all intermediate steps so I can track down errors
#change this behavious later e.g. use os.tmpfile or maybe tempfile.XXX()

# start /usr/local/bin/ipython --colors lightbg # use this python because it knows where pyraf is.

#import matplotlib.pyplot as plt
import pyfits
import shutil
from pyraf import iraf
import glob
import os
import sys
import asciitable
sys.path.append("/data/hguenther/Dropbox/code/python")
import atpyextensions
sys.path.append("/data/hguenther/Dropbox/code/python/atpy")
import atpy
sys.path.append("/swiper.real/kpoppen/Dropbox/MyPython/PAIRITEL/")
from atpyextensions import catalog
import numpy as np
import string
import pickle
from copy import deepcopy

import input_info
reload(input_info)
import photometry
reload(photometry)
import photometry_both
reload(photometry_both)

iraf.noao()
iraf.digiphot()
iraf.daophot()


photometry.set_Pairitel_params()

# get list of all long exposure files, both actual target observations and the weight files.
datapath = input_info.rawdatafolder + '*YSO*/*_long_*coadd*fits'
imlist = glob.glob(datapath)
imlist.sort()

# this finds the observations which are not actually centered on the cluster.
badlist = photometry_both.make_badlist(imlist)

# sanitycheck excludes observations which do not include the cluster from the analysis file list.
imlist = photometry_both.sanitycheck(imlist, badlist)

# this copies all data which will be used in the analysis to the resultfolder (so that you have a copy of the unaltered data in the original folder), and then, in the resultfolder, trims the images and normalizes them by the exposure time mask.
photometry.prepare_files(imlist, input_info.threshold, input_info.min_width, input_info.min_height, input_info.resultfolder)

# get a list of the exposure-normalized and trimmed images.
datapath_sky = input_info.resultfolder + '*YSO*/*_coadd_normed.fits' 
datalist = glob.glob(datapath_sky)
datalist.sort()

# correct coordinate shifts by comparing each image to the 2MASS catalogue, apply necessary coordinate corrections.
photometry.correct_coordinates(datalist, radius=10.)

#os.getcwd()

