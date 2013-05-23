# -*- coding: utf-8 -*-

# the scripts assume that your data is organized like this:
# resultfolder>
#                one folder per night>
#                (usually named YSO.xx.xx)
#                                       containing at least these files per night:
#                                       h_long_xxxxxxx_coadd.fits
#                                       h_long_xxxxxxx_coadd.weight.fits
#                                       j_long_xxxxxxx_coadd.fits
#                                       j_long_xxxxxxx_coadd.weight.fits
#                                       k_long_xxxxxxx_coadd.fits
#                                       k_long_xxxxxxx_coadd.weight.fits

# it is generally a good idea to make a copy of your entire data before you start the analysis, in case something goes wrong and an original data file gets scrambled.

# this is your parent folder in which you stored all your unprocessed data.
rawdatafolder='/swiper.real/kpoppen/IR/I20050/'


# this is your parent folder in which you want all your processes files to be stored.
resultfolder='/swiper.real/kpoppen/IR/I20050/output/'

#paramdir='/swiper.real/kpoppen/IR/I20050'

# this is the folder in which you stored the pairitel scripts.
pairitel_scripts_path = '/data/swolk/kpoppen/Dropbox/MyPython/PAIRITEL'

# these are parameters for the initial cropping of the nightly images. 
# threshold: the images are co-added per night, and the nightly images have lower exposure at the outer edges. You will throw away parts at the edges which have exposure time exp_edge/exp_center < threshold.
threshold = 0.9
# You will furthermore require that a rectangular region with good exposure exists, and the minimal dimensions of that rectangle are in pixels:
min_width = 100
min_height = 100

# these are the approximate central coordinates of your cluster (in degrees).
# this is necessary because, for my cluster, some of the Pairitel observations were pointed somewhere else, and we don't need that data for now.
ra_cluster = 301.75
dec_cluster = 27.5
# this is the approximate width of the Pairitel images in degrees. Observations whose central point is farther away than width_images from the central point of your master image will not be analyzed (because they have no overlap with the master image).
width_images = 0.1

# this sets the H band magnitude at which psf fitting stars become too bright (i.e. saturated) for sensible psf fitting. Set this too a small value (10 or so) if selecting psf stars by hand anyway.
satmag = 10.

# this is the name of (including the path to) your master image.
# you will only know this after running the first pairitel_pyraf script.
masterimage = resultfolder + 'YSO.38.12/h_long_YSO.38.12_wcs.fits'

# this is the path+name of the ds9 region file which contains the stars you selected for psf fitting.
# you will only know this after running the first pairitel_pyraf script and defining the psf fitting stars in ds9 by hand (see PAIRITEL_RAEDME.txt).
psfstarfile = resultfolder + 'psfstars.reg'

masterregfile = resultfolder + 'masterstars.reg'


#mastermagfile = '/swiper.real/kpoppen/IR/I20050/output/YSO.38.15/h_long_YSO.38.15_wcs.fits.mag.1'
#masteralsfile = '/swiper.real/kpoppen/IR/I20050/output/YSO.38.15/h_long_YSO.38.15_wcs.fits.als.3'
#mastercoo = '/swiper.real/kpoppen/IR/I20050/output/YSO.38.15/h_long_YSO.38.15_wcs.mastercoo.1'
#masterimage = '/swiper.real/kpoppen/IR/I20050/output/YSO.38.15/h_long_YSO.38.15_wcs.fits'
