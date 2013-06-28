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

# The pairitel_pyraf_0.py script makes a copy of all the raw data, so in case something gets scrambled, you will still have an unaltered copy of the data.

# this is your parent folder in which you stored all your unprocessed data. Add '/' at the end.
rawdatafolder='/swiper.real/kpoppen/IR/GD1215/'

# this is your parent folder in which you want all your processes files to be stored. Add '/' at the end.
resultfolder='/swiper.real/kpoppen/IR/GD1215/output/'

# this is the folder in which you stored the pairitel reduction scripts.
pairitel_scripts_path = '/data/swolk/kpoppen/Dropbox/MyPython/YSOVAR/PAIRITEL'

# these are parameters for the initial cropping of the nightly images. 
# threshold: the images are co-added per night, and the nightly images have lower exposure at the outer edges. You will throw away parts at the edges which have exposure time exp_edge/exp_center < threshold.
threshold = 0.8
# You will furthermore require that a rectangular region with good exposure exists, and the minimal dimensions of that rectangle are in pixels:
min_width = 100
min_height = 100

# these are the approximate central coordinates of your cluster (in degrees).
# this is necessary because for most clusters some of the Pairitel observations were pointed somewhere else, and we don't need that data for now.
ra_cluster = 92.706
dec_cluster = -6.19
# this is the approximate width of the Pairitel images in degrees. Observations whose central point is farther away than width_images from the central point of your master image will not be analyzed (because they have no overlap with the master image).
width_images = 0.1

# is wcs correcting by hand necessary? The put this value to "yes"
badwcs = "yes"
# where's the file for that manual correction?
wcs_corrbyhand_file = resultfolder + "wcs_by_hand.txt"

# this sets the H band magnitude at which psf fitting stars become too bright (i.e. saturated) for sensible psf fitting. Set this too a small value (15 or so) if selecting psf stars by hand anyway.
satmag = 15.

# set some parameters for the reduction. These defaults should be fine.
fwhm_psf = 4.5
radius_annulus = 4.5
width_annulus = 4.0
fwhm_apertures = 4.5
psf_function = "gauss"
psf_varorder = 0 # gauss + pixeltable
psf_fitrad = 2.2
psfrad = 6.0
fitsky = "yes"
skyannulus = 4.5
width_skyannulus = 4.
groupsky = "yes"
skyvalue = 0.

# define (non-empty) suffix for the magnitude file which will be created by the preliminary aperture photometry.
photfilesuffix = '.magapphot'

# define radius for which faint stars are removed from the wings of the psf fitting stars. 10. or 15. is an ok value.
pfscleaningradius = 15.

# this is the name of (including the path to) your master image.
# you will only know this after running the first two pairitel_pyraf scripts.
masterimage = resultfolder + 'YSO.47.23/h_long_YSO.47.23_wcs.fits'

# this is the path+name of the ds9 region file which contains the stars you selected for psf fitting.
# you will only know this after running the first two pairitel_pyraf scripts and defining the psf fitting stars in ds9 by hand (see PAIRITEL_RAEDME.txt). Do not change the following line.
psfstarfile = resultfolder + 'psfstars.reg'

# if you encounter bad files after running the aperture photometry, list them here as shown and put "bad_files_exist" to "yes". (This could be files in which the good exposure part is so tiny that there are almost no stars in it and no psf can be fitted.)
bad_files_exist = "yes"
bad_exposures = [resultfolder + 'YSO.47.14/k_long_YSO.47.14_wcs.fits', resultfolder + 'YSO.47.23/k_long_YSO.47.23_wcs.fits']

# after running the psf photometry, save aa DS9 region file with the detected sources in the masterimage (see PAIRITEL_RAEDME.txt). Do not change the following line.
masterregfile = resultfolder + 'masterstars.reg'

# provide a catalogue file name here that is readable by astropy.io.ascii (for example a standard electronic data table from a paper) which lists the class 1 and class 2 objects in your cluster. This can be from Rob's paper, for example. 
# if no catclogue exists, put catalogfile = ''
catalogfile = ''
# define what the columns containing RA, DEC, and YSO-class are named in that file:
rakeyword = 'col2'
deckeyword = 'col3'
classkeyword = 'col7'
# and define what kind of idenitifiers the non-YSO sources have in the YSO-class column:
list_of_nonyso_values = ['-1', '0']

# provide a threshold for the preliminarily calibrated light curves. Sources with a 
threshold_lc = 0.15

