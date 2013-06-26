Recipe for Pairitel data reduction
----------------------------------

At CfA, I have to use two different python shells because Astropy and Pyraf cannot be loaded in the same python shell here. If you can load both pyraf and astropy in the same shell on your system, just run all the scripts from the same shell.

Step 1

DO THIS:
- Copy the input_info_sample.py file to input_info.py.
- Modify in input_info.py:
  - rawdatafolder
  - resultfolder
  - pairitel_scripts_path
  - ra_cluster
  - dec_cluster


Step 2

Open a python shell in your home directory. Run 'pairitel_pyraf_0.py' in this shell. This shell now has pyraf loaded and will referred to as the "pyraf shell" from now on.
For my CfA system: use /usr/local/bin/ipython --colors lightbg

DO THIS:
execfile('/pathtofile/pairitel_pyraf_0.py')

The 'pairitel_pyraf_0.py' script will do some basic operations on the files - it will copy everything you need into the specified resultfolder, and will trim, exposure-normalize, and coordinate-shift the images.


Step 3

Run 'pairitel_pyraf_1.py' in the pyraf shell. This script will do the aperture photometry for all nights.

DO THIS:
execfile('/pathtofile/pairitel_pyraf_1.py')

(By the way: you can check the results of the aperture photometry by loading the .coo.1 file of an image into DS9 (load region - selection: all - (click file) - format: X Y - coordinate system: image).)

The script will tell you which night yielded the sharpest images, it also tells you in which image the most sources were found by aperture photometry. Look at those images by hand (for example with ds9) and choose one of them as your master image for the source extraction. Save that filename, including the path, in the input_info.py file.

DO THIS:
- choose master image file
- modify "masterimage" in input_info.py

Then open that image in ds9 (if it's not open yet) and define ca. 10 single, bright, but not saturated stars as your stars for the psf fitting. Do this by placing circular regions centered on those stars (the radius does not matter), and save those regions, using decimal fk5 coordinates. Save the filename of that file, including the path, in the input_info.py file. I usually save mine in the parent folder of the individual nights, so that I can find it easily if I want to check something by hand in ds9.
It can help to open some more images and add more psf stars to that file, especially if the central point of your nightly observations shifts a lot.

DO THIS:
- open master image in ds9
- put regions on ca. 10 stars for psf fitting
- dave region file in standard ds9 format in the "resultfolder" under a name like "psfstars.reg". Use fk5 decimal coordinates.
- load that region into images from other nights and add some more psf stars. Save region file over old filename. Make sure to use decimal coordinates again.
- modify "psfstarfile" in input_info.py


Step 4

Open another python shell in whatever directory you want. Run 'pairitel_astropy_1.py' in this shell; it loads astropy and will be referred to as the "astropy shell" from now on.

DO THIS:
execfile('/pathtofile/pairitel_astropy_1.py')

This script translates the ds9 region file which contains your stars for psf fitting into the image coordinates of each single observation and saves them as files ending in 'pstbyhand'. 


Step 5

Run 'pairitel_pyraf_2.py' in the pyraf shell.

DO THIS:
execfile('/pathtofile/pairitel_pyraf_2.py')

This script performs the psf photometry for your master image. It fits a gaussian plus an empirical lookup-table as the psf (this turned out to be the best choice for Pairitel), does an initial run of psf photometry on the master image, then subtracts all found stars from the image, detects the remaining stars, and then performs a final psf photometry for all those found sources. The final source magnitudes etc. will be saved in a file ending in 'als.3'.

After that, open the masterimage in ds9 and look at the sources which were finally found by psf photometry, stored in the region file 'stars_2.reg'. Add stars by hand if any obvious ones were missed, and save the region file in the resultfolder, naming it something obvious like 'masterstars.reg'; use WCS and decimal coordinates for saving the file. Put the file name into the input_info.py file under 'masterregfile'.

DO THIS:
- open master image in ds9.
- open the *.psf.2.fits file which belongs to the master image in ds9 and check if there are faint stars in the psf wing. If not, good. If yes, check your psf star regions again and throw out psf stars which have close neighbors; repeat steps 4 and 5.
- open the *.sub.2.fits file which belongs to the master image in ds9. All detected sources have been subtracted from the image, using the psf model. Check if there are "holes" at the source positions; usually there are. Check if the holes are in the few % range of the actual source in the master image. They should be.
- Then open the master image in ds9 and load the region file stars2.reg from the same folder. These are all sources which will be extracted from the other nights. Check by hand if there are sources missed, or if some extra sources have been detected which are not really there. Delete or add circular regions for that. Save the cleaned up file in the "resultfolder" under a name like "masterstars.reg".
- modify "masterregfile" in input_info.py


Step 6

Run 'pairitel_astropy_2.py' in the astropy shell.

DO THIS:
execfile('/pathtofile/pairitel_astropy_2.py')

This script takes the masterregfile produced for the master image and transforms it to coordinate files (ending in '.mastercoo.1') for all the individual images. This will be the position list for the source extraction in all images.


Step 7

Run 'pairitel_pyraf_3.py' in the pyraf shell.

This script performs the psf photometry for all other images, using the master source list derived from the master image. After doing this, blink all final psf fits in ds9 and see if everything looks good. Also blink all .sub.2.fits files; those are the files in which all found sources are subtracted, and those files should not have huge holes in them. If there are files in which no stars were subtracted or weird stuff happened, try a smaller radius for input_info.psfcleaningradius.


Step 8

Run 'pairitel_astropy_3.py' in the astropy shell.

This script collects the 2MASS catalogue centered on the masterimage, and will be used in the following steps to calibrate the Pairitel magnitudes.


Step 9

Run 'pairitel_astropy_4.py' in the astropy shell.

This script collects the raw magnitudes of all sources in all nights and combines them into (uncalibrated) light curves for each source. All light curves are stored in a single file, where each source is given a unique identifier (similar to the YSOVAR2 database). This file is by default:
input_info.resultfolder + 'rawlcs.dat'


Step 10

Run 'pairitel_astropy_5.py' in the astropy shell.

This script uses the infrastructure provided by the pYSOVAR package; the uncalibrated light curves are transformed into an atlas object which adds some nice functionalities the script will use.
The matching 2MASS magnitudes are found for sources which have 2MASS counterparts. 
If the user has supplied a catalog file which lists class 1 and class 2 sources, those identifiers are also added (see catalogfile etc. in input_info.py). 
The calibrated the Pairitel data to the 2mass data in two steps:
1. All sources with 2MASS counterparts are used. All magnitudes are shifted linearly to the 2MASS amgnitudes, and the standard deviation of (Pairitel_calibarted - 2MASS) is added to the photometric errors as the systematic error induced by the fit. This will overestimate the true errors: In this fit, there will be intrinsically variable sources which induce a large scatter. 

There will be a figure displayed which shows a histogram of the scatter of the now calibrated light curves. Choose a threshold from this; you want mostly contant light curves for the second step of the calibration. But you also want to choose enough quasi-constant sources to make a meaningful fit, so go for something slightly to the left of the peak. Put this value into input_info.py (threshold_lc).


Step 11

Run 'pairitel_astropy_6.py' in the astropy shell.

This performs the second round of calibration to the 2MASS data:
2. The script looks at all calibrated light curves created in the previous step, takes the ones which are less variable than the user-defined threshold, and re-calibrates the Pairitel data using only those (mostly non-variable) objects. Then the scatter will be much smaller, but still somewhat over-estimated because some of the calibration sources may still have some low-level intrinsic variability.





