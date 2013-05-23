Recipe for Pairitel data reduction
----------------------------------

At CfA, I have to use two different python shells because Astropy and Pyraf cannot be loaded in the same python shell here. If you can load both pyraf and astropy in the same shell on your system, just run all the scripts from the same shell.

Step 1

Define the central coordinates of your cluster and several other constants in the input_info.py file. That file also gives some info on how your data should be stored.


Step 2

Open a python shell in your home directory. Run 'pairitel_pyraf_0.py' in this shell. This shell now has pyraf loaded and will referred to as the "pyraf shell" from now on.
For my CfA system: use /usr/local/bin/ipython --colors lightbg

The 'pairitel_pyraf_0.py' script will do some basic operations on the files - it will copy everything you need into the specified resultfolder, and will trim, exposure-normalize, and coordinate-shift the images.


Step 3

Run 'pairitel_pyraf_1.py' in the pyraf shell. This script will do the aperture photometry for all nights and tell you which night yielded the sharpest images, it also tells you in which image the most sources were found by aperture photometry. Look at those images by hand (for example with ds9) and choose one of them as your master image for the source extraction. Save that filename, including the path, in the input_info.py file.

Then open that image in ds9 (if it's not open yet) and define ca. 10 single, bright, but not saturated stars as your stars for the psf fitting. Do this by placing circular regions centered on those stars (the radius does not matter), and save those regions, using decimal fk5 coordinates. Save the filename of that file, including the path, in the input_info.py file. I usually save mine in the parent folder of the individual nights, so that I can find it easily if I want to check something by hand in ds9.
It can help to open some more images and add more psf stars to that file, especially if the central point of your nightly observations shifts a lot.


Step 3

Open another python shell in whatever directory you want. Run 'pairitel_astropy_1.py' in this shell; it loads astropy and will be referred to as the "astropy shell" from now on.

This script translates the ds9 region file which contains your stars for psf fitting into the image coordinates of each single observation and saves them as files ending in 'pstbyhand'. 
execfile('/data/swolk/kpoppen/Dropbox/MyPython/PAIRITEL/pairitel_astropy_1.py')


Step 4

Run 'pairitel_pyraf_2.py' in the pyraf shell.

This script performs the psf photometry for your master image. It fits a gaussian plus an empirical lookup-table as the psf (this turned out to be the best choice for Pairitel), does an initial run of psf photometry on the master image, then subtracts all found stars from the image, detects the remaining stars, and then performs a final psf photometry for all those found sources. The final source magnitudes etc. will be saved in a file ending in 'als.3'.

After that, open the masterimage in ds9 and look at the sources which were finally found by psf photometry, stored in the region file 'stars_2.reg'. Add stars by hand if any obvious ones were missed, and save the region file in the resultfolder, naming it something obvious like 'masterstars.reg'; use WCS and decimal coordinates for saving the file. Put the file name into the input_info.py file under 'masterregfile'.

Step 5

Run 'pairitel_astropy_2.py' in the astropy shell.

This script takes the masterregfile produced for the master image and transforms it to coordinate files (ending in '.mastercoo.1') for all the individual images. This will be the position list for the source extraction in all images.


Step 6

Run 'pairitel_pyraf_3.py' in the pyraf shell.

This script performs the psf photometry for all other images, using the master source list derived from the master image. 


Step 7

Run 'pairitel_astropy_3.py' in the astropy shell.

This script collects the 2MASS catalogue centered on the masterimage, and will be used in the following steps to calibrate the Pairitel magnitudes.


Step 8

Run 'pairitel_astropy_4.py' in the astropy shell.

This script matches the Pairitel sources to nearby 2MASS sources for each image and applies a linear correction to match the Pairitel magnitudes with the 2MASS ones. The result is written into 3 files with the names 'calibratedmags_(band).dat' in each nightly folder.


Step 9

Run 'pairitel_astropy_5.py' in the astropy shell.

This script uses the coordinates saved in the masterreg file to bundle the objects (necessary because I couldn't get the pyraf object IDs to match); this yields three files, one for each band, which contain [object_id ra dec time band_magnitude band_error] similar to the YSOVAR database, i.e. multiple lines for one object which have the same object_ID.
These files can be added to the YSOVAR_atlas object class as usual.






