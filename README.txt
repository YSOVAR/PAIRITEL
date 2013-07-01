Recipe for Pairitel data reduction
----------------------------------

At CfA, I have to use two different python shells because Astropy and Pyraf cannot be loaded in the same python shell. If you can load both pyraf and astropy in the same shell on your system, just run all the scripts from the same shell.

Step 1

Rename some files and provide some basic information about the paths where your data is and some info about your cluster.

DO THIS:
- Copy the input_info_sample.py file to input_info.py
- Copy the wcs_by_hand_sample.txt file to wcs_by_hand.txt
- If you want to add stuff to this file, rename this file to README_yourname.txt or something. Otherwise it will be overwritten when you pull a new version of the PAIRITEL pipeline.
- Modify in input_info.py:
  - rawdatafolder
  - resultfolder
  - pairitel_scripts_path
  - ra_cluster
  - dec_cluster


Step 2

The script will copy all data files into an output folder, set some iraf parameters, and add a readout noise keyword to the pairitel headers because that is missing by default.

DO THIS:
- open ipython shell in home directory. (This shell will be called "pyraf shell" from now on.) (For my CfA system: use /usr/local/bin/ipython --colors lightbg)
- Run script pairitel_pyraf_0a.py:
execfile('/pathtofile/pairitel_pyraf_0a.py')


Step 3

In this step, the coordinate shifts in the sky images are corrected.

DO THIS:
- In the pyraf shell, execute pairitel_pyraf_0b.py script:
execfile('/pathtofile/pairitel_pyraf_0b.py')
- check manually if the coordinates have been matched successfully:
- open a new shell
- go to your pairitel "output" folder
- type: ds9 -scale log YSO*/*_wcs_unnormed.fits &
- click into the first frame
- click in the menu: Frame > Single Frame
- zoom and center this frame in a way that you find convenient. 
- click in menu: Edit > Crosshair
- put the crosshair on some central, bright star in your cluster that you will easily recognize in the other images.
- click in menu: Frame > Match > Frame > WCS
- click in menu: Frame > Lock > Crosshair > WCS
- click in the "tab" bar on "frame", then "next" to check out all the other frames. If you happen to mis-click somewhere and you shifted a frame, just match the frame/WCS and lock the crosshair/wcs again.
- Are all frames neatly aligned? if YES, skip Step 4 and go to Step 5.
- if NO: go to Step 4. Leave images in ds9 open.


Step 4

In this step we correct some coordinate systems by hand which the automatic coordinate correcting algorithm didn't get right.

DO THIS:
- Open the file wcs_by_hand.txt in a text editor. You see some sample entries which you will substitute with your own entries. You will have one line for every frame that is mismatched.
- In ds9, go to a frame where the coordinates are good. Get the (decimal) coordinates of one star which you will recognize, for example the obe on which you put the crosshairs.
- To get decimal coordinates, click in menu: WCS > Degrees
- Copy the ra and dec of that star into the ra_good and dec_good column in wcs_by_hand.txt. Do this for every line (all the same coordinates).
- In ds9, click frame/next until you reach a frame that is mismatched.
- Find the star you had marked before and get its x and y pixel values. Put those values into a line in wcs_by_hand.txt; put them into the columns x_pixel and y_pixel.
- If you want to try and correct this file (if you think it's not complete garbage), put "no" into the column "beyond_saving" in wcs_by_hand.txt.
- Then put the filename of the mismatched image, including the path, in the column "filename_incl_path" in wcs_by_hand.txt.
- In ds9, go back to the last good image, put the crosshair on your reference star again, and do lock/crosshairs again.
- Click frame/next until you find the next mitmatched image and put the info into wcs_by_hand.txt. If you find a frame that you  think is beyond saving, put its filename in wcs_by_hand.txt and set "beyond_saving" to "yes". The x_pixel and y_pixel values can then be arbitrary.
- When you're done, save the wcs_by_hand.txt file.
- Then go to your pyraf shell and execute the script pairitel_pyraf_0c.py
execfile('/pathtofile/pairitel_pyraf_0c.py')
- Look at all images again in ds9 using: ds9 -scale log YSO*/*_wcs_unnormed.fits &
- If there are still some unmatched files, repeat. But it should be fine now.


Step 5

In this step, the coordinate shifts performed on the images are also applied to the exposure masks. The images are then trimmed so that low-exposure edges are thrown away, and then normalized by the exposure mask.

DO THIS:
- in the pyraf shell, run the pairitel_pyraf_0d.py script
execfile('/pathtofile/pairitel_pyraf_0d.py')
- If you get weird errors like "xxx not an image or a number", close python, open a new shell in yout home directory, start a fresh ipython there and run this script again. (This is a strange IRAF problem with some of the IRAF versions.)


Step 6

This step performs the preliminary aperture photometry on the images.

DO THIS:
- Run the pairitel_pyraf_1.py script in the pyraf shell.
execfile('/pathtofile/pairitel_pyraf_1.py')
- In the pyraf shell, you will see some text which tells you in which images the most sources were found.
- Open a few of those images in ds9. You will need to choose one of them as your "masterimage". It should be an image with many detected sources, which should be roughly centered on the region you are most interested in.
- You can overplot the found sources in ds9 by clicking: region > load > All (bottom right) > choose the .coo.1 file which belongs to that image > Format: X Y > Coordinate System: image > OK.
- The detection threshold in the preliminary run is very low, so you will pick up some spurious sources. But that's okay.
- When you decided on a master image, specify that image in input_info.py in the line "masterimage".

Step 7

Now you pick some psf fitting stars by hand.

DO THIS:
- Open the master image in ds9.
- put some circular regions on ca. 10 moderately bright, not saturated, single stars. Be careful not to pick stars which are somehow elongated because of an unresolved companion. (The radius of the regions doesn't matter.)
- Save the region file in your pairitel "output" folder. Click: region > (go to output folder) > name: "psfstars.reg" or similar > OK > Format: ds9 > Coordinate System: fk5/Degrees.
- Open a few more of the images with many detected sources in ds9, load the psfstars.reg file you just created, add some more sources, save again using the same name and decimal coordinates.
- Put the filename of the psf star region file into input_info.py in the line "psfstarfile" (and save).


Step 8

This script translates the psfstar.reg region file into a pixel position table for the masterfile observation (IRAF needs X/Y coordinates for the next step.)

DO THIS:
- Type "import astropy", "import aplpy", and "import matplotlib.pyplot as plt" into the pyraf shell. If it doesn't complain, proceed with this shell. If it complains, open a new ipython shell and proceed with that shell (which will be called "astropy shell" from now on.)
- Execute the pairitel_astropy_1.py script.
execfile('/pathtofile/pairitel_astropy_1.py')
- It is possible that this script crashes at some point. If it does, look at the shell output which will tell you with which *_wcs.fits file it crashed. Look at that file with ds9. You either did not define enough psf stars so that none are found for this image; if that's the case, just add more psf stars to the psfstars.reg file, save it, and run pairitel_astropy_1.py again. Other possibility: that image is so small that it does not make sense to reduce it at all. If that's the case, go to input_info.py, set "bad_files_exist" to "yes" and put the filename into "bad_exposures" as shown. Repeat Step 8. If it crashed at a different file, repeat these repair steps until Step 8 runs without problems.
- After running successfully, the script also gives you output in the shell about other files that you should add to the "bad_exposures", even if they did not cause the script to crash.


Step 9

In this step you produce the psf photometry for the masterimage.

DO THIS:
- Run 'pairitel_pyraf_2.py' in the pyraf shell.
execfile('/pathtofile/pairitel_pyraf_2.py')
- open the following files (located in the folder of the masterimage) in ds9:
  - (masterimage).psf.1.fits
  - (masterimage).psf.2.fits
  - (masterimage).sub.2.fits
  - and the masterimage itself.
- The psf.1 and psf.2 files are the first and the second try to fit a good psf. The images show the deviation from a gaussian. You should see some black and white regions in the middle, and uniform grey in the outer parts. - If you see faint companion stars in the wings of psf.1, they should be gone in psf.2. If you still see some in psf.2, go back to the definition of the psf stars, throw out stars that have neighbors, and repeat everything up to here.
- The sub.2.fits file is a file in which all detected sources have been fitted with the psf model and subtracted form the image. 
- This file should show no huge holes or bright spots. If it does, compare the depth of the holes or spikes with the actual pixel values for that star in the master image. 
- You can do that by placing a region on it and clicking Analysis > Funtools > Counts in region. If the holes or spikes are only a few percent of the masterimage, then that's okay. If not, you need to fiddle with the psf fitting parameters. Go ask Katja about this.
- Caveat: In regions where the background is very inhomogeneous, your results can be worse in individual cases, with errors of up to 20%. This should be rare, but cannot be avoided completely in this semi-automatic reduction.


Step 10

Define the master source list.

DO THIS:
- open the master image in ds9.
- load the region "stars2.reg" from the same folder.
- These are all the sources for which light curves will be constructed. Sources which are not detected in this image will not be extracted in other images.
- You now have to clean this region file by hand. Look closely at each part of it, zoom in and play with the contrast. 
  - If you think a source was missed, carefully place a new region on this source. 
  - If you think a nonexistent source was detected because of some psf-subtraction artifacts, delete that region. 
  - If you think a background fluctuation has been picked up as a source, leave it in. It won't be detected in other nights and thus not contribute to the light curve in the end.
- Save the cleaned-up region file in the pairitel "output" folder, using a convenient name like 'masterstars.reg'. Use Format=DS9, Coordinates: fk5/Degrees. 
- Put this filename into the input_info.py file in the line "masterregfile".


Step 11

This script takes the masterregfile produced for the master image and transforms it to coordinate files (ending in '.mastercoo') for all the individual images. This will be the position list for the source extraction in all images.

DO THIS:
- Go to astropy shell, run 'pairitel_astropy_2.py' in the astropy shell.
execfile('/pathtofile/pairitel_astropy_2.py')



Step 12

This script performs the psf photometry for all other images, using the master source list derived from the master image. After doing this, blink all final psf fits in ds9 and see if everything looks good. Also blink all .sub.2.fits files; those are the files in which all found sources are subtracted, and those files should not have huge holes in them. If there are files in which no stars were subtracted or weird stuff happened, try a smaller radius for input_info.psfcleaningradius and/or ask Katja about the psf fitting parameters.

DO THIS:
- Run 'pairitel_pyraf_3.py' in the pyraf shell.
execfile('/pathtofile/pairitel_pyraf_3.py')
- go to pairitel "output" folder and open all psf.2 files in ds9:
ds9 YSO*/*psf.2.fits &
- check if a psf looks weird.
- go to pairitel "output" folder and open all sub.2 files in ds9:
ds9 YSO*/*sub.2.fits &
- click through the frames to see if a sub file looks weird.
- There will be some unsubtracted stars in many images, because they were not in the master stars list.


Step 13

DO THIS:
- Run 'pairitel_astropy_3.py' in the astropy shell.
execfile('/pathtofile/pairitel_astropy_3.py')

This script collects the 2MASS catalogue centered on the masterimage, and will be used in the following steps to calibrate the Pairitel magnitudes.


Step 14

DO THIS:
- Run 'pairitel_astropy_4.py' in the astropy shell.
execfile('/pathtofile/pairitel_astropy_4.py')
- This will take a while.

This script collects the raw magnitudes of all sources in all nights and combines them into (uncalibrated) light curves for each source. All light curves are stored in a single file, where each source is given a unique identifier (similar to the YSOVAR2 database). This file is by default:
input_info.resultfolder + 'rawlcs.dat'


Step 15

This script uses the infrastructure provided by the pYSOVAR package; the uncalibrated light curves are transformed into an atlas object which adds some nice functionalities the script will use.
The matching 2MASS magnitudes are found for sources which have 2MASS counterparts. 
In this step, all sources with 2MASS counterparts are used for the calibration. All magnitudes are shifted linearly to the 2MASS magnitudes, and the standard deviation of (Pairitel_calibarted - 2MASS) is added to the photometric errors as the systematic error induced by the fit. This will overestimate the true errors: In this fit, there will be intrinsically variable sources which induce a large scatter. 

DO THIS:
- Run 'pairitel_astropy_5.py' in the astropy shell.
execfile('/pathtofile/pairitel_astropy_5.py')
- Look at the last output figure. It is a histogram of the standard deviation of the first-round calibration light curves. You  want to identify the light curves which are mostly contant, i.e. have a low standrad deviation:
- pick a threshold, somewhere near the peak of the histogram.
- put that value into input_info.py (line "threshold_lc").


Step 16

This step re-does the calibration, using only mostly constant sources. This will yield a better fit and a smaller systematic error introduced by the fit. The script will also display a random light curve from your cluster.

- Run 'pairitel_astropy_6.py' in the astropy shell.
execfile('/pathtofile/pairitel_astropy_6.py')





