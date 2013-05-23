# -*- coding: utf-8 -*-
'''This module contains functions to help with IRAF/DAOPHOT photometry.

It has to be run in PYRAF to call the relevant IRAF commands.
'''


import urllib
import StringIO
from pyraf import iraf
import numpy as np
import glob
import os
import string
import shutil
import sys
import asciitable
#import catalog
from atpyextensions import catalog
import pyfits
from copy import deepcopy

paramdir=''

def n_extensions(fitsfile):
  '''returns number of extension in a fits file
  
  Usage: n = n_extensions(filename)'''
  iraf.mscred(_doprint=0)
  result=iraf.mscextensions(fitsfile, output='none')
  return iraf.mscextensions.nimages

def fits_ext(fitsfile):
    '''returns a list of strings in the form [XXX.fits[1], XXX.fits[2]]
    
    This function retunrs a list of strings, that IRAF can use a filenames to access all
    extension of a given fits file. '''
    return [fitsfile+'['+str(i)+']' for i in range(1,1+n_extensions(fitsfile))]

def get_last_iraf(prefix,suffix=''):
    '''Return the file with highest version number of IRAF/DAOPHOT's numbering scheme
    
    The function returns the filename as string
    Example:
    
    > ls
    0108.Land01.fits2.mag.1.sex  0108.Land01.fits2.mag.2.sex  0108.Land01.fits2.mag.3.sex
    > python
    >>>get_last_iraf('0108.Land01.fits2.mag.',suffix='sex')
    '0108.Land01.fits2.mag.3.sex'
    
    Prefix and suffix can be supplied with or without trailing / leading "." '''
    prefix=prefix.rstrip('.')
    prefix=prefix.replace('[','')
    prefix=prefix.replace(']','')
    suffix=suffix.lstrip('.')
    if len(suffix) > 0: suffix='.'+suffix #make suffix begin with '.' if not empty
    nl=glob.glob(prefix+'*')
    if len(suffix) == 0:
        nl=filter(lambda s:s[len(prefix)+1:].isdigit(),nl)
    else:
        nl=filter(lambda s:s[len(prefix)+1:-len(suffix)].isdigit(),nl)
        nl=map(lambda s:s[:-len(suffix)],nl) #cut off suffixes, so the case with and withour suffix can be treated consostenstly from here on
    if nl == []:
        return None
    else:
        return prefix+'.'+str(max(map(lambda s: int(s[len(prefix)+1:]),nl)))+suffix



def get_next_iraf(prefix,suffix=''):
    '''Return the file with next version number of IRAF/DAOPHOT's numbering scheme
    
    The function returns the filename as string
    Example:
    
    > ls
    0108.Land01.fits2.mag.1.sex  0108.Land01.fits2.mag.2.sex  0108.Land01.fits2.mag.3.sex
    > python
    >>>get_next_iraf('0108.Land01.fits2.mag.',suffix='sex')
    '0108.Land01.fits2.mag.4.sex'

    Prefix and suffix can be supllied with or without trailing / leading "." '''
    prefix=prefix.rstrip('.')
    prefix=prefix.replace('[','')
    prefix=prefix.replace(']','')
    suffix=suffix.lstrip('.')    
    last_iraf_file=get_last_iraf(prefix,suffix=suffix)  #make suffix begin with '.' if not empty
    if len(suffix) > 0: suffix='.'+suffix
    if last_iraf_file == None:
        return prefix+'.1'+suffix
    else:
        if len(suffix) == 0:
            return prefix+'.'+str(int(last_iraf_file[len(prefix)+1:])+1)
        else:
            return prefix+'.'+str(int(last_iraf_file[len(prefix)+1:-len(suffix)])+1)+suffix

def get_sky(image):
    '''return estimate for sky value and standard deviation of sky
    
    implicit assumption is that the sky is constant throughout the image
    and that most of the image actually shows the sky (no very crowded fields)!
    calls iraf.imstat with ncli=10
    this calculates mean and standard dev of mean, clips all values significantly above
    the standard and then recalculates mean and std dev. This loop runs 10 times to 
    clip of values of the stellar counts. 
    stat=iraf.imstat(image,fields='midpt,stddev',ncli=10,Stdout=1)
    return float(stat[1].split()[0]), float(stat[1].split()[1])

    Usage: sky,skydev=get_sky(filename)''' 
    stat=iraf.imstat(image,fields='midpt,stddev',ncli=10,Stdout=1)
    
    return float(stat[1].split()[0]), float(stat[1].split()[1])

def get_catalogue(filename,ra,dec,radius=20.,catalogue='tmc',minmag=100.):
  '''query a cataloge for positons around a specific aimpoint.
  
  This function send a query to a cgi script at NOAO to query the scat tool there
  filename=output
  ra,dec: coordinates can be decimal or sexagesimal
  radius: in arcmin, default=20.
  catalog: identifier (default tmc for 2MASS)
  minmag: to filter out objects fainter than this limit (default=100) 
  '''
  if catalogue !='tmc': raise NotImplementedError
  data=urllib.urlopen('http://archive.tuc.noao.edu/cgi-bin/scat?catalog='+catalogue+'&ra='+str(ra)+'&dec='+str(dec)+'&sys=J2000&mrad='+str(radius)+'&nstar=-1&mag='+str(minmag)).read()
  #mask out comments, change order of columns? so that iraf.msccmatch() understands the file
  webfile=StringIO.StringIO(data)
  for line in range(10):
    dump=webfile.readline()  #discard header
  # Pyhton 2.6 syntax for open/close with open(filename, 'w') as file:
  file=open(filename, 'w')
  for line in webfile: 
    splitline=line.split('\t')
    file.write(splitline[1]+' '+splitline[2]+' '+splitline[0]+"\n")
  file.close()


def correct_wcs(image,clobber=False,minmag=13, radius=20.):
  '''correct wcs of fitsfile to 2MASS coordinates
  
  This method loads the 2MASS catalog in the region covered by the fits file
  and shifts,stretches, and rotates the wcs in the fits file to match the 2MASS coordinates
  The initial wcs needs to be approximately correct (within a few arcsec), otherwise 
  the matching will fail.
  image = filename
  clobber = overwrite existing 2MASS catalog file?
  minmag = mag of faintest object retrieved from 2MASS'''
  iraf.mscred(_doprint=0)
  object = iraf.hselect(image+'[0]','OBJECT', 'yes', Stdout=1)[0]
  filter = iraf.hselect(image+'[0]','FILTER', 'yes', Stdout=1)[0]
  catalogue=os.path.join(os.path.dirname(image),"2mass_"+str(object)+'_'+str(filter)+'_'+str(minmag)+".cat")
  #Gererate catalogue, if it does not exist in this directory
  if clobber or not os.path.isfile(catalogue):
    get_catalogue(catalogue,iraf.hselect(image+'[0]','RA','yes', Stdout=1)[0],iraf.hselect(image+'[0]','DEC','yes', Stdout=1)[0], radius=radius, minmag=minmag)
  #cfrac: maximum fraction of targets which do not center correctly
  #       a large number is no problem here, since often there are a few 100 sources in the catalog to match.
  matchout=iraf.msccmatch(image,catalogue,interactive=False, Stdout=1,cfrac=0.95,rms=3)
  try: print matchout[-5]
  except IndexError: print matchout

def shift_wcs(image, shifttab, number):
  '''Shift wcs in image according to value in shifttab
  
  image: filename
  shifttab: structured record array with cols "from to rashift decshift"
  rashift and decshift in arcsec
  number: ID of image to chose record in shifttab
  
  This method applies the rashift and decshift specified in the record of shifttab with
  from <= number <= to
  
  from to rashift decshift
  90  105  46.7  16.0
  106 111  48.0  21.0
  112 127 -29.6  18.8
  '''
  iraf.mscred(_doprint=0)
  index=(shifttab['from']<=number)*(shifttab['to']>=number)
  iraf.mscwcs(image,ra_shift=shifttab['rashift'].compress(index)[0],dec_shift=shifttab['decshift'].compress(index)[0])

def set_Pairitel_params():
  '''Set some daophot parameters, which are typical for my Pairitel images
  '''
  iraf.digiphot(_doprint=0)
  iraf.daophot(_doprint=0)
  iraf.apphot(_doprint=0)
  iraf.datapars.fwhmpsf = 4.5
  iraf.datapars.datamax=15000.
  iraf.datapars.ccdread=10. # 30. in 2MASS?
  iraf.datapars.gain = "GAIN"
  iraf.datapars.exposure = "EXPTIME"
  iraf.datapars.airmass = "AIRMASS"
  iraf.datapars.filter = "FILTER"
  iraf.datapars.obstime="DATE-OBS"



def set_FLWO_params():
  '''Set some daophot parameters, which are typical for my FLWO images
  '''
  iraf.digiphot(_doprint=0)
  iraf.daophot(_doprint=0)
  iraf.datapars.fwhmpsf = 4.5
  iraf.datapars.datamax=40000.
  iraf.datapars.ccdread = "rdnoise"
  iraf.datapars.gain = "gain"
  iraf.datapars.exposure = "exptime"
  iraf.datapars.airmass = "air"
  iraf.datapars.filter = "filter"
  iraf.datapars.obstime="date-obs"





def aperture_photometry(image,threshold=4.,wcs='logical',mode='small', telescope='Pairitel'):
  '''perform aperture photometry on a given image
  
  This procedure perfrom IRAF/DAOPHOT aperture phtometry on (all extensions of) a fits image.
  It sets some global daophot parameters, then perfroms daophot.daofind and daophot.phot
  The output .mag file is filtered for valid magnitudes and a mag.sex file computed, where
  all coordinates are transformed to sexagesimal coordinates.
  
  input: image: filename
  keywords:  threshold=4 :daofind detection threshold in sigma_sky
             wcs='logical' :working wcs system, can be  'logical','physical' or 'tv'
                        in any case an additional sexagesimal output file will be computed
         mode='small' : 'standard' = 'large' or 'psf' = 'small' choses a pre-defined aperture size for 
                        aperture photometry of standards or just as starting point for psf photometry
  '''
  print 'Performing aperture photometry on '+image
  if telescope == 'Pairitel':
      set_Pairitel_params()
  elif telescope == 'FLWO':
      set_FLWO_params()
  else:
      print 'No parameters found for this telescope!'
  
  iraf.mscred(_doprint=0)
  iraf.digiphot(_doprint=0)
  iraf.daophot(_doprint=0)
  iraf.daophot.verify = False
  iraf.daophot.wcsin = wcs
  iraf.daophot.wcsout = wcs
  iraf.daophot.wcspsf = wcs
  iraf.centerpars.cbox=2.*iraf.datapars.fwhmpsf
  iraf.fitskypars.annulus=100 # usuaaly one would use 5.*iraf.datapars.fwhmpsf
  #annulus must be wide enough to have good sky statistics
  iraf.fitskypars.dannulus=10. # this is standard
  if 'standard'.find(mode.lower()) != -1 or 'large'.find(mode.lower()) != -1:
    iraf.photpars.aperture = 4.*iraf.datapars.fwhmpsf
    iraf.centerpars.calgorithm = 'centroid'
  elif 'small'.find(mode.lower()) != -1 or 'psf'.find(mode.lower()) != -1:
    iraf.photpars.aperture = 4. # close in for crowded regions
    iraf.centerpars.calgorithm = 'none'
  else:
    print '''use: aperture_photometry(image,threshold=4.,wcs="logical",mode="small")
    
    The following mode keywords are implemented: standard, large, psf, small'''
    raise NotImplementedError
  #for imagebi in fits_ext(image):
  sky,skydev=get_sky(image)
  iraf.datapars.datamin=sky-5*skydev
  iraf.datapars.sigma=skydev      
  #this could be changed to use Stdin and Stdout so that fewer files are written (saves space and time), but for debug purposes better keep all that
  iraf.daofind(image, output='default', verbose=False, verify=False, threshold=threshold) 
  iraf.daophot.phot(image, coords='default', output='default', verbose=False, verify=False, interactive=False, wcsout=wcs)
  # MAG!=INDEF is stronger than PERROR=="NoError", because some stars (those with BigShift have no error, but still INDEF as MAG
  # want to gt rid of thoses as well - they are artifacts on the chip boundary anyway 

    
###  find transformation coefficients ###




def ds9_lines(filename, ra1, dec1, ra2, dec2, title):
    '''makes ds9 region files, which contain lines connecting 2 points
    
    All input lists must have the same length, because each match is percieved as
    filename, ra1[i],dec1[i],ra2[i],dec2[i],title[i]
    
    filename:
    ra1, dec1: ra and dec of point one in any format ds9 can read
    ra2, dec2: same for point 2
    tile: string to be printed on the line connecting point 1 and 2
    '''
    file=open(filename, 'w')
    file.write('# Region file format: DS9 version 4.1\n')
    file.write('# Filename ' + filename + '\n')
    file.write('global color=blue dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    file.write('fk5\n')
    for i in range(len(ra1)):
        file.write('line('+str(ra1[i])+','+str(dec1[i])+','+str(ra2[i])+','+str(dec2[i])+') # line=0 0 text={'+str(title[i])+'}\n')
    file.close()


def psf_photometry_pairitel(imagebi, satmag, psfstarfile='', ds9=False):
    '''This method performs psf photometry on a single fits extension
    
    It first sets some global daophot parameters, selects psf stars, computes the psf, cleans
    the psf stars from close neighbours and recomputes te psf and then fits all sources in the field.
    Existing aperture photometry (.coo and .mag files) in the directory is a PREREQUISITE for this routine.
    
    keyword: ds9=False: If true the method writes a files with X Y coords of all succesfully fitted sources in the current WCS 
    '''
    print imagebi
    iraf.daophot.verify = False
    #needs to be set here because it is different for each amplifier
    sky,skydev=get_sky(imagebi)
    iraf.datapars.datamin=sky-5*skydev
    iraf.datapars.sigma=skydev
    #bad pixels and image edges are not allowed in fitrad, but generate a warning only in psfrad
    #to avoid all problems associated with that use fitrad only in the initial selection of psf stars
    if psfstarfile == '':
        # let iraf select stars for psf fitting.
        out=iraf.pstselect(imagebi, 'default', 'default', maxnpsf=10, psfrad=0, fitrad=iraf.daopars.psfrad, Stdout=1)
        pstcat=catalog.CartesianCatalog(get_last_iraf(imagebi+'.pst'), type='daophot', RA='XCENTER', DEC='YCENTER')
        # The following line works because Python short-circuits bolean expressions
        # If the first is true, then pst cat is [], so it is important that pstcat['MAG'] is not called
        print pstcat
        # filter the pst file so that stars which are too bright are deleted. 
        psf_star_checking(get_last_iraf(imagebi+'.pst'), satmag)
        # now fit a psf model to the cleaned list of psf stars.
        iraf.daophot.psf(imagebi,photfile='default',pstfile=get_last_iraf(imagebi+'.pst'),psfimage='default',opstfile='default',groupfile='default', interactive=False,verify=False,nclean=10)
    else:
        # if a list of psf stars was selected by hand, just load that file and do the psf fitting to those stars without further ado.
        iraf.daophot.psf(imagebi,photfile='default',pstfile=imagebi+'.pstbyhand',psfimage='default',opstfile='default',groupfile='default', interactive=False,verify=False,nclean=10)
        psf_star_checking(get_last_iraf(imagebi+'.pst'), satmag)
        
    #use a smaller radius for nstar and substar to substract the "core only" of close neighbours
    original_psfrad=iraf.daopars.psfrad
    iraf.daopars.psfrad = 15
    # find neighbors of psf stars which are in the fitting radius.
    iraf.nstar(imagebi,get_last_iraf(imagebi+'.psg'),'default','default','default', verbose = False)
    # subtract the first-round psf model of those neighboring stars from the psf stars (so that wings of psf stars are star-free).
    out=iraf.substar(imagebi,get_last_iraf(imagebi+'.nst'), get_last_iraf(imagebi+'.pst'), get_last_iraf(imagebi+'.psf',suffix='fits'), 'default', verbose=True, Stdout=1)
    #then reset psfrad to original value
    iraf.daopars.psfrad = original_psfrad
    # now get the second-round psf model by fitting the neighbor-subtracted psf stars.
    iraf.daophot.psf(get_last_iraf(imagebi+'.sub',suffix='fits'), get_last_iraf(imagebi+'.mag'),get_last_iraf(imagebi+'.pst'), get_next_iraf(imagebi+'.psf',suffix='fits'), get_next_iraf(imagebi+'.pst'), get_next_iraf(imagebi+'.psg'), interactive=False,verify=False,nclean=10)
    # now extract all stars using the second-round psf model.
    
    iraf.daophot.allstar(imagebi,photfile=get_last_iraf(imagebi+'.mag'),psfimage=get_last_iraf(imagebi+'.psf',suffix='fits'),allstarfile='default',rejfile='default',subimage='default',verbose=False)
    
    xycenter_to_ds9reg(get_last_iraf(imagebi+'.als'), 'stars_1.reg')
    
    # take the sub.image (produced as a by-product in the previous step) where all the found sources are subtracted, and find all remaining sources which were so far hidden in the glare of brighter stars.
    # find the sources in the sub.image, using somewhat higher significance threshold so that we don't pick up residuals from the psf subtraction:
    iraf.daophot.daofind(get_last_iraf(imagebi+'.sub',suffix='fits'),'default',threshold = 6., verbose=False)
    # get preliminary magnitudes for those sources, using the original image:
    iraf.daophot.phot(imagebi, coords=get_last_iraf(imagebi+'.sub.2.fits.coo'), output=imagebi+'.secondrun.mag.1', verbose=False, verify=False, interactive=False)
    # now merge the preliminary new sources with the sourcelist from the first round:
    iraf.pfmerge(inphotfiles=get_last_iraf(imagebi+'.als')+','+imagebi+'.secondrun.mag.1', outphotfile=get_next_iraf(imagebi+'.als'))
    # and renumber the sources so that they have unique IDs:
    iraf.prenumber(get_last_iraf(imagebi+'.als'))
    #  Now extract all sources, old and new, with PSF fitting. This will shift the sources around a bit and also throw out sources that are not significantly detected.
    
    iraf.daophot.allstar(imagebi,get_last_iraf(imagebi+'.als'),get_last_iraf(imagebi+'.psf',suffix='fits'),get_next_iraf(imagebi+'.als'),'default','default',verbose=False)
    
    xycenter_to_ds9reg(get_last_iraf(imagebi+'.als'), 'stars_2.reg')


def psf_photometry_pairitel_with_coo(imagebi, satmag, psfstarfile='', ds9=False):
    print imagebi
    iraf.daophot.verify = False
    #needs to be set here because it is different for each amplifier
    sky,skydev=get_sky(imagebi)
    iraf.datapars.datamin=sky-5*skydev
    iraf.datapars.sigma=skydev
    # PSF has to be fitted for each night and band individually.
    if psfstarfile == '':
        # let iraf select stars for psf fitting.
        out=iraf.pstselect(imagebi, 'default', 'default', maxnpsf=10, psfrad=0, fitrad=iraf.daopars.psfrad, Stdout=1)
        pstcat=catalog.CartesianCatalog(get_last_iraf(imagebi+'.pst'), type='daophot', RA='XCENTER', DEC='YCENTER')
        # The following line works because Python short-circuits bolean expressions
        # If the first is true, then pst cat is [], so it is important that pstcat['MAG'] is not called
        print pstcat
        # filter the pst file so that stars which are too bright are deleted. 
        psf_star_checking(get_last_iraf(imagebi+'.pst'), satmag)
        # now fit a psf model to the cleaned list of psf stars.
        iraf.daophot.psf(imagebi,photfile='default',pstfile=get_last_iraf(imagebi+'.pst'),psfimage='default',opstfile='default',groupfile='default', interactive=False,verify=False,nclean=10)
    else:
        # if a list of psf stars was selected by hand, just load that file and do the psf fitting to those stars without further ado.
        iraf.daophot.psf(imagebi,get_last_iraf(imagebi+'.mag'),pstfile=imagebi+'.pstbyhand',psfimage='default',opstfile='default',groupfile='default', interactive=False,verify=False,nclean=10)
        psf_star_checking(get_last_iraf(imagebi+'.pst'), satmag)
    
        #use a smaller radius for nstar and substar to substract the "core only" of close neighbours
    original_psfrad=iraf.daopars.psfrad
    iraf.daopars.psfrad = 15
    # find neighbors of psf stars which are in the fitting radius.
    iraf.nstar(imagebi,get_last_iraf(imagebi+'.psg'),'default','default','default', verbose = False)
    # subtract the first-round psf model of those neighboring stars from the psf stars (so that wings of psf stars are star-free).
    out=iraf.substar(imagebi,get_last_iraf(imagebi+'.nst'), get_last_iraf(imagebi+'.pst'), get_last_iraf(imagebi+'.psf',suffix='fits'), 'default', verbose=True, Stdout=1)
    #then reset psfrad to original value
    iraf.daopars.psfrad = original_psfrad
    # now get the second-round psf model by fitting the neighbor-subtracted psf stars.
    iraf.daophot.psf(get_last_iraf(imagebi+'.sub',suffix='fits'), get_last_iraf(imagebi+'.mag'),get_last_iraf(imagebi+'.pst'), get_next_iraf(imagebi+'.psf',suffix='fits'), get_next_iraf(imagebi+'.pst'), get_next_iraf(imagebi+'.psg'), interactive=False,verify=False,nclean=10)
    
    ## Using the mastercoo file, all sources are extracted from the known positions, first only as aperture photometry to get initial magnitudes right, then using the individual psfs fitted for that night and band.
    # Do not shift star positions in these extractions. 
    iraf.daopars.recenter='no'
    iraf.daophot.phot(imagebi,coords=imagebi.replace('.fits','.mastercoo'),output='default',verbose=False)
    iraf.daophot.allstar(imagebi,get_last_iraf(imagebi+'.mag'),get_last_iraf(imagebi+'.psf',suffix='fits'),get_next_iraf(imagebi+'.als'),'default','default',verbose=False)
    iraf.daopars.recenter='yes'
    
    xycenter_to_ds9reg(get_last_iraf(imagebi+'.als'), 'stars_found_' + os.path.basename(imagebi)[0]  + '.reg')



def print_cols(filename,*args,**kargs):
      '''Prints a simple ascii table with the arrays in columns
      
      Optional keyword: header: will be printed first. If header is a list of strings
            it is printed in multiple lines.
      '''
      if len(args) < 1: raise ValueError('Need to pass in at least one column')
      file=open(filename, 'w')
      if 'header' in kargs:
            try:
                  kargs['header']=kargs['header']+''
                  file.write(kargs['header']+'\n')
            except TypeError:
                  for line in kargs['header']: file.write(line+'\n')
      ncols=len(args)
      nrows= min(map(len,args))
      for row in range(nrows):
            for col in range(ncols):
                  file.write(str(args[col][row]))
                  file.write(' ')
            file.write('\n')
      file.close()

def deletefrompst(pstfile, deletelist):
    '''Deletes all sources with given IDs from a pst file
    
    This selectively deletes stars from a pst file. deletelist is a list of
    ID numbers for those star that should be deleted. The header stays unaltered.
    
    Inputs: pstfile - filename
       deletelist: list of integers'''
    fin = open(pstfile)
    fout = open(get_next_iraf(string.rsplit(pstfile,'.',1)[0]), "wt")
    for line in fin:
        try:
            number = int(string.split(line)[0])
            if not number in deletelist: fout.write(line)
        except ValueError:
            fout.write(line)
    fin.close()
    fout.close()
    
def filter_pstfile(pstfile, magfile):
    '''Delete stars from pst file, which have neighbours of similar magnitude close by.
    
    This is a test in addition to iraf.pstselect, which sorts out those stars
    with brighter neighbours, with bad pixels etc.
    '''
    pstcat=catalog.CartesianCatalog(pstfile, type='daophot', RA='XCENTER', DEC='YCENTER')
    '''It would be simpler to read magfile directly into BaseCatalog, but in some cases the fields get so long
    that there is no space in between and asciitable cannot seperate them. txdump on the other hand
    actually read the column width from the file and extracts the fields even if they are not separated by spaces.'''
    magoutput='IMAGE,ID,XCENTER,YCENTER,XAIRMASS,IFILTER,MAG,MERR,PERROR'
    mags=iraf.txdump(magfile,magoutput,'PERROR=="NoError" && MAG!=INDEF && MERR !=INDEF',Stdout=1)
    magcat=catalog.CartesianCatalog(mags, type = 'ascii', Reader=asciitable.NoHeaderReader, names = magoutput.split(','), RA='XCENTER', DEC='YCENTER')

    pstquality = np.zeros(len(pstcat))

    for source in range(len(pstcat)):
        neighbours = magcat.allwithin((pstcat['XCENTER'][source],pstcat['YCENTER'][source]), iraf.daopars.psfrad)
        pstquality[source] = judge_pst_quality(neighbours)

    cutoff = len(pstcat) / 3  # use integers
    dellist = (pstcat['ID'][pstquality.argsort()])[:cutoff]
    
    #star special filtering by hand. This is stuff I could not come by automatically
    filttab = asciitable.read(os.path.join(paramdir,'filttab.dat'), converters={'str1': asciitable.convert_numpy('str'), 'str2': asciitable.convert_numpy('str')} )
    # Go through each line of a table
    # This is not ideal but fast to program
    for line in filttab:
        if (line['str1'] in pstfile) and (line['str2'] in pstfile):
            if line['Action'] == 'del':
                indexinpstcat = pstcat.coords.NNindex((line['X'],line['Y']), 5)
                if indexinpstcat != None:  # It could be None if that object is not in the pstcat for other reasons
                    dellist = np.hstack((dellist,pstcat.ID[indexinpstcat]))
            elif line['Action'] == 'keep': 
                index = np.ones(len(dellist), dtype='bool')
                indexinpstcat = pstcat.coords.NNindex((line['X'],line['Y']), 5)
                if indexinpstcat != None:  # It could be None if that object is not in the pstcat for other reasons
                    index[np.where(dellist == pstcat.ID[indexinpstcat])[0]] = False
                    dellist = dellist[index]
            elif line['Action'] == 'ignore':
                return 0
            else:
                raise ValueError(os.path.join(paramdir,'filttab.dat')+ ' contains unrecognised instructions')
    deletefrompst(pstfile, dellist)
    return len(pstcat)-len(dellist)


def judge_pst_quality(pstcat):
    '''assign a number for the usability of a star for psf fitting
    
    The function computes a number based on the distance and the magnitude difference 
    of the psf star to the brightest neighbour within iraf.daophot.psfrad
     
     A bigger number indicates a better quality, i.e. the star is much brighter 
     than all neighbours and the next brightest neighbour is far away.
    '''
    if len(pstcat) == 0: 
        return 0    #source might be in *.mag file, but has been sorted from the catalog before, so it is likely not a good source
    else:
        if len(pstcat) == 1:
            return (30-float(pstcat['MAG'][0])) * 10 #no neighbours, but make brighter stars better
        else:
            pstcat.sort('MAG')
            dist = pstcat.coords.distto((pstcat['XCENTER'][0],pstcat['YCENTER'][0]))
            return (float(pstcat['MAG'][1])-float(pstcat['MAG'][0])) * np.sqrt(dist[1])

    



def do_aperture_photometry(list_files):
    for datafile in list_files[:]:
        print datafile
        iraf.cd(os.path.dirname(datafile))
        # delete old aperture extractions
        for fl in glob.glob(datafile + '.coo.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.mag.*'):
            os.remove(fl)
        
        # extract sources.
        aperture_photometry(datafile,threshold=4.,wcs='logical',mode='small', telescope='Pairitel')
        found_sources_to_ds9(datafile + '.coo.1', 'apertureresult' )


def do_psf_photometry(list_files, satmag, psfstarfile=''):
    for datafile in list_files[:]:
        print datafile
        iraf.cd(os.path.dirname(datafile))
        # remove old files from previous runs (pst files = psf fitting stars, sub = subtracted images for psf fitting, als and arj = detected sources with psf fitting, nrj and nst = peak finding detections)
        # extract sources.
        for fl in glob.glob(datafile + '.pst.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.sub.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.als.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.arj.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.nst.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.nrj.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.psg.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.psf.*.fits'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.secondrun.mag.*'):
            os.remove(fl)
        
        psf_photometry_pairitel(datafile, satmag, psfstarfile)
        #found_sources_to_ds9(datafile + '.als.1')


def do_psf_photometry_with_coo(list_files, satmag, psfstarlist=''):
    for i in np.arange(0, len(list_files)):
        datafile = list_files[i]
        print datafile
        iraf.cd(os.path.dirname(datafile))
        # remove old files from previous runs (pst files = psf fitting stars, sub = subtracted images for psf fitting, als and arj = detected sources with psf fitting, nrj and nst = peak finding detections)
        # extract sources.
        for fl in glob.glob(datafile + '.pst.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.sub.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.als.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.arj.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.nst.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.nrj.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.psg.*'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.psf.*.fits'):
            os.remove(fl)
        
        for fl in glob.glob(datafile + '.secondrun.mag.*'):
            os.remove(fl)
        
        if (psfstarlist != ''):
            psf_photometry_pairitel_with_coo(datafile, satmag, psfstarfile=psfstarlist[i])
        
        else:
            psf_photometry_pairitel_with_coo(datafile, satmag)





def correct_coordinates(list_files, radius):
    # makes a copy of the original file and corrects the coordinates of the copies, using the 2MASS catalogue.
    for datafile in list_files[:]:
        print datafile
        print datafile.replace('.fits', '_wcs.fits')
        image = datafile.replace('_coadd_normed.fits', '_wcs.fits')
        shutil.copy(datafile,datafile.replace('_coadd_normed.fits', '_wcs.fits'))
        minmag = 15
        object=iraf.hselect(image+'[0]','OBJECT', 'yes', Stdout=1)[0]
        filter = iraf.hselect(image+'[0]','FILTER', 'yes', Stdout=1)[0]
        catalogue=os.path.join(os.path.dirname(image),"2mass_"+str(object)+'_'+str(filter)+'_'+str(minmag)+".cat")
        correct_wcs(datafile.replace('_coadd_normed.fits', '_wcs.fits'), clobber=True, radius=radius, minmag=minmag)
        twomass_to_ds9(catalogue)


def found_sources_to_ds9(sourcefile, suffix):
    #writes a downloaded 2mass catalog to a ds9 region file.
    data = asciitable.read(sourcefile, Reader=asciitable.Daophot)
    radius = 1.
    filename = os.path.dirname(sourcefile) + '/' + os.path.basename(sourcefile)[0:6] + '_' + suffix + '_found.reg'
    file=open(filename, 'w')
    file.write('# Region file format: DS9 version 4.1\n')
    file.write('# Filename ' + 'bla' + '\n')
    file.write('global color=blue dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    file.write('image\n')
    for i in range(len(data)):
	file.write('circle(' + str(data[i]['XCENTER']) + ',' + str(data[i]['YCENTER']) + ',' + str(radius) + '")' + '\n')
    
    file.close()


def twomass_to_ds9(catalogfile):
    #writes a downloaded 2mass catalog to a ds9 region file.
    data = asciitable.read(catalogfile)
    radius = 1.
    file=open(catalogfile.replace('cat','reg'), 'w')
    file.write('# Region file format: DS9 version 4.1\n')
    file.write('# Filename ' + 'bla' + '\n')
    file.write('global color=blue dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    file.write('fk5\n')
    for i in range(len(data)):
	file.write('circle(' + str(data[i]['col1']) + ',' + str(data[i]['col2']) + ',' + str(radius) + '")' + '\n')
    
    file.close()




def prepare_files(list_files, threshold, min_width, min_height, resultpath):
    # copies data to output directory, trim data files and mask files to contain only good exposure parts, normalize masks, normalize images by masks. Yields exposure-normalized, trimmed sky images.
    for datafile in list_files[:]:
	    #prepare directories and copy image to resultpath
	    filename=os.path.basename(datafile)
	    obsid=os.path.basename(os.path.dirname(datafile))
	    imagepath=os.path.join(resultpath,obsid)
	    image=os.path.join(imagepath,filename)
	    if not os.access(imagepath,os.F_OK): os.makedirs(imagepath)
	    print image
	    shutil.copy(datafile,imagepath)
	    os.chdir(imagepath)
	    if 'weight' in filename: # this only works because imlist is sorted and therefore the sky image was already copied to the right directory.
	        im_sky = image.replace('.weight','')
	        print "trimming..."
                apply_rectangle_mask(im_sky, image, threshold, min_width, min_height)
                print "normalizing mask..."
	        mask_norm(image.replace('.fits', '_trimmed.fits'))
	        print "applying mask..."
	        divide_by_mask(im_sky.replace('.fits', '_trimmed.fits'), image.replace('.fits', '_normed.fits'))
	        #photometry.correct_wcs(im_sky.replace('_trimmed.fits', '') + '_normed.fits')





def apply_rectangle_mask(image, mask, threshold, min_width, min_height):
    #finds rectangular middle part with good exposure and trims images to that part.
    hdulist = pyfits.open(image)
    imdata = hdulist[0].data
    #print len(imdata)
    hdulist.close()
    hdulist = pyfits.open(mask)
    immask = hdulist[0].data
    #print len(immask)
    hdulist.close()
    
    good_row = np.zeros(immask.shape[1])
    good_col = np.zeros(immask.shape[0])
    for i in np.arange(0, immask.shape[1]):
        good_col[i] = len(np.where(immask[i,:] > threshold*immask.max())[0])
    
    #print good_row
    for i in np.arange(0, immask.shape[0]):
        good_row[i] = len(np.where(immask[:,i] > threshold*immask.max())[0])
    
    # if sufficient exposure exists (i.e. not a bad frame), trim the image and the mask accordingly:
    if (len(np.where(good_row > min_width)[0]) > min_height) & (len(np.where(good_col > min_height)[0]) > min_width):
        good_width = 0.8 * max(good_row)
        good_height = 0.8 * max(good_col)
	ind_row = [np.where(good_row >= good_width)[0].min(), np.where(good_row >= good_width)[0].max()]
	ind_col = [np.where(good_col >= good_height)[0].min(), np.where(good_col >= good_height)[0].max()]
	
	print ind_row
	print ind_col
	
	outdata = deepcopy(imdata[ind_row[0]:ind_row[1], ind_col[0]:ind_col[1]])
	outmask = deepcopy(immask[ind_row[0]:ind_row[1], ind_col[0]:ind_col[1]])
	
	outimage = image.replace('.fits', '_trimmed.fits')
	outmask = mask.replace('.fits', '_trimmed.fits')
	
	iraf.imcopy(image + '[' + str(ind_row[0]+1) + ':' + str(ind_row[1]+1) + ',' + str(ind_col[0]+1) + ':' + str(ind_col[1]+1) + ']', outimage)
	iraf.imcopy(mask + '[' + str(ind_row[0]+1) + ':' + str(ind_row[1]+1) + ',' + str(ind_col[0]+1) + ':' + str(ind_col[1]+1) + ']', outmask)


def mask_norm(mask):
    # normalizes mask.
    # test if trimmed mask exists:
    if os.path.isfile(mask):
	hdulist = pyfits.open(mask)
	imdata = hdulist[0].data
	hdulist.close()
	outname = mask.replace('_trimmed.fits', '_normed.fits')
	print imdata.max()
	iraf.imarith(mask,"/",imdata.max(),outname)


def divide_by_mask(image, mask):
    # divides trimmed image by trimmed and normed mask.
    # test if trimmed image and trimmed+normed mask exist:
    if (os.path.isfile(mask)) & (os.path.isfile(image)):
	outname = image.replace('_trimmed.fits', '_normed.fits')
	iraf.imarith(image,"/",mask,outname)




def psf_star_checking(filename, satmag):
    # checks if there are stars in the psf fitting list which are brighter than satmag (choose satmag to exclude saturated stars), and discards those stars.
    print satmag
    f = open(filename, 'r')
    header = f.read()
    f.close()
    header = header.split('\n')
    indices = np.array([], int)
    for i in np.arange(0, len(header)-1): # -1 because last line is empty.
	if header[i][0] != '#':
	    if (header[i][29:34] == 'INDEF'):
		indices = np.append(indices,i)
	    else:
	        if (float(header[i][29:35]) < satmag):
	           indices = np.append(indices,i) 
    
    if len(indices) > 0: # tests if indices is empty.
        for i in indices:
	    header[i] = '!' # mark for deletion
	
	for i in np.arange(0, len(indices)):
	    header.remove('!')
	
	f = open(filename, 'w')
	for i in np.arange(0, len(header)):
	    f.write(header[i] + '\n')

	f.close()

def xycenter_to_ds9reg(daofile, outputname):
    
    regions = asciitable.read(daofile, Reader=asciitable.Daophot)
    radius = 1.
    file=open(outputname, 'w')
    file.write('# Region file format: DS9 version 4.1\n')
    file.write('# Filename ' + 'bla' + '\n')
    file.write('global color=blue dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    file.write('image\n')
    for i in range(len(regions)):
	file.write('circle(' + str(regions[i]['XCENTER']) + ',' + str(regions[i]['YCENTER']) + ',' + str(radius) + '")' + '\n')
    
    file.close()

def sort_by_sharpness(datalist, nbest):
    sharpnesses = np.zeros(len(datalist))
    n_found = np.zeros(len(datalist))
    for i in np.arange(0,len(datalist)):
        sharpnesses[i] = np.median(asciitable.read(datalist[i], Reader=asciitable.Daophot, fill_values=[('INDEF'), (np.nan)])['SHARPNESS'])
        n_found[i] = len(asciitable.read(datalist[i], Reader=asciitable.Daophot, fill_values=[('INDEF'), (np.nan)]))
    
    allstuff = zip(sharpnesses, datalist, n_found) # sort by first argument in zip.
    allstuff.sort()
    (sortsharp, sortnames, sortn) = zip(*allstuff)
    bestnames = sortnames[0:nbest]
    bests = sortsharp[0:nbest]
    
    allstuff = zip(n_found, datalist, sharpnesses)
    allstuff.sort()
    (sortn, sortnames, sortsharp) = zip(*allstuff)
    mostn = sortn[::-1][0:nbest]
    mostnames = sortnames[::-1][0:nbest]
    
    return (sharpnesses, bestnames, bests, mostnames, mostn)


