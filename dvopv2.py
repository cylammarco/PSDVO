import numpy as np
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
import subprocess
from inspect import getfile, currentframe
from os import path
from dvo_functions import *

'''

This script is for query of small number of objects. The searching algorithm
is not optimised for large scale query.

File organisation:
cpt - astrometry
cps - average photometry
cpm - epoch photometry and astrometry
cpn - known missing objects (?)
cpx, cpy - lensing shearing
SkyTable - define the skycells
Photcode - define the photometric system

Specifically for Pan-STARRS Desktop Virtual Observatory (DVO) Processing
Version 2 (PV2) because the values in Photcodes.dat are slightly different.

Last updated 13th Feb 2017

'''

# Get the folder path of this script
folderpath = path.dirname(path.abspath(getfile(currentframe())))

# Open the sky cell tessellation table and photcode table
skytable = np.array(pyfits.getdata(folderpath + '/PV2/SkyTable.fits'))
photcode = np.array(pyfits.getdata(folderpath + '/PV2/Photcodes.dat'))

# Set your storage folder, $HOME by default
storagepath = path.expanduser("~/PV2")
# storagepath = "MANUAL PATH HERE"
if not path.exists(storagepath):
    subprocess.call('mkdir ' + storagepath, shell='True')

# Define the download path if file does not exist locallyin the storagepath
PVurl = "http://dvodist.ipp.ifa.hawaii.edu/3pi.pv2.20141215/"

dvopv2 = DvoFunctions(skytable, photcode, PVurl, folderpath, storagepath)
