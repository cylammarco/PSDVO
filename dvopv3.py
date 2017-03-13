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
Version 3 (PV3) because the values in Photcodes.dat are slightly different.

Last updated 13th Feb 2017

'''

# Get the folder path of this script
pv3folderpath = path.dirname(path.abspath(getfile(currentframe())))

# Open the sky cell tessellation table and photcode table
pv3skytable = np.array(pyfits.getdata(folderpath + '/PV3/SkyTable.fits'))
pv3photcode = np.array(pyfits.getdata(folderpath + '/PV3/Photcodes.dat'))

# Set your storage folder, $HOME by default
pv3storagepath = path.expanduser("~/PV3")
# storagepath = "MANUAL PATH HERE"
if not path.exists(pv3storagepath):
    subprocess.call('mkdir ' + pv3storagepath, shell='True')

# Define the download path if file does not exist locallyin the storagepath
pv3url = "http://dvodist.ipp.ifa.hawaii.edu/3pi.pv3.20160422/"

dvopv3 = DvoFunctions(pv3skytable, pv3photcode, pv3url, pv3folderpath, pv3storagepath)
