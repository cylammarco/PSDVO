# PSDVO
To load files from Pan-STARRS DVO (Desktop Virtual Observatory)

Require IfA (or other participating institute) IP address (so use ROE VPN or use the linux workstation).

pyfits is not updated with the newest cfitsio yet, so PV3 data require manual uncompression with the funpack from the cfitsio folder. It is included in the PV3 version of the code.

The PV3 download scripts only works on Linux machine at the moment (if you already have the uncompressed files, it's fine), trying to figure out how to use different formats of cfitsio based on the machine type. If you want to load with a Mac, please manually replace the cfitsio folder with the Mac Version which can be downloaded from http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html . 

Files will be downloaded to $HOME/PV2 or $HOME/PV3, those folders will be created automatically if they don't exist. However, the downloaded files WILL **NOT** BE DELETED AUTOMATICALLY, make sure you check the disk space regularly, because the complete PV2 catalogue is ~15TB and PV3 is ~25TB after uncompression.

The script is not tested with cpn, cpx and cpy files because I have no idea what they are. But in theory ffile can load them properly.


PV2
http://dvodist.ipp.ifa.hawaii.edu/3pi.pv2.20141215

PV3
http://dvodist.ipp.ifa.hawaii.edu/3pi.pv3.20160422
