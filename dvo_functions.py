import numpy as np
from inspect import getfile, currentframe
from os import path
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

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

Last updated 13th Feb 2017

'''

# Class for querying Pan-STARRS Desktop Virtual Observatory (DVO)


class DvoFunctions(object):

    def __init__(self, skytable, photcode, PVurl, folderpath, storagepath):
        self.skytable = skytable
        self.photcode = photcode
        self.PVurl = PVurl
        self.folderpath = folderpath
        self.storagepath = storagepath

        self.rmin = skytable['R_MIN'].astype('float64')
        self.rmax = skytable['R_MAX'].astype('float64')
        self.dmin = skytable['D_MIN'].astype('float64')
        self.dmax = skytable['D_MAX'].astype('float64')
        self.depth = skytable['DEPTH']
        self.filelist = skytable['NAME']  # filenames
        self.index = skytable['INDEX']  # catid

        self.code = photcode['CODE']
        self.clam = photcode['C_LAM']
        self.Kcorr = photcode['K']

    # Finding great circle distance between given points and the list of points
    # given in equatorial coordinates.
    def GreatCircleDistance(self, ra1, dec1, ra_list, dec_list):
        # (ra1, dec1)is the reference point
        # (ra_list, dec_list)is the list of object positions

        ra1 = np.radians(ra1)
        dec1 = np.radians(dec1)
        ra2 = np.radians(ra_list)
        dec2 = np.radians(dec_list)

        # Great circle distance in radians
        dist = (np.arccos(np.sin(dec1) * np.sin(dec2) +
                          np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)))

        return np.degrees(dist) * 3600.

    # Converting from Galactic to Equatorial Coordinate System.
    def gal2equ(self, l, b, intype='degrees', outtype='degrees'):
        """
        Convert Galactic to Equatorial coordinates (J2000.0)
        (use at own risk)
        Source:
        - Book: "Practical astronomy with your calculator" (P. Duffett-Smith)
        - Wikipedia "Galactic coordinates"
        Parameters:
        [ l, b ]: Galactic longitude and latitude in decimal degrees
        Return:
        [ ra, dec ]: RA and DEC in decimal degrees
        """
        if (intype == 'degrees'):
            l = np.radians(l)
            b = np.radians(b)

        # North galactic pole (J2000) -- according to Wikipedia
        p_ra = np.radians(192.859508)
        p_dec = np.radians(27.128336)
        pa = np.radians(122.9319 - 90.)

        # North galactic pole (B1950)
        # p_ra = radians(192.25)
        # p_dec = radians(27.4)
        # pa = radians(123.0-90.0)

        ra = (np.arctan2((cos(b) * cos(l - pa)),
                         (sin(b) * cos(p_dec) - cos(b) * sin(p_dec) *
                          sin(l - pa))) + p_ra)
        dec = np.arcsin(cos(b) * cos(p_dec) * sin(l - pa) +
                        sin(b) * sin(p_dec))

        if (outtype == 'degrees'):
            ra = np.degrees(ra)
            dec = np.degrees(dec)

        return np.array([ra, dec])

    # Download the fits files required from the DVO distribution store with
    # shell script, and uncompress with the funpack from cfitsio
    def DownloadFile(self, filename, filetype):
        # Download file
        subprocess.call('wget ' + self.PVurl + filename + '.' + filetype +
                        ' -O ' + self.storagepath + '/' + filename[6:] + '.' +
                        filetype, shell=True)
        # Uncompress file, -F means writing over the compressed file
        subprocess.call(self.folderpath + '/cfitsio/funpack -F ' +
                        self.storagepath + '/' + filename[6:] + '.' +
                        filetype, shell=True)

    # Query Catalogue ID from a given position, default to be equatorial
    def PositionToCatid(self, x, y, system='equatorial'):

        if system == 'galactic':
            ra, dec = self.gal2equ(x, y)
        else:
            ra = x
            dec = y

        mask = ((ra >= self.rmin) & (ra < self.rmax) &
                (dec >= self.dmin) & (dec < self.dmax) &
                (self.depth == '\x04'))

        return self.index[mask][0]

    # Query raw file name from a given position
    def PositionToFilename(self, ra, dec):

        mask = ((ra >= self.rmin) & (ra < self.rmax) &
                (dec >= self.dmin) & (dec < self.dmax) &
                (self.depth == '\x04'))
        filename = filelist[mask][:14]

        if (filename[:10] != 's8230/pole')or (filename[:10] != 'n8230/pole'):
            filename = filename[:13]

        return filename

    # Query raw file name from a Catalogue ID
    def CatidToFilename(self, catid):
        # Take in CATID and return filepath, e.g. n0000/0001.01
        # Some filenames are not readable at \x04 level

        mask = (self.index == catid) & (self.depth == '\x04')
        filename = filelist[mask][0][:14]

        if (filename[:10] != 's8230/pole')or (filename[:10] != 'n8230/pole'):
            filename = filename[:13]

        return filename

    # Take in file name and file type, return fits table with both
    # header and data
    def LoadFitsFile(self, filename, filetype):
        if not path.isfile(self.storagepath + '/' + filename[6:] + '.' +
                           filetype):
            self.DownloadFile(filename, filetype)
        fitsfile = pyfits.open(self.storagepath + '/' + filename[6:] + '.' +
                               filetype)

        return fitsfile

    # Take in file name and file type, return FITS table header
    def LoadFitsHeader(self, filename, filetype):
        if not path.isfile(self.storagepath + '/' + filename[6:] + '.' +
                           filetype):
            self.DownloadFile(filename, filetype)

        fitsheader = pyfits.open(self.storagepath + '/' + filename[6:] + '.' +
                                 filetype)[1].header

        return fitsheader

    # Take in file name and file type, return the FITS table data
    # cpm and cpt files have conflicted headers, so the names have to be
    # changed before data can be loaded
    # nparray if False return the FITS table format
    def LoadFitsData(self, filename, filetype, nparray=True):
        if not path.isfile(self.storagepath + '/' + filename[6:] + '.' +
                           filetype):
            self.DownloadFile(filename, filetype)

        fitstable = pyfits.open(self.storagepath + '/' + filename[6:] + '.' +
                                filetype)[1]

        if filetype == 'cpm':
            fitstable.header['TTYPE44'] += '_1'
            fitstable.header['TTYPE60'] += '_2'
        if filetype == 'cpt':
            fitstable.header['TTYPE40'] += '_NUM_PARAMETER'
            fitstable.header['TTYPE41'] += '_UNASSIGNED'
        fitstable = fitstable.data

        if nparray:
            return np.array(fitstable)
        else:
            return fitstable

    # Take in position and radius (in arcsec), return list of Object IDs,
    # separations, file names and the cpt table
    def SearchAllObjidsInRadius(self, ra, dec, radius=5.):

        filename = self.PositionToFilename(ra, dec)[0]
        cpt = ftable(filename, 'cpt')
        dist = self.GreatCircleDistance(np.float64(ra), np.float64(dec),
                                        cpt['RA'].astype('float64'),
                                        cpt['DEC'].astype('float64'))
        mask = (dist < radius)

        return cpt['OBJ_ID'][mask], dist[mask], filename, cpt['CAT_ID'][0]

    # Take in CATID and OBJID, return average photometry
    # The sum in quadrature of 0.015 mag is to include systematics
    # Photometric system is available in psf, kron and aper
    def CatidAndObjidToAveragePhotometry(self, catid, objid, system='psf',
                                         allfilter=False):

        data = ftable(self.CatidToFilename(catid), 'cps')
        if allfilter:
            data = data[objid*9:(objid+1)*9]
        else:
            data = data[objid*9:objid*9+5]

        if system == 'kron':
            phot = data['MAG_KRON']
            photerr = np.sqrt(data['MAG_KRON_ERR']**2. + 0.015**2.)
        elif system == 'aper':
            phot = data['MAG_AP']
            photerr = np.sqrt(data['MAG_AP_ERR']**2. + 0.015**2.)
        else:
            phot = data['MAG']
            photerr = np.sqrt(data['MAG_ERR']**2. + 0.015**2.)

        return phot, photerr

    # Take in CATID and OBJID, return individual epoch and photometry
    # The sum in quadrature of 0.015 mag is to include systematics
    #
    # Photometric system: available in psf, kron and aper
    # tformat: available in unix, julian and reduced_julian
    # grouped: if True returns 5 sets of arrays for grizy respectively
    #
    # From Eugene Magnier:
    # mag:cat = MAG - 25.0 + C_LAM*0.001 + K*(AIRMASS - 1.0)
    # mag:rel = mag:cat - M_CAL
    # where MAG, AIRMASS, and M_CAL are from the cpm and C_LAM and K are
    # from the photcode table. K is the atmospheric absorption.
    def CatidAndObjidToEpochData(self, catid, objid, system='psf',
                                 mformat='calibrated', tformat='unix',
                                 grouped=False):

        data = ftable(self.CatidToFilename(catid), 'cpm')
        data = data[(data['OBJ_ID'] == objid)]
        data = data[(data['PHOTCODE'] >= 10000) & (data['PHOTCODE'] < 10500)]

        # astrometry
        ra = data['RA']
        dec = data['DEC']

        # photometry
        if system == 'psf':
            phot = data['MAG']
            photerr = np.sqrt(data['MAG_ERR']**2.0 + 0.015**2.0)
        elif system == 'kron':
            phot = data['MAG_KRON']
            photerr = np.sqrt(data['MAG_KRON_ERR']**2.0 + 0.015**2.0)
        elif system == 'aper':
            phot = data['MAG_AP']
            photerr = np.sqrt(data['MAG_AP_ERR']**2.0 + 0.015**2.0)

        photcodes = data['PHOTCODE']
        epoch = data['TIME']
        if mformat != 'instrumental':
            # calibrated magnitude
            airmass = data['AIRMASS']
            mcal = data['M_CAL']

            photcodes_pos = np.array([np.where(i == code)[0][0] for i
                                     in photcodes])
            clam_data = clam[photcodes_pos]
            Kcorr_data = Kcorr[photcodes_pos]

            # Instrumental zero-point correction
            phot = (phot - 25.0 + clam_data * 0.001 +
                    Kcorr_data * (airmass - 1.0) - mcal)

        # sky background flux, exposure time and db_flag
        sky = np.abs(data['SKY_FLUX'])
        exp = 10.**(data['M_TIME'] / 2.5)
        db_flag = data['DB_FLAGS']

        if tformat == 'julian':
            epoch = (data['TIME'] / 86400.0) + 2440587.5
        elif tformat == 'reduced_julian':
            epoch = (data['TIME'] / 86400.0) + 40587.5
        elif tformat == 'unix':
            pass
        else:
            print 'Unknown date format given, UNIX time is returned.'

        if grouped:
            mask = [(self.photcode >= 10000 + i * 100) &
                    (self.photcode < 10100 + i * 100) for i in range(5)]
            return [(ra[i], dec[i], phot[i], photerr[i], photcode[i], epoch[i],
                    sky[i], exp[i], db_flag[i]) for i in mask]
        else:
            return ra, dec, phot, photerr, photcode, epoch, sky, exp, db_flag

    # Take in CAT_ID and OBJ_ID and return the table for Astrometry
    def CatidAndObjidToAstrometry(self, catid, objid):

        data = ftable(self.CatidToFilename(catid), 'cpt')

        return data[data['OBJ_ID'] == objid]

    # Take in position and return average photometry of all objects within
    # radius (in arcsec)
    # Photometric system is available in psf, kron and aper
    def PositionSearchAveragePhotometry(self, ra, dec, dist=5.0):

        objid, dist, name, catid = self.SearchAllObjectsInRadius(ra, dec,
                                                                 radius=dist)
        data = []

        for i in objid:
            data.append(self.CatidAndObjidToAveragePhotometry(
                            catid, i, system='psf', allfilter=False))

        return np.array(data)

    # Take in position and return the astrometric table of all objects
    # within radius (in arcsec)
    def PositionSearchAstrometry(self, ra, dec, dist=5.0):

        objid, dist, name, catid = self.SearchAllObjectsInRadius(ra, dec,
                                                                 radius=dist)
        data = []

        for i in objid:
            data.append(self.CatidAndObjidToAstrometry(catid, i))

        return np.array(data)
