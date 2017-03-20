import numpy as np
import subprocess
import collections
from inspect import getfile, currentframe
from os import path
from sys import platform
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

Last updated 15th March 2017

'''

# Class for querying Pan-STARRS Desktop Virtual Observatory (DVO)


class DvoFunctions(object):

    def __init__(self, version):
        '''

        Initialize the parameters of a DvoFunctions object.

        Input:
        SkyTable - define the skycells
        Photcode - define the photometric system
        PVurl - the url of the DVO distribtuion store
        folderpath - path of the folder storing SkyTable and Photcode
        storagepath - path of the folder storing data

        Parameters:
        rmin, rmax - range of the RA bounding box of a cell
        dmin, dmax - range of the Dec bounding box of a cell
        depth - level of skycell subdivision (from 0 to 4)
        filelist - the list of file names of the cells
        index - the Catalogue ID for the cell
        code - the ID for the CCD
        clam - instrumental zero point for the CCD
        Kcorr - the atmospheric absorption (multiply by airmass to get the
                total absorption)

        '''
        # Define the download path if file does not exist locally in the
        # storagepath
        if version == 2:
            PVurl = "http://dvodist.ipp.ifa.hawaii.edu/3pi.pv2.20141215/"
        elif version == 3:
            PVurl = "http://dvodist.ipp.ifa.hawaii.edu/3pi.pv3.20160422/"
        else:
            print "Required version not available"
            return None

        # Get the folder path of this script
        folderpath = path.dirname(path.abspath(getfile(currentframe())))

        # Open the sky cell tessellation table and photcode table
        skytable = np.array(pyfits.getdata(folderpath + '/PV' +
                              ("%d" % version) + '/SkyTable.fits'))
        photcode = np.array(pyfits.getdata(folderpath + '/PV' +
                              ("%d" % version) + '/Photcodes.dat'))

        # Set your storage folder, $HOME by default
        storagepath = path.expanduser('~/PV%d' % version)
        # storagepath = "MANUAL PATH HERE"
        if not path.exists(storagepath):
            subprocess.call('mkdir ' + storagepath, shell='True')

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

    def _great_circle_distance(self, ra1, dec1, ra_list, dec_list):
        '''

        Finding the great circle distance between given points and the list of
        points given in equatorial coordinates.

        Input:
        ra1, dec1 - the reference point
        ra_list, dec_list - is the list of object positions

        Output:
        list of distances between the reference point and the list of objects
        in units of arcseconds.

        '''

        ra1 = np.radians(ra1)
        dec1 = np.radians(dec1)
        ra2 = np.radians(ra_list)
        dec2 = np.radians(dec_list)

        # Great circle distance in radians
        dist = (np.arccos(np.sin(dec1) * np.sin(dec2) +
                          np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)))

        return np.degrees(dist) * 3600.

    def _galactic_to_equatorial(self, l, b, intype='degrees',
                                outtype='degrees'):
        '''

        Convert Galactic to Equatorial coordinates (J2000.0)

        Source:
        - Book: "Practical astronomy with your calculator" (P. Duffett-Smith)
        - Wikipedia "Galactic coordinates"

        Input:
        l, b - Galactic longitude and latitude in decimal degrees
        intype - if degrees, convert to radians; else do nothing
        outtype - if degrees, convert back to degrees; else do nothing

        Output:
        ra, dec - RA and Dec in decimal degrees

        '''

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

        ra = (np.arctan2((np.cos(b) * np.cos(l - pa)),
                         (np.sin(b) * np.cos(p_dec) - np.cos(b) * np.sin(p_dec) *
                          np.sin(l - pa))) + p_ra)
        dec = np.arcsin(np.cos(b) * np.cos(p_dec) * np.sin(l - pa) +
                        np.sin(b) * np.sin(p_dec))

        if (outtype == 'degrees'):
            ra = np.degrees(ra)
            dec = np.degrees(dec)

        return np.array([ra, dec])

    def _download_file(self, filename, filetype):
        '''

        Download the FITS files required from the DVO distribution store with
        shell script and uncompress with the funpack from cfitsio

        Input:
        filename - filename (xxxx.yy where xxxx are 0-9 and yy = [00-15])
        filetype - 6 file types available at the store:
                   (1) cpm - epoch measurements
                   (2) cpt - average astrometry
                   (3) cps - average photometry
                   (4) cpn - missing objects
                   (5) cpx - weak lensing
                   (6) cpy - weak lensing

        '''

        # Download file
        subprocess.call('wget ' + self.PVurl + filename + '.' + filetype +
                        ' -O ' + self.storagepath + '/' + filename[6:] + '.' +
                        filetype, shell=True)

        # Uncompress file, -F means writing over the compressed file
        if platform == "darwin":
            subprocess.call(self.folderpath + '/cfitsio_mac/funpack -F ' +
                        self.storagepath + '/' + filename[6:] + '.' +
                        filetype, shell=True)
        else:
            subprocess.call(self.folderpath + '/cfitsio/funpack -F ' +
                        self.storagepath + '/' + filename[6:] + '.' +
                        filetype, shell=True)

    def position_to_catid(self, x, y, system='equatorial'):
        '''

        Query Catalogue ID from a given position, default to be equatorial

        Input:
        x, y - 2D coordinates in decimal degrees
        system - if galactic, convert to equatorial; else do nothing

        Output:
        array of Catalogue ID

        '''

        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        else:
            ra = x
            dec = y

        # Finding the skycell that contained the object
        mask = ((ra >= self.rmin) & (ra < self.rmax) &
                (dec >= self.dmin) & (dec < self.dmax) &
                (self.depth == '\x04'))

        return self.index[mask][0]

    def position_to_filename(self, x, y, system='equatorial'):
        '''

        Query raw file name from a given position, default to be equatorial

        Input:
        x, y - 2D coordinates in decimal degrees
        system - if galactic, convert to equatorial; else do nothing

        Output:
        file path (e.g. n0000/0001.01)
        N.B. Some names are not readable at \x04 level in some older releases

        '''

        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        else:
            ra = x
            dec = y

        mask = ((ra >= self.rmin) & (ra < self.rmax) &
                (dec >= self.dmin) & (dec < self.dmax) &
                (self.depth == '\x04'))
        filename = self.filelist[mask][:14]

        if (filename[:10] != 's8230/pole')or (filename[:10] != 'n8230/pole'):
            filename = filename[:13]

        return filename

    def catid_to_filename(self, catid):
        '''

        Query raw file name from a single Catalogue ID

        Input:
        catid - Catalogue ID

        Output:
        file path (e.g. n0000/0001.01)
        N.B. Some names are not readable at \x04 level in some older releases

        '''

        mask = (self.index == catid) & (self.depth == '\x04')
        filename = self.filelist[mask][0][:14]

        if (filename[:10] != 's8230/pole')or (filename[:10] != 'n8230/pole'):
            filename = filename[:13]

        return filename

    def load_fits_file(self, filename, filetype):
        '''

        Return the FITS file, in form of FITS object

        Input:
        filename - file path (e.g. n0000/0001.01)
        filetype - the 6 file types (m,t,s,n,x,y)

        Output:
        FITS object

        '''

        if not path.isfile(self.storagepath + '/' + filename[6:] + '.' +
                           filetype):
            self._download_file(filename, filetype)
        fitsfile = pyfits.open(self.storagepath + '/' + filename[6:] + '.' +
                               filetype)

        return fitsfile

    def load_fits_header(self, filename, filetype):
        '''

        Return the FITS header

        Input:
        filename - file path (e.g. n0000/0001.01)
        filetype - the 6 file types (m,t,s,n,x,y)

        Output:
        FITS header

        '''

        if not path.isfile(self.storagepath + '/' + filename[6:] + '.' +
                           filetype):
            self._download_file(filename, filetype)

        fitsheader = pyfits.open(self.storagepath + '/' + filename[6:] + '.' +
                                 filetype)[1].header

        return fitsheader

    def load_fits_data(self, filename, filetype, conflited_header=False,
                       nparray=True):
        '''

        Return the FITS table. Some older versions have header issue, namely
        PV2 downloaded before ~late 2016

        Input:
        filename - file path (e.g. n0000/0001.01)
        filetype - the 6 file types (m,t,s,n,x,y)

        Output:
        FITS table, default to be converted to Numpy Array for efficiency gain,
        set to False to retrieve as FITS object table

        '''

        if not path.isfile(self.storagepath + '/' + filename[6:] + '.' +
                           filetype):
            self._download_file(filename, filetype)

        fitstable = pyfits.open(self.storagepath + '/' + filename[6:] + '.' +
                                filetype)[1]

        # Fix header issue
        if conflited_header:
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

    def search_objects_in_radius(self, x, y, radius=5., system='equatorial'):
        '''

        Searching all ojbects within the radius from the given list of
        positions in equatorial coordinates.

        Input:
        x, y - 2D coordinates in decimal degrees
        dist - list of separations between reference point and objects within
               the given search radius, in unit of arcseconds
        system - if galactic, convert to equatorial; else do nothing

        Output:
        Catalogue ID
        list of Object IDs
        dist - list of separations between reference point and objects within
               the given search radius, in unit of arcseconds
        filename - file path (e.g. n0000/0001.01)

        '''

        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        else:
            ra = x
            dec = y

        filename = self.position_to_filename(ra, dec)[0]
        cpt = self.load_fits_data(filename, 'cpt')
        dist = self._great_circle_distance(np.float64(ra), np.float64(dec),
                                           cpt['RA'].astype('float64'),
                                           cpt['DEC'].astype('float64'))
        mask = (dist < radius)

        return cpt['CAT_ID'][0], cpt['OBJ_ID'][mask], dist[mask], filename

    def catid_objid_to_average_photometry(self, catid, objid, system='psf',
                                          allfilter=False):
        '''

        Return the average photometry and the associated uncertainties.

        Input:
        catid - Catalogue ID
        objid - Object ID
        system - the native psf, kron and 3"-aperture magnitude (not a
                 conversion between them)
        allfilter - return clear(sum grizy) and w(sum ugr) filter as well,
                    which are usually empty, mostly only available when the
                    position concides the asteroid and/or planet survey.
                    Default to return 5 filters, the grizy, if allfilter is
                    True, return 9 filters.

        Output:
        phot - numpy array of the magnitudes in the requested fitlers (5 or 9),
               the sum in quadrature of 0.015 mag is to include systematic
               uncertainty of the response system
        photerr - numpy array of magnitude uncertainties in the requested
                  filters

        '''

        data = self.load_fits_data(self.catid_to_filename(catid), 'cps')
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

    def catid_objid_to_epoch_data(self, catid, objid, system='psf',
                                  mformat='calibrated', tformat='unix',
                                  grouped=False):
        '''

        Return some selected epoch-wise measurments

        Input:
        catid - Catalogue ID
        objid - Object ID
        system - the native psf, kron and 3"-aperture magnitude (not a
                 conversion between them)
        mformat - default to be 'calibrated' such that intrumental zero-point
                  and the atmospheric absorption are corrected, if
                  'instrumental' the raw instrumetal magnitudes will be
                  returned
        tformat - default in unix, also available in julian and reduced_julian
        grouped - if True returns 5 sets of arrays for grizy respectively, else
                  1 single array

        Output:
        ra, dec - position in degrees
        phot - numpy array of the magnitudes in the requested fitlers (5 or 9),
               the sum in quadrature of 0.015 mag is to include systematic
               uncertainty of the response system
        photerr - numpy array of magnitude uncertainties in the requested
                  filters
        photcode - the OTA numbering system for different available filters,
                   the last two digits refer to the X and Y coordinate of the
                   CCD at the detector: 01-06, 10-17, 20-27, ...60-67, 71-76
                   g: 10000-10080
                   r: 10100-10180
                   i: 10200-10280
                   z: 10300-10380
                   y: 10400-10480
        epoch - the epoch of the measurements in the tformat requested
        sky - the measured local sky brightness
        exp - the exposure times
        db_flag - the database flags

        N.B.
        From Eugene Magnier:
        mag:cat = MAG - 25.0 + C_LAM*0.001 + K*(AIRMASS - 1.0)
        mag:rel = mag:cat - M_CAL
        where MAG, AIRMASS, and M_CAL are from the cpm and C_LAM and K are
        from the photcode table. K is the atmospheric absorption.

        '''

        data = self.load_fits_data(self.catid_to_filename(catid), 'cpm')
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
            airmass = data['AIRMASS']
            # magnitude calibration
            mcal = data['M_CAL']

            photcodes_pos = np.array([np.where(i == self.code)[0][0] for i
                                     in photcodes])
            clam_data = self.clam[photcodes_pos]
            Kcorr_data = self.Kcorr[photcodes_pos]

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
            return [(ra[i], dec[i], phot[i], photerr[i], photcodes[i],
                    epoch[i], sky[i], exp[i], db_flag[i]) for i in mask]
        else:
            return ra, dec, phot, photerr, photcodes, epoch, sky, exp, db_flag

    def catid_objid_to_astrometry(self, catid, objid):
        '''

        Return the astrometric properties of an object with a given Catalogue
        and object ID.

        Input:
        catid - Catalogue ID
        objid - Object ID

        Output:
        The astrometry of the object requested

        '''

        data = self.load_fits_data(self.catid_to_filename(catid), 'cpt')

        return data[data['OBJ_ID'] == objid]

    def position_search_average_photometry(self, x, y, dist=5.0,
                                           system='equatorial'):
        '''

        Return the average photometry of objects within a given search radius
        from the given reference point (ra, dec). The default search radius is
        5 arcseconds.

        Input:
        x, y - 2D coordinates in decimal degrees
        dist - list of separations between reference point and objects within
               the given search radius, in unit of arcseconds
        system - if galactic, convert to equatorial; else do nothing

        Output:
        The photometry of the objects searched within the radius

        '''

        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        else:
            ra = x
            dec = y

        catid, objid, dist, name = self.search_objects_in_radius(ra, dec,
                                                                 radius=dist)

        if isinstance(objid, collections.Iterable):
            data = []
            for i in np.array(objid):
                data.append(self.catid_objid_to_average_photometry(
                            catid, i, system='psf', allfilter=False))
        else:
            data = self.catid_objid_to_average_photometry(
                       catid, objid, system='psf', allfilter=False)

        return np.array(data)

    def position_search_astrometry(self, x, y, dist=5.0, system='equatorial'):
        '''

        Return the average photometry of objects withint a given search radius
        from the given reference point (ra, dec). The default search radius is
        5 arcseconds.

        Input:
        x, y - 2D coordinates in decimal degrees
        dist - list of separations between reference point and objects within
               the given search radius, in unit of arcseconds
        system - if galactic, convert to equatorial; else do nothing

        Output:
        The astrometry of the objects searched within the radius

        '''

        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        else:
            ra = x
            dec = y

        catid, objid, dist, name = self.search_objects_in_radius(ra, dec,
                                                                 radius=dist)
        if isinstance(objid, collections.Iterable):
            data = []
            for i in np.array(objid):
                data.append(self.catid_objid_to_astrometry(catid, i))
        else:
            data = self.catid_objid_to_astrometry(catid, objid)

        return np.array(data)
