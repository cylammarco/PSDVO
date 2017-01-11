import numpy as np
import pyfits

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


Specifically for Pan-STARRS Desktop Virtual Observatory (DVO)Processing
Version 2 (PV2)because of the following problem:
(1)cpm and cpt have conflicted headers which has to be fixed before
    reading data.
(2)Photcode is bound to the respective PV

Last updated 19th Oct 2016

'''

folderpath = '/disk8/mlam/lap.pv2.20141215/'

skytable = np.array(pyfits.getdata(folderpath+'SkyTable.fits'))
photcode = np.array(pyfits.getdata(folderpath+'Photcodes.dat'))

# RA, Dec need to be at double quadruple precision
rmin = skytable['R_MIN'].astype('float64')
rmax = skytable['R_MAX'].astype('float64')
dmin = skytable['D_MIN'].astype('float64')
dmax = skytable['D_MAX'].astype('float64')
depth = skytable['DEPTH']
filelist = skytable['NAME']  # filenames
index = skytable['INDEX']  # catid

code = photcode['CODE']
clam = photcode['C_LAM']
Kcorr = photcode['K']


def _greatcirdist(ra1, dec1, ra_list, dec_list):
    # (ra1, dec1)is the reference point
    # (ra_list, dec_list)is the list of object positions

    ra1 = np.radians(ra1)
    dec1 = np.radians(dec1)
    ra2 = np.radians(ra_list)
    dec2 = np.radians(dec_list)

    # Great circle distance in radians
    dist = (np.arccos(np.sin(dec1) * np.sin(dec2) +
                      np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)))

    return np.degrees(dist)*3600.


def pos2catid(ra, dec):
    # Take in position and return CAT_ID

    mask = ((ra >= rmin) & (ra < rmax) & (dec >= dmin) & (dec < dmax) &
            (depth == '\x04'))

    return index[mask][0]


def pos2filename(ra, dec):
    # Take in position and return filename

    mask = ((ra >= rmin) & (ra < rmax) & (dec >= dmin) & (dec < dmax) &
            (depth == '\x04'))
    filename = filelist[mask][:14]

    if (filename[:10] != 's8230/pole')or (filename[:10] != 'n8230/pole'):
        filename = filename[:13]

    return filename


def catid2name(catid):
    # Take in CATID and return filepath, e.g. n0000/0001.01
    # Some filenames are not readable at \x04 level

    mask = (index == catid) & (depth == '\x04')
    filename = filelist[mask][0][:14]

    if (filename[:10] != 's8230/pole')or (filename[:10] != 'n8230/pole'):
        filename = filename[:13]

    return filename


def ffile(filename, filetype):
    # Return the original FITS file

    fitsfile = pyfits.open(folderpath + filename + '.' + filetype)

    return fitsfile


def fheader(filename, filetype):
    # Return the FITS table header

    fitsheader = pyfits.open(folderpath + filename + '.' + filetype)[1].header

    return fitsheader


def ftable(filename, filetype, nparray=True):
    # Return the FITS table data
    # cpm and cpt files have conflicted headers, so the names have to be
    # changed before data can be loaded
    # nparray if False return the FITS table format

    fitstable = pyfits.open(folderpath + filename + '.' + filetype)[1]

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


def search_radius_objid(ra, dec, radius=5.):
    # Take in position and return all objects within radius (in arcsec)

    filename = pos2filename(ra, dec)[0]
    cpt = ftable(filename, 'cpt')
    dist = _greatcirdist(np.float64(ra), np.float64(dec),
                         cpt['RA'].astype('float64'),
                         cpt['DEC'].astype('float64'))
    mask = (dist < radius)

    return cpt['OBJ_ID'][mask], dist[mask], filename, cpt['CAT_ID'][0]


def id2avgphot(catid, objid, system='psf', allfilter=False):
    # Take in CATID and OBJID, return average photometry
    # The sum in quadrature of 0.015 mag is to include systematics
    # Photometric system is available in psf, kron and aper

    data = ftable(catid2name(catid), 'cps')
    if allfilter:
        data = data[objid*9:(objid+1)*9]
    else:
        data = data[objid*9:objid*9+5]

    if system == 'kron':
        phot = data['MAG_KRON']
        photerr = np.sqrt(data['MAG_KRON_ERR']**2.0 + 0.015**2.0)
    elif system == 'aper':
        phot = data['MAG_AP']
        photerr = np.sqrt(data['MAG_AP_ERR']**2.0 + 0.015**2.0)
    else:
        phot = data['MAG']
        photerr = np.sqrt(data['MAG_ERR']**2.0 + 0.015**2.0)

    return phot, photerr


def id2epoch_data(catid, objid, system='psf', mformat='calibrated',
    tformat='unix', grouped=False):

    # Take in CATID and OBJID, return individual epoch and photometry
    # The sum in quadrature of 0.015 mag is to include systematics
    # Photometric system is available in psf, kron and aper
    # tformat is available in unix, julian and reduced_julian
    # grouped if True returns 5 sets of arrays for grizy respectively
    # From Eugene Magnier:
    # mag:cat = MAG - 25.0 + C_LAM*0.001 + K*(AIRMASS - 1.0)
    # mag:rel = mag:cat - M_CAL
    # where MAG, AIRMASS, and M_CAL are from the cpm and C_LAM and K are
    # from the photcode table.

    data = ftable(catid2name(catid), 'cpm')
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

    if mformat != 'instrumental':
        # calibrated magnitude
        airmass = data['AIRMASS']
        mcal = data['M_CAL']
        photcode = data['PHOTCODE']
        epoch = data['TIME']

        photcode_pos = np.array([np.where(i == code)[0][0] for i in photcode])
        clam_data = clam[photcode_pos]
        Kcorr_data = Kcorr[photcode_pos]

        # Instrumental zero-point correction
        phot = (phot - 25.0 + clam_data * 0.001 +
                Kcorr_data * (airmass - 1.0) - mcal)

    # sky background flux
    sky = np.abs(data['SKY_FLUX'])
    exp = 10.**(data['M_TIME'] / 2.5)

    if tformat == 'julian':
        epoch = (data['TIME'] / 86400.0) + 2440587.5
    elif tformat == 'reduced_julian':
        epoch = (data['TIME'] / 86400.0) + 40587.5
    elif tformat == 'unix':
        pass
    else:
        print 'Unknown date format given, UNIX time is returned by default.'

    if grouped:
        mask = [(photcode >= 10000 + i * 100) & (photcode < 10100 + i * 100)
                for i in range(5)]
        return [(ra[i], dec[i], phot[i], photerr[i], photcode[i], epoch[i],
                sky[i], exp[i]) for i in mask]
    else:
        return phot, photerr, photcode, epoch, sky, exp


def id2astro(catid, objid):
    # Take in CAT_ID and OBJ_ID and return the table for Astrometry

    data = ftable(catid2name(catid), 'cpt')

    return data[data['OBJ_ID'] == objid]


def pos2avgphot(ra, dec, dist=5.0):
    # Take in position and return average photometry of all objects within
    # radius (in arcsec)
    # Photometric system is available in psf, kron and aper

    objid, dist, name, catid = search_radius_objid(ra, dec, radius=dist)
    data = []

    for i in objid:
        data.append(id2avgphot(catid, i, system='psf', allfilter=False))

    return np.array(data)


def pos2astro(ra, dec, dist=5.0):
    # Take in position and return the astrometric table of all objects within
    # radius (in arcsec)

    objid, dist, name, catid = search_radius_objid(ra, dec, radius=dist)
    data = []

    for i in objid:
        data.append(id2astro(catid, i))

    return np.array(data)
