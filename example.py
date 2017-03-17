from dvo_functions import *

# Initialising the DvoFunctions object
pv2 = DvoFunctions.set_version(2)
print 'Created an object pv2 initialised by DvoFunctions.set_version(2)'
pv3 = DvoFunctions.set_version(3)
print 'Created an object pv3 initialised by DvoFunctions.set_version(3)'

# an M7 dwarf
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=WISEA+J154045.67-510139.3
ra, dec = 162.0525, -11.3355
l, b = 260.9948, 41.2900

catid = 97606

#
# Query Catalogue ID from a given position, default to be equatorial
#
list_of_catid_1 = pv2.position_to_catid(ra, dec, system='equatorial')
list_of_catid_2 = pv2.position_to_catid(l, b, system='galactic')
if list_of_catid_1 != list_of_catid_2:
    print ('Searching with equatorial and galactic cooridinates of the same '
        'position return different list of catalogue ID!')

print 'pv2.position_to_catid(' + str(ra) + ', ' + str(dec) + ') returns:'
print list_of_catid_1
print ''
print ''

#
# Query raw file name from a given position, default to be equatorial
#
list_of_filename_1 = pv2.position_to_filename(ra, dec, system='equatorial')
list_of_filename_2 = pv2.position_to_filename(l, b, system='galactic')
if list_of_filename_1 != list_of_filename_2:
    print ('Searching with equatorial and galactic cooridinates of the same '
        'position return different list of file!')

print 'pv2.position_to_filename(' + str(ra) + ', ' + str(dec) + ') returns:'
print list_of_filename_1
print ''
print ''

#
# Query raw file name from a single Catalogue ID
#
filename = pv2.catid_to_filename(catid)
print 'pv2.catid_to_filename(' + str(catid) + ') returns:'
print filename
print ''
print ''

#
# Return the FITS file, in form of FITS object
#
cpt = pv2.load_fits_file(filename, 'cpt')
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cpt\') loaded astrometric file successfully.'
cpm = pv2.load_fits_file(filename, 'cpm')
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cpm\') loaded epoch file successfully.'
cps = pv2.load_fits_file(filename, 'cps')
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cps\') loaded average photometry file successfully.'
cpt = None
cpm = None
cps = None
print 'cpt, cpm and cps files dereferenced to free memory.'
print ''
print ''

#
# Return the FITS header
#
cpt_header = pv2.load_fits_header(filename, 'cpt')
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cpt\') loaded file header successfully.'
cpm_header = pv2.load_fits_header(filename, 'cpm')
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cpm\') loaded file header successfully.'
cps_header = pv2.load_fits_header(filename, 'cps')
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cps\') loaded file header successfully.'
print ''
print ''

#
# Return the FITS data
#
cpt_data = pv2.load_fits_data(filename, 'cpt', conflited_header=False,
                              nparray=True)
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cpt\') loaded file data successfully.'
cpm_data = pv2.load_fits_data(filename, 'cpm', conflited_header=False,
                              nparray=True)
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cpm\') loaded file data successfully.'
cps_data = pv2.load_fits_data(filename, 'cps', conflited_header=False,
                              nparray=True)
print 'pv2.load_fits_file(' + str(filename) +\
      ', \'cps\') loaded file data successfully.'
print ''
print ''

#
# Searching all ojbects within the radius from the given list of positions in
# equatorial coordinates
#        
list_of_catid, list_of_objid, list_of_angular_distance, list_of_filename =\
    pv2.search_objects_in_radius(ra, dec, radius=5.)
if len(list_of_objid) != len(list_of_angular_distance):
    print 'The size of the list of object IDs and angular separations are ' +\
          'different!'

print 'pv2.search_objects_in_radius(' + str(ra) + ', ' + str(dec) +\
      ', radius=5.) returns:'
print list_of_catid
print list_of_objid
print list_of_angular_distance
print list_of_filename
print ''
print ''
objid = list_of_objid[0]

#
# Return the average photometry and the associated uncertainties
#
phot = pv2.catid_objid_to_average_photometry(catid, objid, system='psf',
                                             allfilter=False)
print 'pv2.catid_objid_to_average_photometry(' + str(catid) + ', ' +\
      str(objid) + ', system=\'psf\', allfilter=False) returns:'
print phot
print ''
print ''

#
# Return some selected epoch-wise measurments
#
epoch = pv2.catid_objid_to_epoch_data(catid, objid, system='psf',
                                  mformat='calibrated', tformat='unix',
                                  grouped=False)
print 'pv2.catid_objid_to_epoch_data(' + str(catid) + ', ' + str(objid) +\
      ', system=\'psf\', mformat=\'calibrated\', tformat=\'unix\', ' +\
      'grouped=False) returns:'
print epoch
print ''
print ''

#
# Return the astrometric properties of an object with a given Catalogue and
# object ID
#
astro = pv2.catid_objid_to_astrometry(catid, objid)
print 'pv2.catid_objid_to_astrometry(' + str(catid) + ', ' + str(objid) +\
      ') returns:'
print astro
print ''
print ''

#
# Return the average photometry of objects within a given search radius from
# the given reference point (ra, dec). The default search radius is 5 arcsec
#
phot_searched = pv2.position_search_average_photometry(ra, dec, dist=5.0)
print 'pv2.position_search_average_photometry(' + str(ra) + ', ' + str(dec) +\
      ', radius=5.) returns:'
print phot_searched
print ''
print ''

#
# Return the average photometry of objects withint a given search radius from
# the given reference point (ra, dec). The default search radius is 5 arcsec
#
astro_searched = pv2.position_search_astrometry(ra, dec, dist=5.0)
print 'pv2.position_search_astrometry(' + str(ra) + ', ' + str(dec) +\
      ', radius=5.) returns:'
print astro_searched
print ''
print ''
