import subprocess
from os import path
import urllib2
from matplotlib.pyplot import *
from dvo_functions import DvoFunctions

class postage_stamp(DvoFunctions):

    def __init__(self, path_location=path.expanduser("~") + "/postage_stamps"):
        self.savepath = path_location
        if not path.isdir(self.savepath):
            subprocess.call('mkdir ' + self.savepath, shell='True')

    def _query_postage_stamp(self, ra, dec, box_size=90, image_dimension=1024,
                             system='equatorial'):
        '''
        Query for postage stamp by wrapping through the PS postage stamp
        service at http://jflkaslfkalld
        which currently only allows single object query

        Input:
        ra, dec - (RA, Dec) in J2000.0
        box_size - the size of the image in arcseconds
        image_dimension - the output image dimension in number of pixels
        system - "equatorial" or "galactic", case sensitive
        formatted - if True, return a 2 x 3 grid of images, with a y/i/g
                    flase colour image on the top left
                    *************
                    * y/i/g * g *
                    *   r   * i *
                    *   z   * y *
                    *************
        saved - if True, save the figure to disk, default path is the home
                directory
        image_format - "jpg", "fits" or "both"
        # the following settings may be useful if user is using python/ipython/
        # python notebook interactively. When using this function in script,
        # stick to the default options to improve runtime
        display - if True, open image with matplotlib.pyplot.show()
        interactive - if True, matplotlib.pyplot.ion() is called. Need to call
                      matplotlib.pyplot.ioff() manually to turn off interative mode

        Output:
        image displayed on screen or saved to the savepath or neither. Depending on
        your input setting.

        '''
        # Conversion for display the PS web interface, don't change the value here
        image_size = box_size * 4

        # the url for returning false colour and total stacked images in each
        # of the grizy filter
        url = "http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=" +\
              str(ra) + "%2C+" + str(dec) + "&filter=color&filter=g&" +\
              "filter=r&filter=i&filter=z&filter=y&filetypes=stack&" +\
              "auxiliary=data&size=" + str(image_size) + "&output_size=" +\
              str(box_size) + "&verbose=0&autoscale=99.500000&catlist="

        # loading the html script of the given url link
        page = urllib2.urlopen(url)

        # load the html script as text
        data = page.read()

        # splitting the text using " as a delimiter
        data = data.split('"')
        html_links = []

        for i in range(len(data)):
            text = data[i]
            # extract the html link for the images
            if ((text[:42] == '//ps1images.stsci.edu/cgi-bin/fitscut.cgi?') &
                (text[-(len(str(box_size))):] == str(box_size))):
                html_links.append('http:' + text)
        return html_links

    def ps_stacked_image(self, x, y, box_size=90, image_dimension=1024,
                         system='equatorial', formatted=True, saved=True,
                         image_format='jpg', display=False, interactive=False):
        '''
        Query for postage stamp by wrapping through the PS postage stamp
        service at http://jflkaslfkalld
        which currently only allows single object query

        Input:
        x, y - either (RA, Dec) or (l, b) depending on the chosen system
        box_size - the output image dimension in number of pixels
        image_dimension - the size of the output image in pixels
        system - "equatorial" or "galactic", case sensitive
        formatted - if True, return a 2 x 3 grid of images, with a y/i/g
                    flase colour image on the top left
                    *************
                    * y/i/g * g *
                    *   r   * i *
                    *   z   * y *
                    *************
        saved - if True, save the figure to disk, default path is the home
                directory
        image_format - "jpg", "fits" or "both"
        # the following settings may be useful if user is using python/ipython/
        # python notebook interactively. When using this function in script,
        # stick to the default options to improve runtime
        display - if True, open image with matplotlib.pyplot.show()
        interactive - if True, matplotlib.pyplot.ion() is called. Need to call
                      matplotlib.pyplot.ioff() manually to turn off interative mode

        Output:
        image displayed on screen or saved to the savepath or neither. Depending on
        your input setting.

        '''
        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        elif system == 'equatorial':
            ra = x
            dec = y
        else:
            print "Unknown coordinate system. Please choose between " +\
                  "equatorial galactic (case sensitive)."

        # get the image links
        html_links = self._query_postage_stamp(ra, dec, box_size, image_dimension,
                                               system)

        if formatted:
            # population the image with the 6 images
            n = -1
            for j in enumerate(html_links):
                if (j[0]%6 == 0):
                    fig, ax = subplots(3, 2, figsize=(7.5,12))
                    suptitle(str(ra) + ", " + str(dec) )
                    tight_layout()
                    n += 1
                val = (j[0]-n*6)
                x = val % 2
                y = val / 2
                f = imread(urllib2.urlopen(j[1]), format='jpg')
                if j[0] == 0:
                    ax[y,x].imshow(f)
                    ax[y,x].get_xaxis().set_ticks([])
                    ax[y,x].get_yaxis().set_ticks([])
                else:
                    ax[y,x].imshow(f, cmap='Greys')
                    ax[y,x].get_xaxis().set_ticks([])
                    ax[y,x].get_yaxis().set_ticks([])
                ax[0,0].set_xlabel('y-i-g False Colour')
                ax[0,1].set_xlabel('g')
                ax[1,0].set_xlabel('r')
                ax[1,1].set_xlabel('i')
                ax[2,0].set_xlabel('z')
                ax[2,1].set_xlabel('y')
            if saved:
                savefig(self.savepath + '/stack_' + ("%f6" % ra) + '_' +\
                        ("%f6" % dec) + '.png')
        else:
            for j in enumerate(html_links):
                f = imread(urllib2.urlopen(j[1]), format='jpg')
                fig, ax = figure(j[0])
                if j[0] == 0:
                    ax.imshow(f)
                else:
                    ax.imshow(f, cmap='Greys')
        if interactive:
            ion()
        if display:
            show()
        fig.clf()
        fig = None
        close()


    def query_eopch_postage_stamp(self, x, y, box_size=90, image_size=1024,
                                  system='equatorial', formatted=True):
        if system == 'galactic':
            ra, dec = self._galactic_to_equatorial(x, y)
        elif system == 'equatorial':
            ra = x
            dec = y
        else:
            print "Unknown coordinate system. Please choose between " +\
                  "equatorial galactic (case sensitive)."

        url = "http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=" +\
              str(ra) + "%2C+" + str(dec) + "&filter=color&filter=g&" +\
              "filter=r&filter=i&filter=z&filter=y&filetypes=stack&" +\
              "filetypes=warp&auxiliary=data&size=" + str(image_size) +\
              "&output_size=" + str(output_size) +\
              "&verbose=0&autoscale=99.500000&catlist="

