#!/usr/bin/env python
import argparse
from astropy.io import fits
from astropy.visualization import LogStretch, LinearStretch, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
import datetime as dt
import glob
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os



parser = argparse.ArgumentParser(description="Create an animation of slices"\
                                             " from a list of fits images")

parser.add_argument("-path",
                    help="Path to a directory containing "\
                         "fits images (must include trailing /)",
                    type=str,
                    default='../test_data/')

parser.add_argument("-suffix",
                    help="Suffix of fits files (e.g. *flt.fits, j*flt.fits)",
                    type=str,
                    default='*flt.fits')

parser.add_argument("-xcenter",
                    help="x value for center of box slice",
                    type=int,
                    default=2850)

parser.add_argument("-ycenter",
                    help="y value for center of box slice",
                    type=int,
                    default=1550)

parser.add_argument("-dx",
                    help="width of box slice",
                    type=int,
                    default=250)

parser.add_argument("-dy",
                    help="height of box slice",
                    type=int,
                    default=350)

parser.add_argument("-ext",
                    help="extension of fits file",
                    type=int,
                    default=1)

parser.add_argument("-keyword",
                    help="keyword to sort file list by",
                    type=str,
                    default='darktime')


class GifObj(object):
    """
    This class is intended to create an animated matplotlib plot given a list 
    of input parameters. There are two ways to use this; the first is from the 
    command line and the second is by importing this into your scripts and 
    instanstiaing the class with your desired inputs, then calling the
    GifObj.mkgif method
    """
    def __init__(self, path, suffix, x_center, y_center, dx, dy, ext, keyword):
        """
        Parameters
        ----------
        path
        suffix
        x_center
        y_center
        dx
        dy
        ext
        keyword
        """

        # Initialize parameters from keyword arguements
        self.path = path
        self.suffix = suffix
        self.x_center = x_center
        self.y_center = y_center
        self.dx = dx
        self.dy = dy
        self.ext = ext
        self.keyword = keyword
        # Initialize parameters for plotting the data
        self.im=None
        self.flist=[]
        self.img_data = []
        self.key_values = []
        self.possible_keywords = {'darktime': 0.,
                                    'exptime': 0.,
                                    'date-obs':'%Y-%m-%d %H:%M:%S',
                                    'useafter':'%b %d %Y %H:%M:%S' }

    def grab_data(self, fname):
        """
        Grab the data from the image
        """
        print(fname)
        hdulist = fits.open(fname)

        chip = hdulist[self.ext].data
        if self.x_center and self.y_center and self.dx and self.dy:
            y1 = self.y_center - self.dy
            y2 = self.y_center + self.dy
            x1 = self.x_center - self.dx
            x2 = self.x_center + self.dx
            # print(y1, y2, x1, x2)
            # TO DO remove try, except
            try:
                chip_slice = chip[y1:y2, x1:x2]
            except IndexError as e:
                raise e
            else:
                self.img_data.append(chip_slice)
                # return data
        else:
            self.img_data.append(chip)
            # return chip
        

    def updatefig(self, i):
        """
        Update the figure to display the ith frame.
        """
        #TODO: fix issue with updating title
        self.im.set_array(self.img_data[i])
        self.ax.set_title('{}, {} = {}'.format(
                                               os.path.basename(self.flist[i]),
                                               self.keyword,
                                               self.key_values[i]))
        return [self.im]

    def quick_sort(self):
        """Sort the list of files, data, and header values

        Zip each of three lists and sort through them
        by the header values (e.g. date-obs, exptime, darktime)
        """
        zipped = zip(self.flist, self.key_values, self.img_data)
        zipped_list = list(zipped)
        zipped_list.sort(key=lambda pair: pair[1]) # makes sure to sort by key
        self.flist, self.key_values, self.img_data = zip(*zipped_list)
        

    def grab_header_value(self, fname):
        """ 
        Grab the header value for the given file
        
        """
        if self.keyword in self.possible_keywords.keys():
            value = fits.getval(fname, keyword=self.keyword)
            if isinstance(value, str):
                try:
                    date = dt.datetime.strptime(value, self.possible_keywords[key])
                except ValueError as e:
                    try:
                        date = dt.datetime.strptime(value,'%B %d %Y %H:%M:%S')
                    except ValueError as e:
                        pass
                    else:
                        value = date
                else:
                    value = date
            self.key_values.append(value)
            print(self.key_values)
        # return value


    def mkgif(self, flist=None):
        """ Class method for creating the actual animated plot.
        
        Parameters
        ----------
        flist

        Returns
        -------

        """
        if not flist:
            print('Searching directory {}'.format(self.path+self.suffix))
            self.flist = glob.glob(self.path+self.suffix)
        else:
            self.flist = flist
        if self.flist:
            idx = int(len(self.flist) / 2) # use this index value to set the normalization and stretch of the image
            print('Found a total of {} images'.format(len(self.flist)))
            try:
                img_units = fits.getval(self.flist[0],keyword='bunit', ext=self.ext)
            except KeyError as e:
                img_units = ''
            #cpu_count = mp.cpu_count()
            # Grab the data 
            # with mp.Pool(processes = cpu_count) as pool:
            #     pool.map(self.grab_data, self.flist)
            # # Grab header information
            # with mp.Pool(processes = cpu_count) as pool:
            #     pool.map(self.grab_value, self.flist)
            for f in self.flist:
            #     # data = self.grab_data(f)
            #     # value = self.grab_value(f)
            #     # self.key_values.append(value)
            #     # self.img_data.append(data)
                self.grab_data(f)
                self.grab_header_value(f)
            
            self.quick_sort() #sort list by key_values passed cmd line arguments

            self.fig = plt.figure(figsize=(6,8))
            self.ax = self.fig.add_subplot(1, 1, 1)
            self.ax.set_title('{}, {} = {}'.format(os.path.basename(self.flist[0]), self.keyword ,self.key_values[0]))
            self.ax.set_xlabel('X [pix]')
            self.ax.set_ylabel('Y [pix]')
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes('right', size='5%', pad=0.1)
            # Get the normalization factor based of data from the midpoint of the list
            norm = ImageNormalize(self.img_data[idx], stretch=LinearStretch(), interval=ZScaleInterval())
            # Draw the first image
            self.im = self.ax.imshow(self.img_data[0], animated=True, origin='lower', norm=norm, cmap='Greys_r')
            cbar = self.fig.colorbar(self.im, cax=cax, orientation='vertical')
            cbar.set_label(img_units)
            ani = animation.FuncAnimation(self.fig, self.updatefig, frames=len(self.flist), interval=500, blit=False)
            plt.show()
        else:
            print('No files were found at {}'.format(self.path+self.suffix))



def parse_cmd_line():
    """
    Parse the cmd line args
    """
    args = parser.parse_args()
    path = args.path
    suffix = args.suffix
    x_center = args.xcenter
    y_center = args.ycenter
    dx = args.dx
    dy = args.dy
    ext = args.ext
    keyword = args.keyword
    return path, suffix, x_center, y_center, dx, dy, ext, keyword


def main():
    # Parse cmd line args, if none provided uses defaults set above
    path, suffix, x_center, y_center, dx, dy, ext, keyword = parse_cmd_line()

    # Create an instance of the gif class, initializing with the cmd line args
    gif_obj = GifObj(path, suffix, x_center, y_center, dx, dy, ext, keyword)

    # Call the method to create the actual gif
    gif_obj.mkgif()
    
if __name__ == '__main__':
    main()

