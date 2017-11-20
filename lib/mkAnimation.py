#!/usr/bin/env python
"""
Script for creating an animated matplotlib plot using user supplied cmd line
arguments or preset default values. 
"""
import argparse
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.visualization import LinearStretch,ZScaleInterval,\
    LogStretch,SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import datetime as dt
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os

mpl.rcParams['animation.writer']='imagemagick'
mpl.rcParams['savefig.bbox']='standard'
parser = argparse.ArgumentParser(description="Create an animation of slices"
                                             " from a list of fits images")

parser.add_argument("-path",
                    help="Path to a directory containing "
                         "fits images (must include trailing /)",
                    type=str,
                    default='../test_data/')

parser.add_argument("-suffix",
                    help="Suffix of fits files (e.g. *flt.fits, j*flt.fits)",
                    type=str,
                    default='*flt.fits')

parser.add_argument("-xcenter",
                    help="x value for center of box slice",
                    type=int)

parser.add_argument("-ycenter",
                    help="y value for center of box slice",
                    type=int)

parser.add_argument("-dx",
                    help="width of box slice",
                    type=int)

parser.add_argument("-dy",
                    help="height of box slice",
                    type=int)

parser.add_argument("-ext",
                    help="extension of fits file",
                    type=int,
                    default=1)

parser.add_argument("-keyword",
                    help="keyword to sort file list by",
                    type=str,
                    default=None)

parser.add_argument("-scale",
                    help="apply min/max scaling (default is z-scale)",
                    action='store_true')

parser.add_argument("-save",
                    help="filename (e.g. animated_fits.mp4) of mp4 to save "
                         "the animated plot",
                    type=str,
                    default=None)


class AnimationObj(object):
    """
    This class is intended to create an animated matplotlib plot given a list 
    of input parameters. There are two ways to use this; the first is from the 
    command line and the second is by importing this into your scripts and 
    instantiating the class with your desired inputs, then calling the
    AnimationObj.animate() method with a set of files.
    """
    def __init__(self, path, suffix,
                 x_center, y_center, dx, dy,
                 ext, keyword, scale, save):
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

        # Initialize parameters from keyword arguments
        self.path = path
        self.suffix = suffix
        self.x_center = x_center
        self.y_center = y_center
        self.dx = dx
        self.dy = dy
        self.ext = ext
        self.keyword = keyword
        self.scale = scale
        self.save = save
        # Initialize parameters for plotting the data
        self.im = None
        self.fig = plt.figure(figsize=(8,6))
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.flist = []
        self.img_data = []
        self.key_values = []
        self.possible_keywords = {'darktime': 0.,
                                  'exptime': 0.,
                                  'flashdur': 0.,
                                  'date-obs': '%Y-%m-%d %H:%M:%S',
                                  'useafter': '%b %d %Y %H:%M:%S'}

    def grab_data(self, fname):
        """
        Grab the data from the image
        """
        hdulist = fits.open(fname)
        chip = hdulist[self.ext].data
        if self.x_center and self.y_center and self.dx and self.dy:
            y1 = self.y_center - self.dy
            y2 = self.y_center + self.dy
            x1 = self.x_center - self.dx
            x2 = self.x_center + self.dx
            chip_slice = 0.5079*chip[y1:y2, x1:x2]
            self.img_data.append(chip_slice)
        else:
            self.img_data.append(chip)

    def updatefig(self, i):
        """Update the figure to display the ith frame.
        """
        
        self.im.set_array(self.img_data[i])
        if self.keyword:
            self.ax.set_title('{}(flashed), {} = {}'.
                format(os.path.basename(self.flist[i]),
                              self.keyword,
                              self.key_values[i]))
        else:
            self.ax.set_title('{}'.format(os.path.basename(self.flist[i])))
        return [self.im]

    def quick_sort(self):
        """Sort the list of files, data, and header values

        Zip each of three lists and sort through them
        by the header values (e.g. date-obs, exptime, darktime)
        """
        zipped = zip(self.flist, self.key_values, self.img_data)
        zipped_list = list(zipped)
        zipped_list.sort(key=lambda pair: pair[1])  # makes sure to sort by key
        self.flist, self.key_values, self.img_data = zip(*zipped_list)

    def grab_header_value(self, fname):
        """ 
        Grab the header value for the given file
        
        """
        if self.keyword in self.possible_keywords.keys():
            value = fits.getval(fname, keyword=self.keyword)
            if isinstance(value, str):
                try:
                    date = \
                        dt.datetime.strptime(value,
                                             self.possible_keywords[
                                                 self.keyword])
                except ValueError:
                    try:
                        date = dt.datetime.strptime(value, '%B %d %Y %H:%M:%S')
                    except ValueError:
                        pass
                    else:
                        value = date
                else:
                    value = date
            self.key_values.append(value)

    def animate(self, flist=None):
        """ Class method for creating the actual animated plot.
        
        Parameters
        ----------
        flist -- User specified list of files, only used when run interactively

        Returns
        -------
        """
        if not flist:
            print('Searching directory {}'.format(self.path+self.suffix))
            self.flist = glob.glob(self.path+self.suffix)
        else:
            self.flist = flist
        if self.flist:

            print('Found a total of {} images'.format(len(self.flist)))
            try:
                img_units = fits.getval(self.flist[0],
                                        keyword='bunit',
                                        ext=self.ext)
            except KeyError:
                img_units = ''
            # TODO: figure how to parallelize the opening of fits files
            '''
            Currently fails because you cannot fork a GUI process. I need
            to break up the GUI and data extraction into to different sections, 
            then communicate via pipes. Will work on later.
            cpu_count = mp.cpu_count()
            Grab the data
            with mp.Pool(processes = cpu_count) as pool:
               print(pool.map_async(self.grab_data, self.flist))
            # Grab header information
            with mp.Pool(processes = cpu_count) as pool:
                pool.map(self.grab_header_value, self.flist)
            '''

            for f in self.flist:
                self.grab_data(f)
                if self.keyword:
                    self.grab_header_value(f)
            
            if self.keyword:
                self.quick_sort()  # Sort the data
                # Initialize some plot things
                self.ax.set_title('{}(flashed), {} = {}'.
                                  format(os.path.basename(self.flist[0]),
                                         self.keyword,
                                         self.key_values[0]))

                
            else:
                self.ax.set_title('{}'.format(os.path.basename(self.flist[0])))
            self.ax.set_xlabel('X [pix]')
            self.ax.set_ylabel('Y [pix]')
            # Offset to account for overscan columns in raw images.
            # 1522, 1773 blob location
            offset=0
            if 'raw' in self.suffix:
                offset=24
            rect = patches.Rectangle((3490+offset,340),40,40,
                                     linewidth=1.5,
                                     edgecolor='r',
                                     facecolor='none')
            self.ax.add_patch(rect)
            # Grab index value to set the normalization/stretch
            idx = int(len(self.flist)/2)


            # Set the scaling to match ds9 z-scale with linear stretch
            if self.scale:
                vmin = input('vmin for scaling: ')
                vmax = input('vmax for scaling: ')
                norm = ImageNormalize(self.img_data[idx],
                                      stretch=LinearStretch(),
                                      vmin=float(vmin),
                                      vmax=float(vmax))
            else:
                norm = ImageNormalize(self.img_data[idx],
                                      stretch=SqrtStretch(),
                                      # interval=ZScaleInterval())
                                      vmin=0,
                                      vmax=175)
            # Draw the first image
            if self.dx and self.dy and self.x_center and self.y_center:
                self.im = self.ax.imshow(self.img_data[0],
                                         animated=True, origin='lower',
                                         norm=norm, cmap='gray',
                                         extent=[self.x_center - self.dx,
                                                 self.x_center + self.dx,
                                                 self.y_center - self.dy,
                                                 self.y_center + self.dy])
            else:
                self.im = self.ax.imshow(self.img_data[0],
                                         animated=True, origin='lower',
                                         norm=norm, cmap='gray')

            # Add a nice colorbar with a label showing the img units
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes('right', size='5%', pad=0.1)
            cbar = self.fig.colorbar(self.im, cax=cax, orientation='vertical')
            cbar.set_label(img_units)
            # Start the animation
            ani = animation.FuncAnimation(self.fig, self.updatefig,
                                          frames=len(self.flist),
                                          interval=500, blit=False)

            if self.save and'gif' in self.save:
                ani.save(filename=self.save, writer='imagemagick',
                         fps=2, bitrate=300)
            elif self.save and 'mp4' in self.save:
                ani.save(filename=self.save, writer='ffmpeg',
                         dpi=180, bitrate=300)
            else:
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
    scale = args.scale
    save = args.save
    return path, suffix, x_center, y_center, dx, dy, ext, keyword, scale, save


def main():
    """ Create the animation using user supplied cmd-line args (or defaults)
    """
    # Parse cmd line args, if none provided uses defaults set above
    path, suffix, x_center, y_center, dx, dy, ext, keyword, scale, save = \
        parse_cmd_line()

    # Create an instance of the animation class, initializing with the cmd line args
    animation_obj = AnimationObj(path, suffix, x_center,
                     y_center, dx, dy, ext, keyword, scale, save)

    # Call the method to create the actual animation
    animation_obj.animate()
    
if __name__ == '__main__':
    main()
