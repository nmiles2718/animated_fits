#!/usr/bin/env python
"""
Script for creating an animated matplotlib plot using user supplied cmd line
arguments or preset default values. 
"""
# TODO: Add functionaltiy for plotting with astropy WCS alginment.

import argparse
from astropy.io import fits
from astropy.wcs import WCS
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
from PIL import Image


mpl.rcParams['animation.writer']='imagemagick_file'
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
                    help="apply custom scaling vmin, vmax (default is z-scale)",
                    action='store_true')

parser.add_argument("-save",
                    help="filename (e.g. animated_fits.mp4) of mp4 to save "
                         "the animated plot",
                    type=str,
                    default=None)

parser.add_argument("-fps",
                    help="Frame rate for animation",
                    type=float,
                    default=2)


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
                 ext, keyword, scale, save, fps):
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
        self.fps = fps
        # Initialize parameters for plotting the data
        self.im = None
        self.fig = plt.figure(figsize=(5,4))
        self.ax = self.fig.add_subplot(1,1,1)
        self.flist = []
        self.img_data = []
        self.key_values = []
        self.possible_keywords = {'darktime': 0.,
                                  'exptime': 0.,
                                  'flashdur': 0.,
                                  'CRSIGMAS':'',
                                  'expstart':0,
                                  'filter1':'',
                                  'date-obs': '%Y-%m-%d %H:%M:%S',
                                  'useafter': '%b %d %Y %H:%M:%S'}

    def grab_data(self, fname):
        """
        Grab the data from the image
        """
        if 'fits' in fname:
            hdulist = fits.open(fname)
            chip = hdulist[self.ext].data
        else:
            # cheap way to handle png or jpeg
            chip = Image.open(fname)
        if self.x_center and self.y_center and self.dx and self.dy:
            y1 = self.y_center - self.dy
            y2 = self.y_center + self.dy
            x1 = self.x_center - self.dx
            x2 = self.x_center + self.dx
            chip_slice = chip[y1:y2, x1:x2]
            self.img_data.append(chip_slice)
        else:
            self.img_data.append(chip)

    def updatefig(self, i):
        """
        Update the figure to display the ith frame.
        """

        self.im.set_array(self.img_data[i])
        if self.keyword:
            self.ax.set_title('{}, {} = {}'.
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
            try:
                value = fits.getval(fname, keyword=self.keyword, ext=0)
            except KeyError as e:
                print('Keyword not in primary header, checking science header')
                value = fits.getval(fname, keyword=self.keyword, ext=1)
            if isinstance(value, str):
                try:
                    time_obs = fits.getval(fname, keyword='time-obs', ext=0)
                except KeyError:
                    time_obs = fits.getval(fname, keyword='time-obs', ext=1)
                try:
                    date = \
                        dt.datetime.strptime(value+' '+time_obs,
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
                print('Keyword for image units not found.')
                img_units = input('Please enter units: ')

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
                if self.keyword and 'fits' in f:
                    self.grab_header_value(f)



            if self.keyword and 'fits' in self.flist[0]:
                self.quick_sort()  # Sort the data
                # Initialize some plot things
                self.ax.set_title('{}, {} = {}'.
                                  format(os.path.basename(self.flist[0]),
                                         self.keyword,
                                         self.key_values[0]))


            else:
                self.ax.set_title('{}'.format(os.path.basename(self.flist[0])))
            self.ax.set_xlabel('X [pix]')
            self.ax.set_ylabel('Y [pix]')

            # Grab index value to set the normalization/stretch
            # idx = int(len(self.flist)/2)
            idx = 0

            # Set the scaling to match ds9 z-scale with linear stretch
            # Not to avoid issues of import
            if self.scale:
                vmin = input('vmin for scaling: ')
                vmax = input('vmax for scaling: ')
                norm = ImageNormalize(self.img_data[idx],
                                      stretch=LinearStretch(),
                                      vmin=float(vmin),
                                      vmax=float(vmax))
            else:
                norm = ImageNormalize(self.img_data[idx],
                                      stretch=LinearStretch(),
                                      interval=ZScaleInterval())
            # Draw the first image
            if self.dx and self.dy and self.x_center and self.y_center:
                central_point = patches.Rectangle((patch_coords[0] ,
                                                   patch_coords[1] ),
                                                  width=1, height=1,
                                                  alpha=1.0, fill=False,
                                                  linewidth=1.75, color='red')
                self.ax.add_patch(central_point)
                self.im = self.ax.imshow(self.img_data[0],
                                         animated=True, origin='lower',
                                         norm=norm, cmap='gray',
                                         extent=[self.x_center - self.dx - 0.5,
                                                 self.x_center + self.dx + 0.5,
                                                 self.y_center - self.dy - 0.5,
                                                 self.y_center + self.dy + 0.5])
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
                                          interval=1000/self.fps, blit=False)

            # Either save the output or display it.
            if self.save and 'gif' in self.save:
                ani.save(filename=self.save, writer='imagemagick_file',
                         fps=self.fps, bitrate=150)
            elif self.save and 'mp4' in self.save:
                ani.save(filename=self.save, writer='ffmpeg',
                         dpi=180, fps=self.fps,bitrate=150)
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
    fps = args.fps
    return path, suffix, x_center, y_center, \
           dx, dy, ext, keyword, scale, save, fps


def main():
    """ Create the animation using user supplied cmd-line args (or defaults)
    """
    # Parse cmd line args, if none provided uses defaults set above
    path, suffix, x_center, y_center, dx, dy, ext, keyword, scale, save, fps =\
        parse_cmd_line()

    # Create an instance of the animation class
    # Initialize with the cmd line args
    animation_obj = AnimationObj(path, suffix, x_center,
                     y_center, dx, dy, ext, keyword, scale, save, fps)

    # Call the method to create the actual animation
    animation_obj.animate()
    
if __name__ == '__main__':
    main()
