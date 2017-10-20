#!/usr/bin/env python
import argparse
from astropy.io import fits
from astropy.visualization import LogStretch, LinearStretch, ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
import datetime as dt
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

# set up command line arguments for directory containing images
# NOTE: If the option isnt used the hotpix variable is defined to be False.
# If it is used then it is stored as true.
parser = argparse.ArgumentParser(\
            description="Create an animation of slices from a list of fits images")

parser.add_argument("-path", \
                    help="Path to a directory containing fits images (must include trailing /)",\
                    type=str, \
                    default='../test_data')

parser.add_argument("-suffix", \
                    help="Suffix of fits files (e.g. *flt.fits, j*flt.fits)", \
                    type=str, \
                    default='*flt.fits')

parser.add_argument("-xcenter", \
                    help="x value for center of box slice", \
                    type=int, \
                    default=2850)

parser.add_argument("-ycenter", \
                    help="y value for center of box slice", \
                    type=int, \
                    default=1550)

parser.add_argument("-dx", \
                    help="width of box slice", \
                    type=int, \
                    default=250)

parser.add_argument("-dy", \
                    help="height of box slice", \
                    type=int, \
                    default=350)

parser.add_argument("-ext", \
                    help="extension of fits file", \
                    type=int, \
                    default=1)

parser.add_argument("-keyword", \
                    help="keyword to sort file list by", \
                    type=str, \
                    default='darktime')

# parser.add_argument("--save", action="store_true", help="save the figure in the path provided.")

# Define some variables in the global scope
i=0
im=None
flist=[]
darktime=[]
fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(1, 1, 1)


def parse_cmd_line():
    """Parse the cmd line args
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

def quick_sort(list1, list2):
    """Sort list1 and list2 pairs in asecnding order by list2

    Keyword args:
    list1 -- list of filenames
    list2 -- correponding list of dates retrieved from given files

    Returns:
    file_list -- original list sorted by date in ascending order
    dates -- original list sorted by date sorted in ascending order
    """
    zipped = zip(list1, list2)
    zipped_list = list(zipped)
    zipped_list.sort(key=lambda pair: pair[1])
    list1, list2 = zip(*zipped_list)
    return list1, list2


def grab_data(fname, ext, x_center, dx, y_center, dy):
    """Grab the data from the image
    """
    try:
        hdulist = fits.open(fname)
    except FileNotFoundError as e:
        raise FileNotFoundError('{} not found'.format(fname))
    else:
        chip = hdulist[ext].data
        if x_center and y_center and dx and dy:
            y1 = y_center - dy
            y2 = y_center + dy
            x1 = x_center - dx
            x2 = x_center + dx
            # print(y1, y2, x1, x2)
            # TO DO remove try, except 
            try:
                chip_slice = chip[y1:y2, x1:x2]
            except IndexError as e:
                raise IndexError('{}, {}'.format((x1, x2), (y1, y2)))
            else:
                return chip_slice
        else:
            return chip
    # finally:
    #     return chip[1300:1800,2500:3200]


def updatefig(i, ext, x_center, dx, y_center, dy):
    """Update the figure to display the ith frame.
    """
    im.set_array(grab_data(flist[i], ext, x_center, dx, y_center, dy))
    ax.set_title('{}, darktime = {}'.format(os.path.basename(flist[i]), key_values[i]))
    return [im]


def grab_value(f, key):
    if key == 'darktime':
        value = fits.getval(f, keyword=key)
    elif key == 'exptime':
        value = fits.getval(f, keyword=key)
    elif key == 'date-obs':
        value_str = fits.getval(f, keyword=key)
        fmt = '%Y-%m-%d %H:%M:%S'
        value = dt.datetime.strptime(value_str, fmt)
    elif key == 'useafter':
        value_str = fits.getval(f, keyword=key)
        fmt = '%b %d %Y %H:%M:%S'
        try:
            value = dt.datetime.strptime(value_str, fmt)
        except ValueError as e:
            try:
                value = dt.datetime.strptime(value_str,'%B %d %Y %H:%M:%S')
            except ValueError as e:
                pass

    return value


def mkgif():
    global im 
    global flist
    global key_values
    global ax
    global i

    path, suffix, x_center, y_center, dx, dy, ext, keyword = parse_cmd_line()
    print('Searching directory {}'.format(path+suffix))
    flist = glob.glob(path+suffix)
    if flist:
        key_values =[]
        print('Found a total of {} images'.format(len(flist)))
        try:
            img_units = fits.getval(flist[0],keyword='bunit', ext=1)
        except KeyError as e:
            img_units = ''
        for f in flist:
            value = grab_value(f, keyword)
            key_values.append(value)

        flist, key_values = quick_sort(flist,key_values) #sort list by darktime
        ax.set_title('{}, {} = {}'.format(os.path.basename(flist[0]), keyword ,key_values[0]))
        ax.set_xlabel('X [pix]')
        ax.set_ylabel('Y [pix]')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        # im_data = grab_data(flist[0],ext, x_center, dx, y_center, dy)
        idx = int(len(flist)/2)
        norm = ImageNormalize(grab_data(flist[idx], ext, x_center, dx, y_center, dy), stretch=LinearStretch(), interval=ZScaleInterval())
        im = ax.imshow(grab_data(flist[0], ext, x_center, dx, y_center, dy), extent=[x_center-dx, x_center+dx, y_center-dy, y_center+dy],animated=True, origin='lower', norm=norm, cmap='Greys_r')
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        cbar.set_label(img_units)
        ani = animation.FuncAnimation(fig, updatefig, fargs=(ext, x_center, dx, y_center, dy), frames=len(flist), interval=250, blit=False)
        plt.show()
    else:
        print('No files were found at {}'.format(path+suffix))

if __name__=="__main__":
    mkgif()