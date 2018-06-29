# Class-based python script for making animated plots using FITS data
####  In order to save mp4 and gifs, the following steps must be taken to install the FFMEPG MovieWrtier:
   1) Add conda-forge to your .condarc file
        - `conda config --add channels conda-forge`
   2) Install ffmpeg and its dependencies into your active environment 
        - `conda install ffmpeg`
### Usage
The code is intended to be used from the command line, but it may be used within an `iPython` session or other scripts as well. 
There are a variety of command line arguments that may be used and their help documentation may be accessed as follows:
```console
user@home:~$ python mkAnimation.py -h
usage: mkAnimation.py [-h] [-path PATH] [-suffix SUFFIX] [-xcenter XCENTER]
                      [-ycenter YCENTER] [-dx DX] [-dy DY] [-ext EXT]
                      [-keyword KEYWORD] [-scale] [-save SAVE] [-fps FPS]

Create an animation of slices from a list of fits images

optional arguments:
  -h, --help        show this help message and exit
  -path PATH        Path to a directory containing fits images (must include
                    trailing /)
  -suffix SUFFIX    Suffix of fits files (e.g. *flt.fits, j*flt.fits)
  -xcenter XCENTER  x value for center of box slice
  -ycenter YCENTER  y value for center of box slice
  -dx DX            width of box slice
  -dy DY            height of box slice
  -ext EXT          extension of fits file
  -keyword KEYWORD  keyword to sort file list by
  -scale            apply min/max scaling (default is z-scale)
  -save SAVE        filename (e.g. animated_fits.mp4) of mp4 to save the
                    animated plot
  -fps FPS          Frame rate for animation
```
Some things to keep in mind are:
- If the `-xcenter` and `-ycenter` arguments are supplied, then `-dx` and `-dy` must also be supplied.
   - The central pixel in the slice (xcenter, ycenter) will be highlighted with a red box.
- `-save fname.gif` will save the output as a gif with a default frame rate of 1 frame/second
- `-save fname.mp4` will save the output as an mp4 with a default frame rate of 1 frame/second

## Example Output with FITS data:
![](https://github.com/nmiles2718/animated_fits/blob/master/saturn_transit1.gif "Titan Transiting Saturn")
## Example Output with PNG:
![](https://github.com/nmiles2718/animated_fits/blob/master/FalconHeavyTakeOff.gif "Falcon Heavy Takeoff")

