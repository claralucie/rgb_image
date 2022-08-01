#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 13:29:42 2022

@author: s1929920
"""
import numpy as np
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import imread

import aplpy


from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename

from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D

from photutils.detection import DAOStarFinder, find_peaks
#from detection import find_peaks

import os

#create path to where fits files are stored
path = "/storage/teaching/SummerProjects2022/s1929920/derek_ceers_210722"
os.chdir(path)

def cutout(ra, dec, mos, size=5.):
    wcs = WCS(mos[0].header)
    wcs.sip = None
    
    if "CD1_1" in list(mos[0].header):
        cdelt = np.abs(mos[0].header["CD1_1"]*3600.)
        
    elif "CDELT1" in list(mos[0].header):
        cdelt = np.abs(mos[0].header["CDELT1"]*3600.)
        
    coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    
    cutout = Cutout2D(mos[0].data, coord, size/cdelt, wcs=wcs)
    
    return cutout

#filters used
filters = ["F115W", "F150W", "F200W", "F277W", "F356W", "F410M", "F444W"]

#images for each filter
images1 =  ["jwst_ceers_first_nircam_f115w_microJy_swarped.fits",
            "jwst_ceers_first_nircam_f150w_microJy_swarped.fits",
            "jwst_ceers_first_nircam_f200w_microJy_swarped.fits",
            "jwst_ceers_first_nircam_f277w_microJy_swarped.fits",
            "jwst_ceers_first_nircam_f356w_microJy_swarped.fits",
            "jwst_ceers_first_nircam_f410m_microJy_swarped.fits",
            "jwst_ceers_first_nircam_f444w_microJy_swarped.fits"]


#change to the coordinates of required object 
ra = 215.039074 
dec = 53.002678    
size = 2 #size of image

all_axes = []

axes = plt.subplot(1, 1, 1)

#change to the filter required from list above
mos = fits.open(images1[1]) 
cut = cutout(ra, dec, mos, size)

axes.imshow(np.flipud(cut.data), cmap='binary_r',
            norm = Normalize(vmin=np.percentile(cut.data, 0.5),
                             vmax=np.percentile(cut.data, 99.9))) #can change autocut here

plt.show()
plt.close()

#or to loop through all filters
"""
for i in range(len(images1)):
    mos = fits.open(images1[i])
    fig = plt.figure(figsize=(10,10))
    axes = plt.subplot(1,1,1)
    
    cut = cutout(ra, dec, mos, size)
    
    axes.imshow(np.flipud(cut.data), cmap='binary_r',
                norm = Normalize(vmin=np.percentile(cut.data, 0.5),
                                 vmax=np.percentile(cut.data, 99)))

plt.show()
plt.close()
"""

#save as fits file
hdu = fits.PrimaryHDU(data=cut.data, header=mos[0].header)
hdu.header.update(cut.wcs.to_header())
#choose output file name
hdu.writeto('passive3_small_150.fits', overwrite=True)

#----------------------------------------------------------------------
