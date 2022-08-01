#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 16:16:34 2022

@author: s1929920
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import aplpy
from astropy.wcs import WCS
from aplpy import make_rgb_image

#define file path for where the fits files are stored
path = "/storage/teaching/SummerProjects2022/s1929920/derek_ceers_210722"
os.chdir(path)

#filters used
filters = ["WFC3 F105W", "WFC3 F125W", "WFC3 F140W", "WFC3 F160W", "ACS F606W", "ACS F814W"]

#create cube from 3 fits files, 3 different filters
#save cube as a fits file
cube = aplpy.make_rgb_cube(['passive3_444.fits', 'passive3_200.fits', 'passive3_150.fits'],
                           'passive3_small_cube.fits')

#create an rgb image from the cube
""" 
stretch_r (red): choose scale - can be linear, power, log, sqrt
exponent_g (green): change power exponent
pmin_b (blue): change the minimum value pixel for the colour scale. Can help to reduce blue background noise
"""
im = aplpy.rgb.make_rgb_image('passive3_cube.fits', 'passive3_rgb.png', stretch_r = 'linear', stretch_g = 'linear', stretch_b = 'linear', pmin_g=40, pmin_b=60)

f = aplpy.FITSFigure('passive3_rgb.png')

f.show_rgb()
