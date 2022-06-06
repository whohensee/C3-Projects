# Filename: c3_lightcurve_functions.py
# Author : William Hohensee
# 
"""
Contains functions used in analyzing the lightcurve of PTF and iPTF transients



"""

# Imports
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import log10, floor
# from sympy import *
from scipy.special import factorial
import glob
import math
import pandas as pd


# ####
# Functions
# ####

def form_transient_names(filename):
    name = filename
    name = name[12:-7]
    name = 'PTF' + name
    name_PTF = name
    name = 'i' + name
    name_iPTF = name

    return name_PTF, name_iPTF

def analyze_lightcurve(filename, plot=False):
    """
    Filename: path to a .csv file
    plot: True to display a graph


    returns:
    risetime: time from halfway rising to peak
    fadetime: time from peak to halfway falling
    """
    name_PTF, name_iPTF = form_transient_names(filename)
    sn_name = filename[12:-7]

    data = np.genfromtxt(filename, skip_header=1, skip_footer=1,usecols=(0,2,3))

    if np.size(data) < 5:
        return -1, -1, -1
    
    # plt.figure(dpi=150)


    time_values = data[:,0]
    magnitudes = data[:,1]
    magnitude_uncertainties = data[:,2]

    x = time_values - time_values[0]
    y = magnitudes
    y_errors = magnitude_uncertainties

    # calculate and plot the interpolation
    interp_xvals = np.linspace(0, np.max(x), 10000)
    interp_yvals = np.interp(interp_xvals, x, y)

    peak_magnitude = np.min(y)
    max_index = np.argwhere(y == np.min(y))
    max_index = max_index[0]        #to eliminate cases where the peak occurs twice
    half_brightness = peak_magnitude + 2.5 * np.log10(2)

    good_xvals = np.array([])
    for i, xval in enumerate(interp_xvals):
        if np.abs(interp_yvals[i] - half_brightness) < 0.05:
            good_xvals = np.append(good_xvals, xval)
    # print(good_xvals)         # for testing

    
    if any(good_xvals[np.argwhere(good_xvals<x[max_index])]):
        half_rising_JD = np.average(good_xvals[np.argwhere(good_xvals<x[max_index])])
    else:
        half_rising_JD = -2
    if any(good_xvals[np.argwhere(good_xvals>x[max_index])]):
        half_falling_JD = np.average(good_xvals[np.argwhere(good_xvals>x[max_index])])
    else:
        half_falling_JD = -3

    risetime = float(x[max_index]) - half_rising_JD
    falltime = half_falling_JD - float(x[max_index])

    risetime = round(risetime, 2)
    falltime = round(falltime, 2)

    if plot:
        # plot the data
        plt.figure(dpi=150)
        fig, ax = plt.subplots(dpi = 150)
        ax.invert_yaxis()           # because lower magnitude is brighter

        # plot the points and error bars
        plt.errorbar(x,y, yerr = y_errors, marker = "o", linestyle="none", markersize = 2, capsize = 1.5, elinewidth = 0.5)

        #plot the half brightness line
        plt.axhline(half_brightness, color="g")

        #plot the interpolated points
        plt.errorbar(interp_xvals, interp_yvals, color="r", linewidth=0.5)

        #label the graph
        titletext = f"Lightcurve for transient {sn_name}"
        plt.title(titletext)
        plt.xlabel("Elapsed Time (days)")
        plt.ylabel("Magnitude")
        textstr = f'risetime: {risetime} \n falltime: {falltime}'
        ax.text(0.5, 0.1, textstr, ha='center', transform=ax.transAxes)

    return risetime, falltime, peak_magnitude

