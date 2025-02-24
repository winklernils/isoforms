#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:37:45 2019

Routines for saving and loading data, interpolation

@author: Oliver M. Drozdowski
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def SAVEDATA(listofarrays, listofnames, filepath, header=None):
    """
    If listofarrays is a list of arrays and the listofnames is a list of names, the data will be saved with
    corresponding names. If listofnames is just a single string, the index of the array is added to
    differentiate. If listofarrays is a single array, it is saved with the corresponding name.
    The file is saved with the given header (as string starting with #n for every line) at
    " *filepath*/*name*.txt"
    """
    if isinstance(listofarrays, list):
        listofarrays = listofarrays 
        if not isinstance(listofnames, list):
            # create listofnames of equal size
            names = [listofnames + str(i) for i,j in enumerate(listofarrays)]
        else:
            # assume equal size if both lists
            names = listofnames
    else: 
        # create 1-lists for single arrays
        listofarrays = [listofarrays]
        names = [listofnames]
    if np.array(listofarrays[0]).ndim > 2:
        print("ERROR! Trying to save array with >2 dimensions to human-readable format.")
        return -1
        
        
    for index, array in enumerate(listofarrays):
        # Write the array to disk
        name = names[index]
        filepath_temp = filepath + '/' + str(name) + '.txt'
        with open(filepath_temp, 'w') as outfile:
            # Write the header Any line starting with "#" will be ignored by numpy.loadtxt
            if header is not None:
                outfile.write(header)
            #print(array)
            #print(name)
            #print(name[0])
            np.savetxt(outfile, array)
        print("Saved: " + str(filepath_temp))
    return 0

def LOADDATA(listoffilepaths):
    """
    If listofpaths is a list of paths, every file will be read successively and returned
    as a numpy array in a list in the same order
    If listofpaths is a single path the array is returned as a single-element list
    """
    returnlist = []
    if not isinstance(listoffilepaths, list):
        listoffilepaths = [listoffilepaths]
    for filepath in listoffilepaths:
        data = np.loadtxt(filepath)
        returnlist.append(data)
        print("Loaded: " + str(filepath))
    return returnlist


def INTERPOLATE(positions, Xs_vel, velocityfield, plot_debug = False):
    """
    For a numerical "velocityfield" given at x-coordinates "Xs_vel" 
    we need to interpolate the velocity of the position given in "positions".
    We take the linear interpolation of the velocity field. If we are outside
    the range of [Xs_vel[0],Xs_vel[-1]] the value is set to 0. Ideally this
    should not happen. Xs_vel has to be in ascending order
    """    
    # Chose velocityfield at time-index
    vel_at_t = velocityfield
    # Make positions and Xs_vel broadcastable to get a boolean matrix at which
    # indices "positions" is smaller than "Xs_vel"
    pos = np.reshape(positions, (positions.shape[0],1))
    xs = np.reshape(Xs_vel, (1,Xs_vel.shape[0]))
    pos_smaller_than_xs = np.less(pos,xs)
    # Find the first index at which position is smaller than X-value
    # For pos < Xs[0], we get 0 and for pos >= Xs[-1] we get 0
    index = pos_smaller_than_xs.argmax(axis=1)
    # Calculate the fraction in this interval
    # For pos < Xs[0]: fraction > 1
    # for pos > Xs[0]: fraction < 0
    fraction = (positions - Xs_vel[index-1]) / (Xs_vel[index]-Xs_vel[index-1])
    # Interpolate velocity for positions
    vel = (1.0 - fraction) * vel_at_t[index-1] + fraction * vel_at_t[index]
    # Identify indices outside the Xs_vel interval and set to 0
    index_problem = np.logical_or(fraction < 0.0, fraction > 1.0)
    vel[index_problem] = 0.0
    if plot_debug:    
    #if True:    
        plt.plot(Xs_vel,vel_at_t)
        plt.scatter(positions,vel)
        plt.show()
    return vel

