#!/usr/bin/python3
# -*- coding: utf-8 -*-
# calculates flux-weighted optical depth scale

import sys
import re
import os
import math
import numpy as np
from matplotlib import rc
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import PoWR
rc('text', usetex=True)

#------------------------------------------------------------------------
def main(arg1 = '', arg2 = 3600., arg3 = 7000.):


    wstart = 3600.
    wend = 7000.
    modname = 'MODEL'
    if (len(arg1) > 0):
        modname=arg1

    if (len(arg2) > 0):
        delimiter=arg2
        if (delimiter[0:1].isdigit()):
            wstart = float(delimiter)

    if (len(arg3) > 0):
        delimiter=arg3
        if (delimiter[0:1].isdigit()):
            wend = float(delimiter)

    if modname[0] not in ['~','/']:
        modpath = os.getcwd()+'/'
        modFilename = modpath+modname
    else:
        modFilename = os.path.expanduser(modname)

    modpath, modname = os.path.split(modFilename)

    #Path to the plot file we want to read
    cformpath = modpath + '/formal.plot'
    cformFileHandle = open(cformpath, 'r')
    #Read in the actual dataset into python variables
    # The second argument is an offset inside the plot file
    # In this example, we want the first "PLOT:" dataset after finding "* OPT" 
    xfdata, yfdata = PoWR.readWRPlotDataset(cformFileHandle, '* OPT')
    cformFileHandle.close()

    #Create the figure and axjes objects
    #Note: We don't need to use addWRPlotFigure, but it cared about a bunch of nice formatting stuff
    # this command will create a figure of 15cm x 10cm
    f1 = PoWR.addWRPlotFigure(plt, 1, 15, 10)
    ax = f1.add_axes([0.,0.,1.0,1.0])

    myfs=14
    ax.tick_params(which='both', left=True, bottom=True, right=True, top=True)
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=12)
    ax.tick_params(which='minor', length=6)
    ax.xaxis.set_minor_locator(plt.MultipleLocator(100))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(r'$\lambda$\,[\AA]', fontsize=myfs)
    ax.set_ylabel(r"normalized Flux", fontsize=myfs)

    ax.plot(xfdata, yfdata, color='b')

    #Read an array from the MODEL mass storage file
    taurcont = PoWR.getMSREADVar(modpath, 'TAURCONT', 'float', modname)
    ND = len(taurcont)
    taurmax = taurcont[ND-1]

    #Read a few other values from the modinfo.kasdefs file
    #tstar: Teff at innermost tauross_cont (typically 20)
    tstar   = PoWR.getModinfoVar(modpath, 'TEFF',   'float')
    #t23: Teff at tauross = 2/3
    t23     = PoWR.getModinfoVar(modpath, 'T23',    'float')
    logl    = PoWR.getModinfoVar(modpath, 'LOGL',   'float')

    ax.text(0.05, 0.9, r'$T_{2/3}='+'{:.2f}'.format(t23/1000.)+r'\,$kK, '+
                       r'$T(\tau_{\mathrm{R,c}}='+'{:.0f}'.format(taurmax)+r') = '+'{:.2f}'.format(tstar/1000.)+r'\,$kK, '+
                       r'$\log\,(L/L_\odot) = '+'{:.2f}'.format(logl)+'$',
             horizontalalignment='left', verticalalignment='bottom',
             transform = ax.transAxes, fontsize=myfs)

    return f1.savefig('fig-example-plotspec.pdf', bbox_inches = 'tight', pad_inches = 0.05, transparent=True)

#-----------------------------------------------------------------------

if __name__ == "__main__":
    carg1 = ''
    carg2 = 3600. 
    carg3 = 7000. 
    if (len(sys.argv) > 1):
        carg1 = sys.argv[1]
    if (len(sys.argv) > 2):
        carg2 = sys.argv[2]
    if (len(sys.argv) > 3):
        carg3 = sys.argv[3]
    res = main(carg1, str(carg2), str(carg3))
    print ('Created normalized spectrum plot for wavelengths between '+str(carg2)+' and '+str(carg3)+' Angstroem.')
