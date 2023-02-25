#!/usr/bin/python3
# -*- coding: utf-8 -*-
# calculates flux-weighted optical depth scale

import sys
import re
import os
import math
import numpy as np

import PoWR

#------------------------------------------------------------------------
def main(arg1 = '', arg2 = 1., arg3 = 128.):

#    modpath=os.path.expanduser(curpath)

#   Default X-ray integration range in Angstroem
    wstart=1.
    wend=128.

    modname = 'MODEL'
#    modname = ''
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

#    print ('MODEL', modpath, modname)

    emcoli  = PoWR.getMSREADVar(modpath, 'EMCOLI',  'float', modname)
    xlambda = PoWR.getMSREADVar(modpath, 'XLAMBDA', 'float', modname)
    fweight = PoWR.getMSREADVar(modpath, 'FWEIGHT', 'float', modname)
    rstar   = PoWR.getMSREADVar(modpath, 'RSTAR',   'float', modname)
    NF = len(xlambda)

    c = 2.99792458E5
#   Constant hc = h * c in erg * Angstroem
    hc = 1.98648E-8
    XLSUN = 3.85E33
    #one solar mass per year in cgs:
    XMSUNPYR = 6.303E25

    kmax = 0
    intflux = 0.
    bolflux = 0.
    for k in range (NF):
        fluxk = max(emcoli[k], 0.)
        if (xlambda[k] >= wstart) and (xlambda[k] <= wend):
            kmax = k
            intflux = intflux + fluxk * fweight[k]
        if (xlambda[k] >= wstart):
            kmax = k
            bolflux = bolflux + fluxk * fweight[k]

    xoverbol = intflux / bolflux

    Lx = 4 * math.pi * (rstar)**2 * math.pi * intflux
    Lbol = 4 * math.pi * (rstar)**2 * math.pi * bolflux

    logxoverbol = np.log10(xoverbol)
    logLx = np.log10(Lx / XLSUN)
    logLbol = np.log10(Lbol / XLSUN)

    return float('{:.3}'.format(logLx)), float('{:.3}'.format(logLbol)), float('{:.3}'.format(logxoverbol))

#-----------------------------------------------------------------------

if __name__ == "__main__":
    carg1 = ''
    carg2 = 1. 
    carg3 = 128. 
    if (len(sys.argv) > 1):
        carg1 = sys.argv[1]
    if (len(sys.argv) > 2):
        carg2 = sys.argv[2]
    if (len(sys.argv) > 3):
        carg3 = sys.argv[3]
    res1, res2, res3 = main(carg1, str(carg2), str(carg3))
    print ('Integrating X-ray flux between '+str(carg2)+' and '+str(carg3)+' Angstroem:')
    print ('log Lx/Lsun   = '+str(res1))
    print ('log Lbol/Lsun = '+str(res2))
    print ('log Lx/Lbol   = '+str(res3))
