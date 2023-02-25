#!/usr/bin/python3
# -*- coding: utf-8 -*-
#Functions to use PoWR output in python scripts

import os
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pdb
from scipy import interpolate

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Read a variable from the PoWR configuration file
def getPoWRConfig(var = 'POWR_WORK'):
#    configfile = os.path.expanduser("~/.powrconfig")
    varval = os.environ.get(var)
    if (len(varval) == 0):
        return False
    else:
        return varval

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Read the last n lines of a file via the tail os command
def tail(f, n, offset = 0):
    num = str(n + offset)
    proc = subprocess.Popen(['tail', '-n', num, f], stdout=subprocess.PIPE)
    lines = proc.stdout.readlines()
    return lines[:n]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#fixes small values from PoWR/WRPlot datasets missing an E in the
# format, i.e. 0.578754-108 is converted to 0.578754E-108
def fixreadsmall(numstring):
    if (len(numstring) > 4 and numstring[-4:-3] == '-'):
        newstring = numstring[:-4] + 'E-' + numstring[-3:]
    elif (len(numstring) > 4 and numstring[-4:-3] == '+'):
        newstring = numstring[:-4] + 'E+' + numstring[-3:]
    else:
        newstring = numstring

    return newstring

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Read variable from modinfo.kasdefs
def getModinfoVar(modelname, variable, format = 'string'):
    modinfoFile = open(modelname+'/modinfo.kasdefs')
    for aline in modinfoFile:
        cline = re.findall(r'\\VAR +'+variable+' +', aline)
        if cline:
            value = aline.split('=')[1].strip()
            if format == 'float':
                value = float(value)
            if format == 'int':
                value = int(value)
            return value

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Read ionizing flux from wruniq.out
def getModIonFlux(modelname, ion, format = 'string'):
    #Convert H I into wruniq label:
    if (ion == 'H I'): 
        ion = 'Ly Alpha'
    #Check if requested ion is in list of ions
    ions = ['Ly Alpha', 'He I', 'He II', 'O II']
    if (ion not in ions):
        print ('Error in getModIonFlux -- Invalid ion requested: ', ion)
        return false

    wruniqLastLines = tail(modelname+'/wruniq.out', 100)
    needle1 = 'Emergent Flux from COLI:'
    needle2 = 'EDGE'
    #We need to search for a line with needle1:
    #Then we need to find the first line after a line with needle2 follows the line for needle1


    doAnalyze = -1
    for aline in wruniqLastLines:
        searchLab = aline.strip()
        if (doAnalyze == -1):
            if (searchLab[:len(needle1)] == needle1):
                doAnalyze = 0
        if (doAnalyze == 0):
            if (searchLab[:len(needle2)] == needle2):
                doAnalyze = 1 
        if (doAnalyze == 1):
            if (searchLab[:len(ion)] == ion):
                #first number column after 35 characters
                #third column contains loq Q
                cline = aline[34:]
                value = cline.split()[2].strip() 
                if format == 'float':
                    value = float(value)
                if format == 'int':
                    value = int(value)
                return value

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Read all magnitudes from wruniq.out
#Returns a dictionary with all magnitudes identified by their PoWR string 
def getMagnitudes(modelname):

    wruniqLastLines = tail(modelname+'/wruniq.out', 100)
    needle1 = 'ABSOLUTE MONOCHROMATIC MAGNITUDES'
    needle2 = 'X-RAY VIEW'
    #We need to search for a line with needle1:
    #Then we need to read all lines with magnitudes until needle2
    #Magnitude lines are identifed by having 'A ' or 'MU' in a certain position (30th and 31st character) 

    doAnalyze = -1
    maglist = {}
    for aline in wruniqLastLines:
        searchLab = aline.strip()
        if (doAnalyze == -1):
            if (searchLab[:len(needle1)] == needle1):
                doAnalyze = 1
        if (doAnalyze == 1):
            if (searchLab[:len(needle2)] == needle2):
                break 
            if (aline[29:31] in ['A ','MU']):
                #current line contains a magnitude
                #third column contains loq Q
                clabel = aline[10:18].strip()
                cvalue = float(aline[32:42].strip())
                maglist[clabel] = cvalue

    return maglist

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Helper function: Read variable via msread from MODEL file
def getMSREADVar(modelname, variable, format = 'string', modfilename = 'MODEL'):
    modelfile = modelname+'/'+modfilename
    powrdir = getPoWRConfig('POWR_WORK')
    if (powrdir is False):
        print ('Fatal error: No PoWR installation found!')
        quit()
    msread = powrdir + 'proc.dir/msread.com'
    command = msread + ' ' + variable + ' ' + modelfile
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    p_status = p.wait()

    dataset=[]
    extract=-1
    for aline in output.decode().split('\n'):
        if extract >= 0:
            extract = extract + 1
        if aline.strip() == 'N=?':
            extract = 0
        if aline.strip() == 'FINISH':
            extract = -1
        if extract > 0:
            dataset.append(aline.strip())

    if format != 'string':
        for i in range(len(dataset)):
            if format == 'float': 
                test = fixreadsmall(dataset[i])
                dataset[i] = float(fixreadsmall(dataset[i])) 
            if format == 'int':
                dataset[i] = int(dataset[i]) 

    #for a single entry array, return the pure value instead
    if len(dataset) == 1:
        dataset = dataset[0]

    return dataset

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def cm2inch(value):
    return value/2.54

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def splinpo(xval, ydata, xdata):

    ipfunc = interpolate.interp1d(xdata, ydata, kind='cubic')
    yval = ipfunc(xval)
    return yval

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def lipo(xval, ydata, xdata):

    ipfunc = interpolate.interp1d(xdata, ydata, kind='linear')
    yval = ipfunc(xval)
    return yval

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Helper function: Set layout of figures to TeX style (WRPlot Times)
def setWRPlotEnv(plt):
    plt.rcParams["font.family"] = 'serif'
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = r'''
    %\usepackage{mathtools}
    \usepackage{amsmath}

    \usepackage{txfonts}
    \renewcommand{\familydefault}{\rmdefault}
    % more packages here
    '''

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Helper function: Add new figure with WRPlot-like cm dimensions
def addWRPlotFigure(plt,mynum=-99,xcm=15,ycm=10):
    #automatic number increasing if mynum is not given
    if (mynum == -99):
        mynum = plt.gcf().number + 1
    return plt.figure(num=mynum,figsize=(cm2inch(xcm),cm2inch(ycm)))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Helper function: Set layout of a (sub-)plot to WRPlot style
def setWRPlotStyle(ax, xminor='auto', xmajor='auto', yminor='auto', ymajor='auto'):
    ax.minorticks_on()
    ax.tick_params(which='both', left=True, bottom=True, right=True, top=True)
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
    if (xminor != 'auto'):
        ax.xaxis.set_major_locator(plt.MultipleLocator(xmajor))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(xminor))
        ax.yaxis.set_major_locator(plt.MultipleLocator(ymajor))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(yminor))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spline interpolation with maxima/minima only at the original datapoints
def powrsplinpo(x, ydata, xdata, calcdfdx=False):

    N = len(xdata)

    dn = xdata[-1] - xdata[0]
    for i in range(1,N):
        dx = xdata[i] - xdata[i-1]
        if (dx*dn <= 0):
            print ("FATAL ERROR: x values not in strictly monotonic order")
            return False

# Find the interval
    ivfound = False
    for i in range(1,N):
        if ((x - xdata[i-1]) * (x - xdata[i]) <= 0.):
            lip = i
            ivfound = True
            break

    if (not ivfound):
        print ("FATAL ERROR: X outside interpolation range")
        print ("xmin, xmax, x = ", min(xdata), max(xdata), x)
        return False

# Determination of the coefficients p1, p2, p3, p4 
# Set up the coefficient matrix
    d1=1./(xdata[lip] - xdata[lip-1])
    d2=d1*d1
    d3=d1*d2
    d23=d2/3.
    h11=d3
    h12=-d3
    h13=d23
    h14=2.*d23
    h21=-d1
    h22=2.*d1
    h23=-0.333333333333333
    h24=-0.666666666666666
    h31=-d3
    h32=d3
    h33=-2.*d23
    h34=-d23
    h41=2.*d1
    h42=-d1
    h43=0.666666666666666
    h44=0.333333333333333

# FOR THE BOUNDARY INTERVALS THE DERIVATIVE CANNOT EXTEND OVER THE BOUNDARY
    la=max(lip-2,0)
    lb=min(lip+1,N-1)

# FUNCTION TO BE INTERPOLATED: ydata
    f1 = ydata[lip-1]
    f2 = ydata[lip]

# Standard version: zentrierte Ableitungen an den Stuetzstellen. Das 
#  ist bei nicht aequidistanten Stuetzstellen fragwuerdig, verringert 
#  aber andererseits das Ueberschwingen insbesondere wenn man nicht MONO verwendet 
    f3 = (ydata[lip] - ydata[la]) / (xdata[lip] - xdata[la])
    f4 = (ydata[lb] - ydata[lip-1]) / (xdata[lb] - xdata[lip-1])

    s4 = (ydata[lip] - ydata[lip-1]) / ( xdata[lip] - xdata[lip-1] )

#   We are not in the first interval:
    if (la != lip-2):
        s3 = s4
    else:
        s3 = ( ydata[lip-1] - ydata[lip-2] ) / ( xdata[lip-1] - xdata[lip-2] )

#   We are not in the last interval:
    if (lb != lip+1):
        s5 = s4
    else:
        s5 = ( ydata[lip+1] - ydata[lip] ) / ( xdata[lip+1] - xdata[lip] )

    f3 = (np.sign(s3)+np.sign(s4)) * min(abs(s3),abs(s4),0.5*abs(f3))
    f4 = (np.sign(s4)+np.sign(s5)) * min(abs(s4),abs(s5),0.5*abs(f4))

#   Calculate polynomial coefficients: P(vector) = h (maxtrix) * f(vector)
    p1 = h11 * f1 + h12 * f2 + h13 * f3 + h14 * f4
    p2 = h21 * f1 + h22 * f2 + h23 * f3 + h24 * f4
    p3 = h31 * f1 + h32 * f2 + h33 * f3 + h34 * f4
    p4 = h41 * f1 + h42 * f2 + h43 * f3 + h44 * f4

#   Evaluation of the interpolation polynomial
    dxm = x - xdata[lip-1]
    dx = xdata[lip] - x
    y = (p1 * dxm * dxm + p2) * dxm + (p3 * dx * dx + p4) * dx

    if (calcdfdx):
        dfdx = 3. * p1 *  dxm * dxm + p2 - 3. * p3 * dx * dx - p4
        return y, dfdx
    else:
        return y

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Helper function: Smooth a curve by performing a spline interpolation
def smoothCurve(Xdata, Ydata):
    Xfine = np.linspace(np.array(Xdata).min(), np.array(Xdata).max(), 300)    
    Yfine = []
    for cX in Xfine:
        cY = powrsplinpo(cX, Ydata, Xdata)
        Yfine.append(cY)
    return Xfine, Yfine

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#Helper function: Slice out the range xstart to xend from an XY dataset
# This roughly corresponds to the WRPlot COMMAND X-CUT
def dataXcut(xstart, xend, Xdata, Ydata):
    xmin = min( Xdata )
    xmax = max( Xdata )
    icopys = 0
    Xnew = []
    Ynew = []
    if (xstart > xmin):
        ystart = powrsplinpo(xstart, Ydata, Xdata)
        Xnew.append( xstart )
        Ynew.append( ystart )
        while (Xdata[icopys] <= xstart):
            icopys = icopys + 1     
        
    icopye = len(Xdata)
    addend = False
    if (xend < xmax):
        yend = powrsplinpo(xend, Ydata, Xdata)
        addend = True
        while (Xdata[icopye-1] >= xend):
             icopye = icopye - 1

    for ic in range(icopys, icopye):
        Xnew.append( Xdata[ic] ) 
        Ynew.append( Ydata[ic] )
    
    if (addend):
        Xnew.append( xend ) 
        Ynew.append( yend )

    return Xnew, Ynew

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# read out a dataset in the WRPlot format and return it as x and y lists    
# Datasets start with "N=" and end with "FINISH", "END", "ENDE" or the next "N="    
# This method is particularly intended for reading out spectra from formal.plot
def readWRPlotDataset(lineset, keyword = "", dataset = 1):
    startkey = "N=" 
    endkeys = ["N=", "FINISH", "END"]

    xdata = []
    ydata = []    
    readindex = -1
    nskip = dataset - 1
    foundkey = (keyword == "")
    for curline in lineset:  
        if (not foundkey): 
            keypos = curline.find(keyword) 
            if (keypos == -1):
                continue
            else:
                foundkey = True

        if ((readindex == -1) and foundkey):
            readindex = 0

        if (readindex == 0):
            if (curline.strip().startswith(startkey)):
                if (nskip > 0):
                    nskip = nskip - 1
                else:
                    readindex = 1
            continue
        elif (readindex == 1 or readindex == 2):
#           now we need to read out pairs of lines

            rawline = curline.rstrip()
            for endkey in endkeys:
                if (rawline.strip().startswith(endkey)):
                    if (readindex == 2):
                        print ("FATAL ERROR: odd number of xy-lines")
                        return False
                    readindex = 99
                    break
            xynewset = rawline.split()
            if (readindex == 1):
                for xynew in xynewset:
                    xdata.append(float(xynew))
                readindex = 2
            elif (readindex == 2):
                for xynew in xynewset:
                    ydata.append(float(xynew))
                readindex = 1
        elif (readindex > 9):
            break

    if (not foundkey):
        print ('Could not find keyword ', keyword)

    if (readindex == 0):
        print ('Could not find dataset ', dataset)

    return xdata, ydata

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#HelperFunction: Read lines from a WRPLOT-like file
#                and ignore all comments and WRPlot-data like lines 
def filterWRPlotFile(filehandle):
    #return only lines with actual data from a file
    #but skip comment lines or WRPLOT statements
    commentchars=['*','#']
    commentwords=['N=?','END','FINISH','COMMAND']

    datalines = []

    for line in filehandle:
        #Do not store empty lines or comment lines
        ccheck = line.lstrip()
        if len(ccheck) == 0:
            continue
        if ccheck[:1] in commentchars:
            continue
        wcheck = line.split()
        if wcheck[0] in commentwords:
            continue
        datalines.append(line)

    return datalines

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
