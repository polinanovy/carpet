#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from datetime import datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

'''
Plotting delta alpha and signal-to-noise ratio for each object of the LHAASO catalog.
'''

res = 4.7 # Carpet resolution, deg
circle = np.pi * res**2 #Signal circle
DA_c = np.cos(4.7 * np.pi / 180) #Signal circle converted to delta alpha

class LHAASO():
    def __init__(self):
        self.name = None #Object name
        self.RA  = None #Right Ascension
        self.DEC = None #Declination
        self.unc95 = None #95% uncertainty
        self.detect = None #LHAASO detector type

class Carpet():
    def __init__(self):
        self.date = None #Event date in dd.mm.yy format
        self.time = None #Event time in hh:mm:ss format
        self.RA  = None #Right Ascension
        self.DEC = None #Declination
        self.Ne = None #relativistic particle number
        self.mjd = None #date+time
        self.mu = None #Muon detector counter

def CreateDir(save_dir):
    """
    Function for creating directory for saving the results.

    Parameters
    ----------
    save_dir : string
        Directory for saving the results.

    Returns
    -------
    None.

    """
    if not(os.path.exists(save_dir)):
        key = True
        while key:
            #ans = input("No path %s. Do you want to create path? (y/n) \n" %save_dir)
            #if (ans == "y"):
                subprocess.run("mkdir -p %s" %save_dir, shell=True)  
                key = not(key)
            #elif(ans == "n"):
            #    print("Exit")
            #    key = not(key)
            #    return
    return

def GetCarpetData(file):
    '''
    Getting Carpet data, assigning some column values as class attributes.
    '''
    data = []
    fop = open(file)
    for line in fop.readlines():
        sl = line.split()
        if sl[8] != '.0':
            continue
        c = Carpet()
        year = '20' + sl[0]
        day = sl[1]
        day.rjust(3 + len(day), '0')
        d = datetime.strptime(year + "-" + day, "%Y-%j").strftime("%d.%m.%y")
        c.date = d
        t = sl[3] + ':' + sl[4] + ':' + sl[5]
        c.time = t
        dt = datetime.strptime(d + ' ' + t, '%d.%m.%y %H:%M:%S')
        c.mjd = dt 
        if file == '2024_001-102_n.txt':
            c.RA = float(sl[13])
            c.DEC = float(sl[14])
            c.Ne = float(sl[17])
            c.mu = float(sl[8])
        else:
            c.RA = float(sl[12])
            c.DEC = float(sl[13])
            c.Ne = float(sl[16])
            c.mu = float(sl[7])
        data.append(c)
    fop.close()
    return data
        
def equ2ga(equ):
    '''
    Converts Equatorial to Galactic coordinates (J2000.0)
    Source: https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_coords.php

    Parameters
    ----------
    [ra, dec] : array
        Right Ascension and Declination in decimal degrees
    
    Returns
    -------
    [l, b] : array
        Galactic longitude and galactic latitude in decimal degrees
    '''
    ra = np.radians(equ[0])
    dec = np.radians(equ[1])

    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = np.radians(192.7666)
    pole_dec = np.radians(27.13)
    posangle = np.radians(122.93314)
    
    # North galactic pole (B1950)
    #pole_ra = radians(192.25)
    #pole_dec = radians(27.4)
    #posangle = radians(123.0-90.0)
    
    l=posangle-np.arctan((np.cos(dec)*np.sin(ra-pole_ra))/(np.sin(dec)*np.cos(pole_dec)-np.cos(dec)*np.sin(pole_dec)*np.cos(ra-pole_ra)))
    b=np.arcsin(np.sin(dec)*np.sin(pole_dec)+np.cos(dec)*np.cos(pole_dec)*np.cos(ra-pole_ra))
    
    return np.array([np.degrees(l), np.degrees(b)])

def GetLHAASOData():
    '''
    Getting LHAASO data from the catalogue.
    Downloaded at: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/271/25
    '''
    data = []
    fop = open('catalog.dat')
    n = ''
    for line in fop.readlines():
        sl = line.split('|')
        if sl[2] != '      ':
            lb = equ2ga([float(sl[2]), float(sl[3])])
            if abs(lb[1]) > 10:
                continue
        elif sl[2] == '      ':
            continue
        if sl[0] == n:
            continue
        l = LHAASO()
        l.name = sl[0]
        n = sl[0]
        l.detect = sl[1]
        l.RA  = float(sl[2])
        l.DEC = float(sl[3])
        l.unc95 = float(sl[4])
        data.append(l)
    fop.close()
    return data

def Angle(carpet_data, RA, DEC):
    '''
    Angle between the LHAASO object and a Carpet event.

    Parameters
    ----------
    carpet_data : array
        Array of Carpet events.
    RA : float
        Right Ascension of the LHAASO object.
    DEC : float
        Declination of the LHAASO object.

    Returns
    -------
    da : array
        Angle between the LHAASO object and the Carpet event.

    '''
    RA_b = RA * np.pi / 180 #LHAASO Right Ascension
    DEC_b = DEC * np.pi / 180 #LHAASO Declination
    ra = [o.RA * np.pi / 180 for o in carpet_data]
    dec = [o.DEC * np.pi / 180 for o in carpet_data]
    
    da = [] #delta alpha values
    for i in range(len(ra)):
        da.append(np.sin(DEC_b)*np.sin(dec[i]) + np.cos(DEC_b)*np.cos(dec[i])*np.cos(RA_b-ra[i]))
        
    return da

def Area(carpet_data, lhaaso, da, RA, DEC):
    '''
    Calcutaling signal and background.

    Parameters
    ----------
    carpet_data : array
        Array of Carpet events.
    lhaaso : array
        Array of LHAASO objects.
    da: array
        Angle between the LHAASO object and the Carpet event.
    RA : float
        Right Ascension of the LHAASO object.
    DEC : float
        Declination of the LHAASO object.

    Returns
    -------
    signal : float
        Signal.
    background : float
        Noise.

    '''
    DEC_b = DEC * np.pi / 180 #LHAASO Declination converted to rad
    decl = [o.DEC * np.pi / 180 for o in carpet_data] #Carpet Declination converted to rad
    
    signal = 0 #Signal event counter. Signal is calculated from the Carpet resolution circle around the LHAASO object.
    for i in range(len(da)):
        if da[i] >= DA_c:
            signal += 1
   
    background = 0 #Background event counter. Background is calculated from a strip, Declination (strip width) = signal circle diameter for all RA.
    for i in range(len(decl)):
        if abs(decl[i] - DEC_b / np.pi * 180) <= 4.7:
            background += 1

    for lh in lhaaso:
        dec = lh.DEC
        theta1 = 90 - dec - res
        theta2 = 90 - dec + res
        omega1 = 2 * np.pi * (1 - np.cos(theta1 * np.pi / 180))
        omega2 = 2 * np.pi * (1 - np.cos(theta2 * np.pi / 180))
        omega = abs(omega1 - omega2) * 180**2 / (np.pi)**2
        omega -= circle
    
    signal -= background #Subtracting signal from background
    background /= omega #Background density
    signal /= circle #Signal density
    
    return signal, background

def PlotAngle(carpet_data, da, files_dir):
    '''
    Plotting the angle between the LHAASO object and the Carpet event.

    Parameters
    ----------
    carpet_data : array
        Array of Carpet events.
    da : array
        Angle between the LHAASO object and the Carpet event.
    files_dir : str
        Name of the directory where the plots will be saved.
    '''
    t = [o.mjd for o in carpet_data] #datetime
    m = [o.mu for o in carpet_data] #muon counter
    ne = [o.Ne for o in carpet_data] #relativistic particles counter 
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    ax.set_xlabel('Date')
    ax.set_ylabel('Cos(d)')
    
    for i in range(len(da)):
        if da[i] >= 0.99:
            if m[i] == 0.0 and ne[i] > 35000:
                ax.plot(t[i], da[i], 'r*', ms=20, label=r'$\gamma$-candidate, $N_{\mu}$=0, Ne>35000')
            else:
                ax.plot(t[i], da[i], 'bo', label='Carpet events')
    
    ax.axhline(DA_c, color='k', label='Ang.res(4.7deg)')
          
    plt.savefig('/Users/darinamustakimova/Desktop/БНО практика/' + files_dir + '/' + 'Anglenew.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return

def PlotSN(carpet_data, RA, DEC, name_dir):
    '''
    Plotting the S/N.

    Parameters
    ----------
    carpet_data : array
        Array of Carpet events.
    RA : float
        Right Ascension of the LHAASO object.
    DEC : float
        Declination of the LHAASO object.
    name_dir : str
        Name of the directory where the plots will be saved.
    '''
    t = [o.mjd for o in carpet_data]
    time_s = list(set([m.strftime('%m.%y') for m in t]))
    time_s = sorted([datetime.strptime(d, '%m.%y') for d in time_s])
    
    sn = []
    for j in range(len(time_s)):
        idx = []
        for i in range(len(t)):
            if t[i].strftime('%m.%y') == time_s[j].strftime('%m.%y'):
                idx.append(i)
        da = Angle(carpet_data[idx[0]:idx[-1]], lh.RA, lh.DEC)
        c, p = Area(carpet_data[idx[0]:idx[-1]], da, lh.RA, lh.DEC)
        if p!=0: #checking if background is 0 so that we don't divide by 0
            sn.append(c / p)
        else:
            sn.append(0)
        
      
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    ax.set_xlabel('Date')
    ax.set_ylabel('S/N')
    #plt.scatter(time_s, sn)
    plt.step(time_s, sn)
    ax.axhline(1, color='k', linestyle='--')
    plt.title(name_dir)
    plt.savefig('/Users/darinamustakimova/Desktop/БНО практика/' + name_dir + '/' + 'SN.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    
    return

if __name__ == "__main__":
    files = ['2018_097-357_n.txt', '2019_001-365_n.txt', '2020_001-366_n.txt', \
             '2021_001-365_n.txt', '2022_001-365_n.txt', '2023_001-365_n.txt', \
             '2024_001-102_n.txt'] #files from 2018 to 2024
    lhaaso = GetLHAASOData()
    
    carpet = []
    for file in files: #looping through all the files
        c = GetCarpetData(file)
        carpet += c
        year = file.split('_')[0]
            
    DA = Angle(carpet, lh.RA, lh.DEC)
        
    namedir=lh.name.split(' ')[0]
    CreateDir(namedir)
    PlotAngle(carpet, DA, namedir)
    PlotSN(carpet, lh.RA, lh.DEC, namedir)
    
    
        
    
    
