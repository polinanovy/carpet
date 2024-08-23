#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt

RA_cc = 307.18
DEC_cc = 41.31

RES = 4.7 #resolution of Carpet
DA_c = np.cos(4.7 * np.pi / 180)
circle = np.pi * RES**2

class Carpet():
    def __init__(self):
        self.date = None
        self.year = None
        self.day = None
        self.msec = None
        self.RA  = None
        self.DEC = None
        self.Ne = None
        self.mu = None #n_mu
        self.b = None #galactic latitude
        self.theta = None #zenith angle
        self.phi = None
        self.da = None
        
def equ2ga(equ):
    """
    Convert Equatorial to Galactic coordinates (J2000.0)
    Input: [ra,dec] in decimal degrees
    Returns: [l,b] in decimal degrees
    Source: https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_coords.php
    """
    ra = np.radians(equ[0])
    dec = np.radians(equ[1])

    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = np.radians(192.7666)
    pole_dec = np.radians(27.13)
    posangle = np.radians(122.93314)
    
    l=posangle-np.arctan((np.cos(dec)*np.sin(ra-pole_ra))/(np.sin(dec)*np.cos(pole_dec)-np.cos(dec)*np.sin(pole_dec)*np.cos(ra-pole_ra)))
    b=np.arcsin(np.sin(dec)*np.sin(pole_dec)+np.cos(dec)*np.cos(pole_dec)*np.cos(ra-pole_ra))
    
    return np.array([np.degrees(l), np.degrees(b)])
        
def GetCarpetData(file, Photon=True):
    if file == '2024_001-102_n.txt': #for 2024 data we have different coloumns in file
        idx = 8 #index of coloumn with n_mu in different files
    else:
        idx = 7
    data = []
    fop = open(file)
    for line in fop.readlines():
        sl = line.split()
        if Photon: #Read only photon-like candidates
            if ((float(sl[idx]) + 0.1) / float(sl[idx+9])) > 10**(-5.90688): #float(sl[idx]) != 0:
                continue
        c = Carpet()
        year = '20' + sl[0]
        c.year = float(year)
        day = sl[1]
        c.day = float(day)
        day.rjust(3 + len(day), '0')
        d = datetime.strptime(year + "-" + day, "%Y-%j").strftime("%d.%m.%y")
        t = sl[3] + ':' + sl[4] + ':' + sl[5]
        dt = datetime.strptime(d + ' ' + t, '%d.%m.%y %H:%M:%S')
        c.date = dt 
        c.msec = float(sl[2])
        c.RA = float(sl[idx+5])
        c.DEC = float(sl[idx+6])
        c.Ne = float(sl[idx+9])
        c.mu = float(sl[idx])
        lb = equ2ga([float(sl[idx+5]), float(sl[idx+6])])
        c.b = lb[1]
        c.phi = float(sl[idx+3])
        c.theta = float(sl[idx+4])
        data.append(c)
    fop.close()
    return data

def Angle(carpet_data):
    RA_b = RA_cc * np.pi / 180
    DEC_b = DEC_cc * np.pi / 180
    ra = [o.RA * np.pi / 180 for o in carpet_data]
    dec = [o.DEC * np.pi / 180 for o in carpet_data]
    da = []
    for i in range(len(ra)):
        da.append(np.sin(DEC_b)*np.sin(dec[i]) + np.cos(DEC_b)*np.cos(dec[i])*np.cos(RA_b-ra[i]))
    return da

def Events_match(carpet_data, Photon):
    DA = Angle(carpet_data)
    m = 0
    if Photon:
        fop = open('CC_match_photons.txt', 'w')
    else:
        fop = open('CC_match.txt', 'w')
    years = [o.year for o in carpet_data]
    days = [o.day for o in carpet_data]
    msecs = [o.msec for o in carpet_data]
    Nes = [o.Ne for o in carpet_data]
    nmus = [o.mu for o in carpet_data]
    ras = [o.RA for o in carpet_data]
    decs = [o.DEC for o in carpet_data]
    fop.write('#Year  Day  mSec  Ne  n_mu  RA  DEC\n')
    for i in range(len(DA)):
        if DA[i] >= DA_c:
            m += 1
            fop.write(str(years[i]) + ' ' + str(days[i]) + ' ' + str(msecs[i]) +\
                      ' ' + str(Nes[i]) + ' ' + str(nmus[i]) + ' ' +\
                      str(ras[i]) + ' ' + str(decs[i]) + '\n')
    fop.close()
    return m

def PlotAngle(carpet_data, Photon=True):
    da = Angle(carpet_data)
    t = [o.date for o in carpet_data]
    m = [o.mu for o in carpet_data]
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel(r'Cos($\Delta \alpha$)')
    
    if Photon:
        for i in range(len(da)):
            if da[i] >= 0.99:
                ax.plot(t[i], da[i], 'bo')
        filename = 'CC_dalpha_photons'
        title = r'$\Delta \alpha$ between Carpet photon-like events and Cygnus Cocoon'
    else:
        for i in range(len(da)):
            if da[i] >= 0.99:
                if m[i] == 0.0:
                    ax.plot(t[i], da[i], 'r*', ms=20)
                else:
                    ax.plot(t[i], da[i], 'bo')
        filename = 'CC_dalpha'  
        title = r'$\Delta \alpha$ between Carpet events and Cygnus Cocoon'
        
    ax.axhline(DA_c, color='k')
    plt.title(title) 
    plt.savefig(filename + '.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return

def Area(carpet_data):
    da = Angle(carpet_data)
    DEC_b = DEC_cc * np.pi / 180
    decs = [o.DEC for o in carpet_data]
    
    theta1 = 90 - DEC_cc - RES
    theta2 = 90 - DEC_cc + RES
    omega1 = 2 * np.pi * (1 - np.cos(theta1 * np.pi / 180))
    omega2 = 2 * np.pi * (1 - np.cos(theta2 * np.pi / 180))
    omega = abs(omega1 - omega2) * 180**2 / (np.pi)**2
    omega -= circle
    
    c = 0 #Количество событий из кружка
    for i in range(len(da)):
        if da[i] >= DA_c:
            c += 1
    p = 0 #Количество событий из полосы
    for i in range(len(decs)):
        if abs(decs[i] - DEC_b / np.pi * 180) <= 4.7:
            p += 1
    p -= c
    p /= omega
    c /= circle
    return c, p

def PlotSN(carpet_data, Photon):
    t = [o.date for o in carpet_data]
    time_s = list(set([m.strftime('%m.%y') for m in t]))
    time_s = sorted([datetime.strptime(d, '%m.%y') for d in time_s])
    
    sn = []
    cc = []; cp = []
    for j in range(len(time_s)):
        idx = []
        for i in range(len(t)):
            if t[i].strftime('%m.%y') == time_s[j].strftime('%m.%y'):
                idx.append(i)
        c, p = Area(carpet_data[idx[0]:idx[-1]])
        cc.append(c); cp.append(p)
        if p != 0:
            sn.append(c / p)
        else:
            sn.append(0)
    sn_full = sum(cc) / sum(cp)
        
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel('S/N')
    plt.scatter(time_s, sn)
    ax.axhline(1, color='k', linestyle='--')
    plt.title('(Full SN = ' + str(sn_full) + ')') 
    if Photon:
        plt.savefig('CC_SN_photons.png',bbox_inches='tight')
    else:
        plt.savefig('CC_SN.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return

if __name__ == "__main__":
    Photon = True
    files = ['2018_097-357_n.txt', '2019_001-365_n.txt', '2020_001-366_n.txt', \
             '2021_001-365_n.txt', '2022_001-365_n.txt', '2023_001-365_n.txt', \
             '2024_001-102_n.txt']
    carpet_obs = []
    for file in files:
        c = GetCarpetData(file, Photon)
        carpet_obs += c
    match = Events_match(carpet_obs, Photon)
    print(match)
    PlotAngle(carpet_obs, Photon)
    PlotSN(carpet_obs, Photon)