#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 10:11:45 2024

@author: darinamustakimova
"""

import os
import subprocess
from datetime import datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path


'''
Plotting signal-to-noise for all LHAASO data, not just one object
'''

res = 4.7
circle = np.pi * res**2
DA_c = np.cos(4.7 * np.pi / 180)

class LHAASO():
    def __init__(self):
        self.name = None
        self.RA  = None
        self.DEC = None
        self.unc95 = None
        self.detect = None

class Carpet():
    def __init__(self):
        self.date = None
        self.time = None
        self.RA  = None
        self.DEC = None
        self.Ne = None
        self.mjd = None
        self.mu = None

def GetCarpetData(file):
    data = []
    fop = open(file)
    if file =='2024_001-102_n.txt':
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
            c.mjd = dt #julian.to_jd(dt, fmt='mjd')
            c.RA = float(sl[13])
            c.DEC = float(sl[14])
            c.Ne = float(sl[17])
            c.mu = float(sl[8])
            data.append(c)
    else: 
        for line in fop.readlines():
            sl = line.split()
            if sl[7] != '.0':
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
            c.mjd = dt #julian.to_jd(dt, fmt='mjd')
            c.RA = float(sl[12])
            c.DEC = float(sl[13])
            c.Ne = float(sl[16])
            c.mu = float(sl[7])
            data.append(c)
    fop.close()
    return data

        
'''def GetCarpetData(file):
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
    return data'''
        
def GetLHAASOData():
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

def equ2ga(equ):
    """
    Convert Equatorial to Galactic coordinates (J2000.0)
    
    Input: [ra,dec] in decimal degrees
    Returns: [l,b] in decimal degrees
    
    Source: 
    - https://www.atnf.csiro.au/people/Tobias.Westmeier/tools_coords.php
    
    """
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

def Angle(carpet_data, RA, DEC):
    RA_b = RA * np.pi / 180
    DEC_b = DEC * np.pi / 180
    ra = [o.RA * np.pi / 180 for o in carpet_data]
    dec = [o.DEC * np.pi / 180 for o in carpet_data]
    
    da = []
    for i in range(len(ra)):
        da.append(np.sin(DEC_b)*np.sin(dec[i]) + np.cos(DEC_b)*np.cos(dec[i])*np.cos(RA_b-ra[i]))
        
    return da

def Area(carpet_data, da, RA, DEC):
    DEC_b = DEC * np.pi / 180
    decl = [o.DEC for o in carpet_data]
    
    dec = DEC
    theta1 = 90 - dec - res
    theta2 = 90 - dec + res
    omega1 = 2 * np.pi * (1 - np.cos(theta1 * np.pi / 180))
    omega2 = 2 * np.pi * (1 - np.cos(theta2 * np.pi / 180))
    omega = abs(omega1 - omega2) * 180**2 / (np.pi)**2
    omega -= circle   
    
    c = 0 #Количество событий из кружка
    for i in range(len(da)):
        if da[i] >= DA_c:
            c += 1
   
    p = 0 #Количество событий из полосы
    for i in range(len(decl)):
        if abs(decl[i] - DEC_b / np.pi * 180) <= 4.7:
            p += 1
    
    p -= c
    p /= omega
    c /= circle
    
    return c, p

def PlotAngle(carpet_data, da, name_dir, name):
    t = [o.mjd for o in carpet_data]
    m = [o.mu for o in carpet_data]
    ne = [o.Ne for o in carpet_data]
    
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
          
    plt.savefig('/Users/darinamustakimova/Desktop/БНО практика/' + name_dir + '/' + 'Anglenew.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return

def PlotSN(carpet_data, RA, DEC, name_dir, name):
    t = [o.mjd for o in carpet_data]
    time_s = list(set([m.strftime('%m.%y') for m in t]))
    time_s = sorted([datetime.strptime(d, '%m.%y') for d in time_s])
    
    sn = []
    cc=[]; cp=[]
    for j in range(len(time_s)):
        idx = []
        for i in range(len(t)):
            if t[i].strftime('%m.%y') == time_s[j].strftime('%m.%y'):
                idx.append(i)
        da = Angle(carpet_data[idx[0]:idx[-1]], RA, DEC)
        c, p = Area(carpet_data[idx[0]:idx[-1]], da, RA, DEC)
        cc.append(c); cp.append(p)
        if p!=0: 
            sn.append(c / p)
        else:
            sn.append(0)
    sn_full=sum(cc)/sum(cp) #общий signal to noise
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    ax.set_xlabel('Date')
    ax.set_ylabel('S/N')
    #plt.scatter(time_s, sn)
    plt.step(time_s, sn)
    ax.axhline(1, color='k', linestyle='--')
    plt.title(name + '(Full SN =' + str(sn_full)+ ')')
    plt.savefig('/Users/darinamustakimova/Desktop/БНО практика/' + name_dir + '/' + 'SN.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

def PlotSN_full(carpet_data, lhaaso_data):
    t = [o.mjd for o in carpet_data]
    time_s = list(set([m.strftime('%m.%y') for m in t]))
    time_s = sorted([datetime.strptime(d, '%m.%y') for d in time_s])
    sn = []
    for j in range(len(time_s)):
        idx = []
        for i in range(len(t)):
            if t[i].strftime('%m.%y') == time_s[j].strftime('%m.%y'):
                idx.append(i)
        cc=[]; cp=[] #общие счетчики для всех кружков и всех полосок за месяц
        for lh in lhaaso_data:
            da = Angle(carpet_data[idx[0]:idx[-1]], lh.RA, lh.DEC)
            c, p = Area(carpet_data[idx[0]:idx[-1]], da, lh.RA, lh.DEC)
            cc.append(c); cp.append(p)
        sn.append(sum(cc)/sum(cp))
        
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    ax.set_xlabel('Date')
    ax.set_ylabel('S/N')
    #plt.scatter(time_s, sn)
    plt.step(time_s, sn)
    ax.axhline(1, color='k', linestyle='--')
    #plt.title(name_dir)
    plt.savefig('SN_full.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
            
            
if __name__ == "__main__":
    files = ['2018_097-357_n.txt', '2019_001-365_n.txt', '2020_001-366_n.txt', \
             '2021_001-365_n.txt', '2022_001-365_n.txt', '2023_001-365_n.txt', \
             '2024_001-102_n.txt']
    lhaaso = GetLHAASOData()
    carpet = []
    for file in files:
        c = GetCarpetData(file)
        carpet += c
    
    
    for lh in lhaaso: 
        DA = Angle(carpet, lh.RA, lh.DEC)
        
        namedir=lh.name.split(' ')[0]
        CreateDir(namedir)
        PlotAngle(carpet, DA, namedir, lh.name)
        PlotSN(carpet, lh.RA, lh.DEC, namedir, lh.name)
        print('done')
        
        
    

        
        
    
    
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            