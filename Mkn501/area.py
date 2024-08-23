#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt

dec = 39.76
res = 4.7
theta1 = 90 - dec - res
theta2 = 90 - dec + res

omega1 = 2 * np.pi * (1 - np.cos(theta1 * np.pi / 180))
omega2 = 2 * np.pi * (1 - np.cos(theta2 * np.pi / 180))

omega = abs(omega1 - omega2) * 180**2 / (np.pi)**2
print(omega)

circle = np.pi * res**2
print(circle)
omega -= circle



RA_b = 253.47 * np.pi / 180
DEC_b = 39.76 * np.pi / 180
DA_c = np.cos(4.7 * np.pi / 180)

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
    if file == '2024_001-102_n.txt':
        data = []
        fop = open(file)
        for line in fop.readlines():
            sl = line.split()
            #if sl[7] != '.0':
            #    continue
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
        fop.close()
    else:
        data = []
        fop = open(file)
        for line in fop.readlines():
            sl = line.split()
            #if sl[7] != '.0':
            #    continue
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

def Angle(carpet_data):
    ra = [o.RA * np.pi / 180 for o in carpet_data]
    dec = [o.DEC * np.pi / 180 for o in carpet_data]
    
    da = []
    for i in range(len(ra)):
        da.append(np.sin(DEC_b)*np.sin(dec[i]) + np.cos(DEC_b)*np.cos(dec[i])*np.cos(RA_b-ra[i]))
        
    return da

def Area(carpet_data, da):
    m = [o.mu for o in carpet_data]
    ne = [o.Ne for o in carpet_data]
    decl = [o.DEC for o in carpet_data]
    '''
    c = 0
    for i in range(len(da)):
        if da[i] >= DA_c:
            c += 1
   
    p = 0
    for i in range(len(decl)):
        if abs(decl[i] - DEC_b / np.pi * 180) <= 4.7:
            p += 1
    
    p -= c
    p /= omega
    c /= circle
            
    '''
    cr_c = 0; cb_c = 0 #Cчетчки красных и голубых точек для кружка радиусом 4.7 градуса
    for i in range(len(da)):
        if da[i] >= DA_c:
            if m[i] == 0.0 and ne[i] > 35000:
                cr_c += 1
            else:
                cb_c += 1
    
    cr_p = 0; cb_p = 0 #Cчетчки красных и голубых точек для полной полосы
    for i in range(len(decl)):
        if abs(decl[i] - DEC_b / np.pi * 180) <= 4.7:
            if m[i] == 0.0 and ne[i] > 35000:
                cr_p += 1
            else:
                cb_p += 1
               
    cr_p -= cr_c #Вычитаем из полосы кружок
    cb_p -= cb_c 
    
    cr_p = cr_p / omega
    cb_p = cb_p / omega
    cr_c = cr_c / circle
    cb_c = cb_c / circle
    
    
    
    return cr_p, cb_p, cr_c, cb_c

def Plot(carpet_data):
    #carpet_data = carpet
    t = [o.mjd for o in carpet_data]
    #time_s = sorted(t)
    time_s = list(set([m.strftime('%m.%y') for m in t]))
    time_s = sorted([datetime.strptime(d, '%m.%y') for d in time_s])
    #time_s = sorted(time_s)
    
    #time_s = sorted(list(set([m.strftime('%m.%y') for m in t])))
    #time_s = [datetime.strptime(d, '%m.%y') for d in time_s]
    #print(time_s)
    #counts_on = []; counts_off = []
    sn = []
    for j in range(len(time_s)):
        idx = []
        for i in range(len(t)):
            if t[i].strftime('%m.%y') == time_s[j].strftime('%m.%y'):
                idx.append(i)
        da = Angle(carpet_data[idx[0]:idx[-1]])
        #cr_p, cb_p, cr_c, cb_c = Area(carpet_data[idx[0]:idx[-1]], da)
        c, p = Area(carpet_data[idx[0]:idx[-1]], da)
        sn.append(c / p)
      
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel('S/N')
    plt.scatter(time_s, sn)
    #plt.step(time_s, sn)
    ax.axhline(1, color='k', linestyle='--')
    plt.savefig('SN.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    
    return

if __name__ == "__main__":
    files = ['2018_097-357_n.txt', '2019_001-365_n.txt', '2020_001-366_n.txt', \
             '2021_001-365_n.txt', '2022_001-365_n.txt', '2023_001-365_n.txt', \
             '2024_001-102_n.txt']
    carpet = []
    for file in files:
        c = GetCarpetData(file)
        carpet += c
    DA = Angle(carpet)
    cr_p, cb_p, cr_c, cb_c = Area(carpet, DA)
    print(cr_p, cb_p, cr_c, cb_c)
    #print(omega/circle)
    #Plot(carpet)