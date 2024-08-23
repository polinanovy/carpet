#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime, timedelta
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import julian
from scipy.optimize import curve_fit

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

def Plot(carpet_data, da):
    t = [o.mjd for o in carpet_data]
    m = [o.mu for o in carpet_data]
    ne = [o.Ne for o in carpet_data]
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel('Cos(d)')
    
    for i in range(len(da)):
        if da[i] >= 0.99:
            if m[i] == 0.0 and ne[i] > 35000:
                ax.plot(t[i], da[i], 'r*', ms=20, label=r'$\gamma$-candidate, $N_{\mu}$=0, Ne>35000')
            else:
                ax.plot(t[i], da[i], 'bo', label='Carpet events')
    
    ax.axhline(DA_c, color='k', label='Ang.res(4.7deg)')
    
    #plt.legend(loc=4)        
    plt.savefig('Angle.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return

def Model(x, a, b):
    return a*x + b

def Background(carpet_data, da, Sr, Sb):
    t = [o.mjd for o in carpet_data]
    m = [o.mu for o in carpet_data]
    ne = [o.Ne for o in carpet_data]
    
    time = []; counts_b = []; counts_r = []; flag = []
    for i in range(len(da)):
        if da[i] >= DA_c:
            if m[i] == 0.0 and ne[i] > 35000:
                flag.append('r')
            else:
                flag.append('b')
            time.append(t[i])
    
    time_s = sorted(list(set([m.strftime('%m.%y') for m in time])))
    time_s = [datetime.strptime(d, '%m.%y') for d in time_s]
    for j in range(len(time_s)):
        cr = 0; cb = 0
        for i in range(len(time)):
            if time[i].strftime('%m.%y') == time_s[j].strftime('%m.%y'):
                if flag[i] == 'r':
                    cr += 1
                else:
                    cb += 1
        counts_r.append(cr); counts_b.append(cb)
        print(cr)
        print(cb)
        print(time_s[j])
    
    #popt_r, pcov_r = curve_fit(Model, time_s, counts_r)    
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel('Counts')
    plt.scatter(time_s, counts_r, color='r')#, label=Sr)
    #plt.plot(time_s, Model(counts_r, *popt_r), 'r-')
    #plt.legend(loc=1)
    plt.savefig(Sr + '.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel('Counts')
    plt.scatter(time_s, counts_b, color='b')#, label=Sb)
    #plt.legend(loc=1)
    plt.savefig(Sb + '.png',bbox_inches='tight')
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
        da = Angle(c)
        Sr = file.split('_')[0] + '_red'
        Sb = file.split('_')[0] + '_blue'
        #Background(c, da, Sr, Sb)
    
    DA = Angle(carpet)
    #Plot(c, da)
    
    Background(carpet, DA, 'r_full', 'b_full')