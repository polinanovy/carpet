#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p-value for all LHAASO objects - final version
"""

import numpy as np
from datetime import datetime, timedelta
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import cmath
import random
from scipy.special import erfinv
import time

st = time.time()

RES = 4.7 #resolution of Carpet
DA_c = np.cos(4.7 * np.pi / 180)

longuitude = 42.684872 # Carpet
latitude = 43.273004 # Carpet
Phi = latitude*np.pi/180
sinB = np.sin(Phi)
cosB = np.cos(Phi)

N_SIM = 10**4 #Number of MC simulations

class LHAASO():
    def __init__(self):
        self.name = None
        self.RA  = None
        self.DEC = None
        self.unc95 = None
        self.detect = None #type of detector from LHAASO

class Carpet():
    def __init__(self):
        self.date = None
        self.RA  = None
        self.DEC = None
        self.Ne = None
        self.mu = None #n_mu
        self.b = None #galactic latitude
        self.theta = None #zenith angle
        self.phi = None

def GetLHAASOData():
    data = []
    fop = open('catalog.dat.txt')
    n = '' #string for objects names
    for line in fop.readlines():
        sl = line.split('|')
        if line.startswith("#"):
            continue
        if sl[2] != '      ':
            lb = equ2ga([float(sl[2]), float(sl[3])])
            if abs(lb[1]) > 10: #LHAASO objects from the Galactic plane
                continue
        elif sl[2] == '      ':
            continue
        #we take data only from one detector (with the less resolution) 
        #in the case there are two detectors for one object
        if sl[0] == n: 
            continue
        l = LHAASO()
        l.name = sl[0]
        n = sl[0]
        l.detect = sl[1]
        l.RA  = float(sl[2])
        l.DEC = float(sl[3])
        l.unc95 = float(sl[4])
        if l.DEC>5 and l.DEC<76:
            data.append(l)
    fop.close()
    return data

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
            if float(sl[idx]) != 0:
                continue
        c = Carpet()
        year = '20' + sl[0]
        day = sl[1]
        day.rjust(3 + len(day), '0')
        d = datetime.strptime(year + "-" + day, "%Y-%j").strftime("%d.%m.%y")
        t = sl[3] + ':' + sl[4] + ':' + sl[5]
        dt = datetime.strptime(d + ' ' + t, '%d.%m.%y %H:%M:%S')
        c.date = dt 
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

def Angle(carpet_data, RA, DEC):
    RA_b = RA * np.pi / 180
    DEC_b = DEC * np.pi / 180
    ra = [o.RA * np.pi / 180 for o in carpet_data]
    dec = [o.DEC * np.pi / 180 for o in carpet_data]
    
    da = []
    for i in range(len(ra)):
        da.append(np.sin(DEC_b)*np.sin(dec[i]) + np.cos(DEC_b)*np.cos(dec[i])*np.cos(RA_b-ra[i]))
        
    return da

def Angle_one(RA_carp, DEC_carp, lhaaso_obj):
    """
    Function for angle between the source and one carpet event
    """
    RA_obj = lhaaso_obj.RA
    DEC_obj = lhaaso_obj.DEC
    da = np.sin(DEC_obj)*np.sin(RA_carp) + np.cos(DEC_obj)*np.cos(DEC_carp)*np.cos(RA_obj-RA_carp)
    return da

def Events_match(carpet_data, lhaaso_obj, Photon=True):
    obj = lhaaso_obj.name
    RA_obj = lhaaso_obj.RA
    DEC_obj = lhaaso_obj.DEC
    Detector = lhaaso_obj.detect
    m = 0
    DA_carpet=[]
    
    DA_res = np.cos(4.7 * np.pi / 180) #delta alpha of the Carpet resolution
    
    DA_carpet=Angle(carpet_data, RA_obj, DEC_obj)
    for DA in DA_carpet:
            if DA>=DA_res:
                m+=1
    return m
'''
    
    if Photon:
        fop = open('Match_photon.txt', 'w')
    else:
        fop = open('Match_allcarpet.txt', 'w')
    fop.write('#LHAASO name      Detect.   RA   DEC   DateC       TimeC   CRA  CDEC   Ne \n')
    for carp in carpet_data:
        #c_carpet = SkyCoord(ra=carp.RA*u.degree, dec=carp.DEC*u.degree, frame='icrs')
        #c_obj = SkyCoord(ra=RA_obj*u.degree, dec=DEC_obj*u.degree, frame='icrs')
        #check separation between two objects and Carpet resolution
        #sep=c_obj.separation(c_carpet).degree
        da_obj=Angle(carp, lhaaso_obj) #delta alpha of the Carpet event
        if  da_obj>=DA_c:
            m += 1
            fop.write(obj + ' ' + Detector + ' ' + str(RA_obj) + ' ' + str(DEC_obj) +\
                      ' ' + str(carp.date) + ' ' + str(carp.RA) + ' ' + str(carp.DEC) + ' ' + str(carp.Ne) + '\n')
    fop.close()
    '''


def PlotAngle(carpet_data, da, obj):
    t = [o.date for o in carpet_data]
    m = [o.mu for o in carpet_data]
    
    fig, ax = plt.subplots()
    fig.set_size_inches(40,10)
    font = {'size'   : 30, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Date')
    ax.set_ylabel('Cos(d)')
    
    for i in range(len(da)):
        if da[i] >= 0.99:
            if m[i] == 0.0:
                ax.plot(t[i], da[i], 'r*', ms=20, label=r'$\gamma$-candidate, $N_{\mu}$=0, Ne>35000')
            else:
                ax.plot(t[i], da[i], 'bo', label='Carpet events')
    
    ax.axhline(DA_c, color='k', label='Ang.res(4.7deg)')
     
    plt.title(obj) 
    plt.savefig('Angle.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    return

def sidtime(date): #Sidereal time for AltAz coordinates
    y = date.year
    m = date.month
    d = date.day
    UTh = date.hour
    UTmin = date.minute
    UTsec = date.second
    d0 = 367*y - math.modf(7/4*(y + math.modf((m + 9)/12)[1] ))[1] + math.modf(275*m/9)[1] + d - 730531.5 + (UTh + UTmin/60. + UTsec/3600.)/24
    sdeg = math.fmod(280.46061837 + 360.98564736629*d0 + longuitude, 360)
    return sdeg #s in degrees

def equatorial(A, Z, s):
    #RA and DEC in degrees from azimuth and zenith in degrees
    A = A * np.pi / 180
    sinZ = math.sin(Z * np.pi / 180)
    sinA = math.sin(A)
    cosA = math.cos(A)
    cosZ = (1-sinZ**2)**0.5
    sinDelta = sinB*cosZ - cosB*sinZ*cosA  
    cosDelta = (1 - sinDelta**2)**0.5
    delta = math.asin(sinDelta)
    cosT = (cosZ*cosB + sinZ*sinB*cosA)/cosDelta
    if sinA > 0 :
       t = cmath.acos(cosT)
    else : 
       t = (math.pi + cmath.acos(-cosT))
    alpha = (s - t).real;
    Delta = 180.0*delta/math.pi
    if alpha > 0 :
       Alpha = 180.0*alpha/math.pi
    else: 
       Alpha = 180*(alpha + 2*math.pi)/math.pi
    return [Alpha, Delta]
   
def shuffleCarpet(carpet_data):
    #Create MC data according to S.Troitsky
    """
    Из каталога фотонных кандидатов берутся случайным образом зенитный угол от одного события
    и азимутальный угол от другого события. Приписывается случайное время прихода (звездное время,
    равномерное распределение от 0 до 24 ч). Восстанавливаются экваториальные координаты 
    случайного события(RA, DEC). Таким образом, составляется каталог случайных событий той же 
    длины, что реальное число фотонных кандидатов в данных.
    """
    zenith = [o.theta for o in carpet_data]
    azimuth = [o.phi for o in carpet_data]
    data = []
    for carp in carpet_data:
        c = Carpet()
        c.theta = random.choice(zenith)
        c.phi = random.choice(azimuth)
        c.date = carp.date + timedelta(hours=random.uniform(0, 24)) 
        stime = sidtime(c.date) * np.pi / 180
        alphadelta = equatorial(c.phi, c.theta, stime)
        c.RA = alphadelta[0]
        c.DEC = alphadelta[1]
        c.Ne = carp.Ne
        c.mu = carp.mu
        data.append(c)
    return data

if __name__ == "__main__":
    Photon = True
    files = ['2018_097-357_n.txt', '2019_001-365_n.txt', '2020_001-366_n.txt', \
             '2021_001-365_n.txt', '2022_001-365_n.txt', '2023_001-365_n.txt', \
             '2024_001-102_n.txt']
    lhaaso = GetLHAASOData()
    carpet = []
    fop = open('allmatches', 'w') 
    for file in files:
        c = GetCarpetData(file, Photon)
        carpet += c
    
    for lh in lhaaso:
        fop = open(lh.name + 'matches_pval.txt', 'w') 
        match_obs = Events_match(carpet, lh, Photon)
        if match_obs == 0:
            continue
        print("Number of matches with ", lh.name, " in real data = ", match_obs) 
        fop.write("Number of matches with " + lh.name + " in real data = " + str(match_obs) + '\n')  
    
        M = 0 # successful trials
        n = []
        for i in range(N_SIM):
            carpet_new = shuffleCarpet(carpet)
            match = Events_match(carpet_new, lh, Photon)
            print("Number of matches in", i, "simulated data = ", match)
            if match_obs <= match:
                M += 1
            n.append(match)
        pval = (M + 1)/(N_SIM + 1)
        print("p-value = ", pval)
        fop.write("p-value = " + str(pval) + '\n')
        xsigmas = erfinv(1-pval) * np.sqrt(2.)
        print("in sigmas this is =",xsigmas)
        fop.write("in sigmas this is =" + str(xsigmas) + '\n')
        n_expected = sum(n) / N_SIM
        print('n_expected = ', n_expected)
        fop.write('n_expected = ' + str(n_expected) + '\n')
        
        fop.close()


print("Time: %.03f s" % (time.time() - st))







