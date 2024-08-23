#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime
import os
import numpy as np
import ephem
import math
import cmath
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl

longuitude = 42.684872 # Carpet
latitude = 43.273004 # Carpet
Phi = latitude*np.pi/180
sinB = np.sin(Phi)
cosB = np.cos(Phi)

class Carpet():
    def __init__(self):
        self.date = None
        self.RA  = None
        self.DEC = None
        self.mu = None #n_mu
        self.b = None #galactic latitude
        self.theta = None #zenith angle
        self.phi = None
        
def GetCarpetData(file):
    data = []
    YEAR = file.split('-')[0]
    fop = open('Carpet2_TrgCygnus_AllEvents/20' + str(YEAR) + '/' + file)
    for line in fop.readlines():
        sl = line.split()
        c = Carpet()
        year = '20' + sl[0]
        day = sl[1]
        day.rjust(3 + len(day), '0')
        d = datetime.strptime(year + "-" + day, "%Y-%j").strftime("%d.%m.%y")
        if float(sl[3]) == 24:
            hours = '0'
        else:
            hours = sl[3]
        t = hours + ':' + sl[4] + ':' + sl[5]
        dt = datetime.strptime(d + ' ' + t, '%d.%m.%y %H:%M:%S')
        c.date = dt 
        c.RA = float(sl[11])
        c.DEC = float(sl[12])
        c.mu = float(sl[8])
        c.b = float(sl[13])
        c.phi = float(sl[9])
        c.theta = float(sl[10])
        data.append(c)
    fop.close()
    return data

def MoonPosition(carpet_data):
    '''
    Returns the azimuth and altitude of the Moon at the time of the Carpet event.

    Parameters
    ----------
    carpet_data : list
        List of Carpet events data. Each event should have datetime as an attribute (carpet_data.date).
    
    Returns
    -------
    moon_alt : float
        altitude of the Moon at the time of the Carpet event, deg

    moon_az : float
        azimuth of the Moon at the time of the Carpet event, deg
    '''
    observer = ephem.Observer()
    observer.lat = str(latitude)
    observer.lon = str(longuitude)
    observer.date = str(carpet_data.date)
    moon = ephem.Moon()
    moon.compute(observer)
    moon_alt=float(moon.alt)/np.pi*180
    moon_az=float(moon.az)/np.pi*180
    return moon_alt, moon_az #return in deg

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
    '''
    Converting AltAz to Equatorial coordinates

    Parameters
    ----------
    A : float
        Altitude in degrees

    Z : float
        Azimuth in degrees

    s : float
        Sidereal time in radians

    Returns
    -------
    [ra, dec] : array
        Right Ascension and Declination in decimal degrees

    '''
    A = A * np.pi / 180 #converting to rad
    sinZ = math.sin(Z * np.pi / 180) #converting to rad
    sinA = math.sin(A)
    cosA = math.cos(A)
    cosZ = (1-sinZ**2)**0.5
    sinDelta = sinB*cosZ - cosB*sinZ*cosA  
    cosDelta = (1 - sinDelta**2)**0.5
    delta = math.asin(sinDelta)
    cosT = (cosZ*cosB + sinZ*sinB*cosA)/cosDelta
    if sinA > 0 :
       t = cmath.acos(cosT) #t in radians
    else : 
       t = (math.pi + cmath.acos(-cosT)) #t in radians
    alpha = (s - t).real; 
    Delta = 180.0*delta/math.pi 
    if alpha > 0 :
       Alpha = 180.0*alpha/math.pi
    else: 
       Alpha = 180*(alpha + 2*math.pi)/math.pi
    return [Alpha, Delta]

def Angle(carpet_data, alt, phi):
    '''
    Returns the angle between the object and the Carpet event azimuthal coordinates

    Parameters
    ----------
    carpet_data : list
        List of Carpet events data. Each event should have theta and phi as attributes (o.theta, o.phi)
    alt : float
        Altitude of the object (e.g. the Moon), deg
    phi : float
        Azimuth of the object (e.g. the Moon), deg

    Returns
    -------
    psi : float
        angle between the object and the Carpet event, deg

    '''
    s = sidtime(carpet_data.date)
    alpha_obj, delta_obj = equatorial(phi, 90-alt, s * np.pi /180) #returns RA and DEC of the Moon in degrees
    #converting to radians
    alpha_obj = alpha_obj * np.pi /180 
    delta_obj = delta_obj * np.pi /180
    alpha_c = carpet_data.RA * np.pi /180
    delta_c = carpet_data.DEC * np.pi /180
    #calculating the angle between the Moon and the Carpet event
    psi_rad=np.arccos(np.sin(delta_obj)*np.sin(delta_c) + np.cos(delta_obj)*np.cos(delta_c)*np.cos(alpha_obj-alpha_c))
    psi=psi_rad*180/np.pi 
    return psi

def Binning(DA):
    '''
    Returns the number of events in each bin

    Parameters
    ----------
    DA : list
        List of angles between the object and the Carpet event

    Returns
    -------
    N : array
        2D array of the number of events in each bin
    N_omega : array
        2D array of the density of events in each bin
    
    '''
    bins = list(np.arange(0, 15*0.5, 0.5)) #bins from 0 to 15 deg with step 1 deg
    angles = []
    for i in range(len(bins)-1):
        a = []
        for angle in DA:
            if bins[i] <= abs(angle) < bins[i+1]: #assigning event angles to the bins
                a.append(angle)
        angles.append(a)

    psi_0 = 0.5 #step size, deg
    omega = []
    omega_0 = np.pi * (psi_0)**2 #area of the circle
    omega.append(omega_0)

    # Calculating N (the the number of events in each bin)
    N = [] 
    for i in range(len(angles)):
        N.append(len(angles[i]))

    # Calculating omega (the area of the annulus)
    for i in range(len(bins)-2):
        omega.append(np.pi*(bins[i+2]**2 - bins[i+1]**2)) #area of the annulus
    
    # Calculating N_omega (the event density in each bin) as N/omega*omega_0
    N_omega = [] 
    for i in range(len(omega)):
        N_omega.append(N[i]/omega[i]*omega_0)

    return N, N_omega

def Model_Gaussian(psi, sigma_psi):
    '''
    Returns the model of the Gaussian function

    Parameters
    ----------
    psi : array, deg
        Array of angles between the Moon and the Carpet event
    sigma_psi : float, deg
        Carpet resolution

    Returns
    -------
    result : array
        Deficit of cosmic rays (delta N / N_background)
    '''
    psi_moon = 0.26 #Moon angular radius
    result = -(psi_moon**2)/(2*sigma_psi**2) * np.exp(-(psi**2)/(2*sigma_psi**2))
    return result*100

if __name__ == "__main__":
    files19 = os.listdir('Carpet2_TrgCygnus_AllEvents/2019')
    files20 = os.listdir('Carpet2_TrgCygnus_AllEvents/2020')
    files21 = os.listdir('Carpet2_TrgCygnus_AllEvents/2021')
    files = files21 # files19 + files20 + files21   
    
    #N = []; N_m18 = []; N_p18 = []; N_m24 = []; N_p24 = []; N_m30 = []; N_p30 = []
    N = [0 for i in range(1,15)]; N_m18 = [0 for i in range(1,15)]; N_p18 = [0 for i in range(1,15)]
    N_m24 = [0 for i in range(1,15)]
    N_p24 = [0 for i in range(1,15)]; N_m30 = [0 for i in range(1,15)]; N_p30 = [0 for i in range(1,15)]
    for i in range(0, len(files), 10):
        if i+10 <= len(files):
            delta = 10
        else:
            delta = len(files) - i
        print(i, i+delta)
        carpet = []
        for file in files[i:i+delta]:
            c = GetCarpetData(file)
            carpet += c
        
        #Calculating angles for real and fake (background) Moon 
        #Fake Moons are calculated as real Moon azimuth ±18, ±24, ±30. Altitude stays the same
        DA = []; DA_p18 = []; DA_m18 = []; DA_p24 = []; DA_m24 = []#; DA_p30 = []; DA_m30 = []
        for carp in carpet:
            alt, az = MoonPosition(carp)
            DA.append(Angle(carp, alt, az))
            DA_p18.append(Angle(carp, alt, az+14))
            DA_m18.append(Angle(carp, alt, az-14))
            DA_p24.append(Angle(carp, alt, az+28))
            DA_m24.append(Angle(carp, alt, az-28))
            #DA_p30.append(Angle(carp, alt, az+30))
            #DA_m30.append(Angle(carp, alt, az-30))
        
        N_h, N_omega = Binning(DA)
        N_m18_h, N_omega_m18 = Binning(DA_m18)
        N_m24_h, N_omega_m24 = Binning(DA_m24)
        #N_m30_h, N_omega_m30 = Binning(DA_m30)
        N_p18_h, N_omega_p18 = Binning(DA_p18)
        N_p24_h, N_omega_p24 = Binning(DA_p24)
        #N_p30_h, N_omega_p30 = Binning(DA_p30)
        
        N = list(map(sum, zip(N, N_h)))
        N_m18 = list(map(sum, zip(N_m18, N_m18_h)))
        N_p18 = list(map(sum, zip(N_p18, N_p18_h)))
        N_m24 = list(map(sum, zip(N_m24, N_m24_h)))
        N_p24 = list(map(sum, zip(N_p24, N_p24_h)))
        #N_m30 = list(map(sum, zip(N_m30, N_m30_h)))
        #N_p30 = list(map(sum, zip(N_p30, N_p30_h)))
               
    #Calculating deficit of cosmic rays from the direction of the Moon relative to the background
    n = 4
    rel_dif = []
    sigma = []
    for i in range(len(N)):
        #Mean of each bin for fake (background) moons
        N_f_mean = (N_m18[i]+N_m24[i]+N_p18[i]+N_p24[i])/n
        print(N[i], N_f_mean)
        #Deficit
        rel_dif.append((N[i]-N_f_mean)/N_f_mean * 100)
        
        #Uncertainty
        sigma.append(N[i]/N_f_mean * np.sqrt(1/N[i] + 1/n/N_f_mean)*100)
    
    #space_ang = np.arange(1, 15, 1)
    space_ang = np.arange(0.5, 15*0.5, 0.5)
    popt, pcov = curve_fit(Model_Gaussian, space_ang, rel_dif)
    print(popt)
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    font = {'size'   : 20, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel('Angular distance from the Moon (degree)')
    ax.set_ylabel(r'Cosmic ray flux deficit ($\%$)')
    ax.plot(space_ang, rel_dif, 'ko')
    ax.plot(space_ang, Model_Gaussian(space_ang, popt), '-r', label = 'Gaussian fit, ang.res.='+str(round(popt[0],3))+' deg')
    plt.errorbar(space_ang, rel_dif, yerr = sigma, fmt = 'o', color = 'k')
    plt.legend(loc='lower right')
    plt.savefig('Deficit_2019_1.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()