#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy import units as u
from astropy.coordinates import SkyCoord
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import julian

class Neutrino():
    def __init__(self):
        self.date = None
        self.rev = None
        self.time = None
        self.event = None
        self.type = None
        self.RA  = None
        self.DEC = None
        self.unc90 = None
        self.unc50 = None
        self.En    = None        

class Carpet():
    def __init__(self):
        self.date = None
        self.time = None
        self.RA  = None
        self.DEC = None
        self.Ne = None

def GetCarpetData(file):
    data = []
    fop = open('2023_001-365_n.txt')
    for line in fop.readlines():
        sl = line.split()
        if sl[7] != '.0':
            continue
        c = Carpet()
        year = '20' + sl[0]
        day = sl[1]
        day.rjust(3 + len(day), '0')
        c.date = datetime.strptime(year + "-" + day, "%Y-%j").strftime("%d.%m.%y")
        c.time = sl[3] + ':' + sl[4] + ':' + sl[5]
        c.RA = float(sl[12])
        c.DEC = float(sl[13])
        c.Ne = float(sl[16])
        data.append(c)
    fop.close()
    return data

def GetNeutrinoData(Cascade=True):
    data = []
    if Cascade:
        fop = open('IceCubeCASCADE.csv')
        for line in fop.readlines():
            if line.startswith(";"):
                continue
            sl = line.split(";")
            n = Neutrino()
            n.event = sl[1]
            n.type = sl[6]
            n.date = sl[2]
            n.time = sl[3]
            n.RA  = float(sl[7])
            n.DEC = float(sl[8])
            n.unc90 = float(sl[9])
            n.unc50 = float(sl[10])
            n.En    = float(sl[11])   
            data.append(n)
        fop.close()
    else:
        fop = open('IceCubeGoldBronze.csv')
        for line in fop.readlines():
            if line.startswith(";"):
                continue
            sl = line.split(";")
            if  float(sl[2]) == 0:
                continue
            n = Neutrino()
            n.event = sl[1]
            n.type = sl[5]
            n.date = sl[3]
            n.time = sl[4]
            n.RA  = float(sl[6])
            n.DEC = float(sl[7])
            n.unc90 = float(sl[8]) / 60
            n.unc50 = float(sl[9]) / 60
            n.En    = float(sl[10])   
            data.append(n)
        fop.close()
    return data

def Plot(carpet_data, data_neutrino):  
    ra = [o.RA for o in carpet_data]
    dec = [o.DEC for o in carpet_data]
    
    date_c = [o.date for o in carpet_data]
    time_c = [o.time for o in carpet_data]
    dt_c = [datetime.strptime(date_c[i] + ' ' + time_c[i], '%d.%m.%y %H:%M:%S') for i in range(len(carpet_data))]
    mjd_c = [julian.to_jd(dt_c[i], fmt='mjd') for i in range(len(dt_c))]
    
    date_n = [o.date for o in data_neutrino]
    time_n = [o.time for o in data_neutrino]
    dt_n = [datetime.strptime(date_n[i] + ' ' + time_n[i], '%y/%m/%d %H:%M:%S.%f') for i in range(len(data_neutrino))]
    mjd_n = [julian.to_jd(dt_n[i], fmt='mjd') for i in range(len(dt_n))]
    
    plt.figure(figsize=(18.6,11.5))
    plt.subplot(111)

    m = Basemap(projection='moll',lon_0=0, resolution='c')

    m.drawparallels([-80,-60,-30,0,30,60,80],labels=[False,True,True,False],linewidth=0.5,dashes=[1, 3],labelstyle='+/-',fmt='%d')
    m.drawmeridians(np.arange(-180.,181.,60.),labels=[True,False,False,True],linewidth=0.5,dashes=[1, 4],labelstyle='+/-',fmt='%d')

    ax = plt.gca()
    ax.xaxis.set_label_coords(0.5, -0.05)
    
    r5s_rads = list(map(lambda x: x.unc90, data_neutrino))
    
    r5s_ras = list(map(lambda x: x.RA, data_neutrino))
    r5s_decs = list(map(lambda x: x.DEC, data_neutrino))
    phi = list(np.linspace(0, 2.*np.pi, 72))
    
    for c in range(len(mjd_c)):
        for n in range(len(mjd_n)):
            if mjd_n[n]-1 <= mjd_c[c] <= mjd_n[n]+1:
                r = r5s_rads[n]
                c_neut = SkyCoord(ra=r5s_ras[n]*u.degree, dec=r5s_decs[n]*u.degree, frame='icrs')
                c_source = SkyCoord(ra=ra[c]*u.degree, dec=dec[c]*u.degree, frame='icrs')
                sep = c_source.separation(c_neut).degree
                if sep < r:
                    print(sep)
                    m.plot(ra[c], dec[c], 'bo', latlon=True, markersize=5, alpha=0.9)
                    
                    i = n
                    r = r5s_rads[i]
                    x = r5s_ras[i] + r*np.cos(phi)
                    y = r5s_decs[i] + r*np.sin(phi)
                    
                    xx1 = []
                    xx2 = []
                    yy1 = []
                    yy2 = []
                    
                    for l in range(len(phi)):
                        if x[l] < 0:
                            x[l] += 360
                        if x[l] >= 180:
                            xx1.append(x[l])
                            yy1.append(y[l])
                        else:
                            xx2.append(x[l])
                            yy2.append(y[l])

                    m.plot(xx1, yy1, color="g", latlon=True, zorder=3)
                    m.plot(xx2, yy2, color="g", latlon=True, zorder=3)
    
    plt.savefig('sky_mollweide.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()
    
def Events_match(carpet_data, data_neutrino, Cascade=True):
    date_c = [o.date for o in carpet_data]
    time_c = [o.time for o in carpet_data]
    dt_c = [datetime.strptime(date_c[i] + ' ' + time_c[i], '%d.%m.%y %H:%M:%S') for i in range(len(carpet_data))]
    mjd_c = [julian.to_jd(dt_c[i], fmt='mjd') for i in range(len(dt_c))]
    
    date_n = [o.date for o in data_neutrino]
    time_n = [o.time for o in data_neutrino]
    dt_n = [datetime.strptime(date_n[i] + ' ' + time_n[i], '%y/%m/%d %H:%M:%S.%f') for i in range(len(data_neutrino))]
    mjd_n = [julian.to_jd(dt_n[i], fmt='mjd') for i in range(len(dt_n))]
    
    r5s_ras = list(map(lambda x: x.RA, data_neutrino))
    r5s_decs = list(map(lambda x: x.DEC, data_neutrino))
    r5s_rads = list(map(lambda x: x.unc90, data_neutrino))
    
    ras = list(map(lambda x: x.RA, carpet_data))
    decs = list(map(lambda x: x.DEC, carpet_data))
    Ne = list(map(lambda x: x.Ne, carpet_data))
    
    name_n = list(map(lambda x: x.event, data_neutrino))
    energy_n = list(map(lambda x: x.En, data_neutrino))
    type_n = list(map(lambda x: x.type, data_neutrino))
    
    if Cascade:
        fop = open('match_cascade.txt', 'w')
        fop.write('#RunEventNum    Date     Time      NeutrinoRA  NeutrinoDEC Energy DateC TimeC   CRA  CDEC   Ne \n')
        for c in range(len(mjd_c)):
            for n in range(len(mjd_n)):
                if mjd_n[n]-1 <= mjd_c[c] <= mjd_n[n]+1:
                    r = r5s_rads[n]
                    c_neut = SkyCoord(ra=r5s_ras[n]*u.degree, dec=r5s_decs[n]*u.degree, frame='icrs')
                    c_source = SkyCoord(ra=ras[c]*u.degree, dec=decs[c]*u.degree, frame='icrs')
                    sep = c_source.separation(c_neut).degree
                    if sep < r:
                        #print('match')
                        fop.write(name_n[n] + ' ' + date_n[n] + ' ' + time_n[n] + ' ' + str(r5s_ras[n]) + ' ' + str(r5s_decs[n]) + ' ' + str(energy_n[n] ) + ' ' + date_c[c] + ' ' + time_c[c] + ' ' + str(ras[c]) + ' ' + str(decs[c]) + ' ' + str(Ne[c]) + '\n')
        fop.close()
    else:
        fop = open('match_BG.txt', 'w')
        fop.write('#RunEventNum    Date     Time        Type NeutrinoRA  NeutrinoDEC Energy DateC TimeC   CRA  CDEC   Ne \n')
        for c in range(len(mjd_c)):
            for n in range(len(mjd_n)):
                if mjd_n[n]-1 <= mjd_c[c] <= mjd_n[n]+1:
                    r = r5s_rads[n]
                    c_neut = SkyCoord(ra=r5s_ras[n]*u.degree, dec=r5s_decs[n]*u.degree, frame='icrs')
                    c_source = SkyCoord(ra=ras[c]*u.degree, dec=decs[c]*u.degree, frame='icrs')
                    sep = c_source.separation(c_neut).degree
                    if sep < r:
                        #print(sep)
                        #print(r)
                        #print(mjd_n[n])
                        #print(c_l.separation(c_source).degree)
                        #print('match')
                        fop.write(name_n[n] + ' ' + date_n[n] + ' ' + time_n[n] + ' ' + type_n[n] + ' ' + str(r5s_ras[n]) + ' ' + str(r5s_decs[n]) + ' ' + str(energy_n[n] ) + ' ' + date_c[c] + ' ' + time_c[c] + ' ' + str(ras[c]) + ' ' + str(decs[c]) + ' ' + str(Ne[c]) + '\n')
        fop.close()
        
    return


if __name__ == "__main__":
    carpet23 = GetCarpetData('2023_001-365_n.txt')
    carpet24 = GetCarpetData('Car24_alert')
    
    neutrino_cascades = GetNeutrinoData(Cascade=True)
    neutrino_gb = GetNeutrinoData(Cascade=False)
    
    Plot(carpet23, neutrino_cascades)
    Events_match(carpet23, neutrino_cascades, Cascade=True)
    Events_match(carpet23, neutrino_gb, Cascade=False)
