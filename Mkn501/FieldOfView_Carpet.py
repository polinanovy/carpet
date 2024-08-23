import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time
import astropy.units as u
import astropy.coordinates as coord

loc = EarthLocation(lat = 43.27*u.deg, lon = 42.6886*u.deg, height = 1700*u.m)
#time = Time('2024-03-07 07:44:57.38')
time = Time('2024-03-08 03:55:52')
az_deg = np.linspace( 0 ,360 , 100 )#np.linspace( -np.pi , np.pi , 2000 )
el_deg = 90-40#
direc = AltAz(location=loc, obstime=time, az=az_deg * u.deg, alt=el_deg * u.deg)
sky = SkyCoord(direc.transform_to(ICRS()))
                 
RA = sky.ra.wrap_at(180*u.degree)*0.01745329252
DEC = sky.dec* u.deg *0.01745329252


plt.subplot(projection="mollweide")
plt.grid(True)
#plt.title("Carpet-2 Field Of View \n at 14:13:30 ")
plt.xlabel('Right ascension')
plt.ylabel('Declination')

line1 = plt.plot(RA, DEC, linestyle = 'dashed', c='r', linewidth =  1.0)
line2 = plt.scatter(RA, DEC)

plt.legend()
plt.savefig('FieldOfView.png', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', format='png', transparent=False, bbox_inches=None, pad_inches='tight', metadata=None)
plt.show()
