from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import AltAz

observing_location = EarthLocation(lat='12.9716d', lon='77.5946d')  
observing_time = Time('2018-06-05 14:40:00')  
aa = AltAz(location=observing_location, obstime=observing_time)

coord = SkyCoord('23.39056d', '58.8d')
coord.transform_to(aa)
print aa