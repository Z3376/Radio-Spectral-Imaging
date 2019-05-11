#!/usr/bin/python

import sys
from os.path import getsize
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

c			=	299792458
pi			=	3.14159265359
JD2000		=	2451558
d2r			=	pi/180

Lat			=	13.6029845
#Long		=	77.4279978
Long		=	77.0

wvlngth		=	input('Wavelength:')


Baseline	=	np.zeros((7,7))
Baseline[0,1]	=	82.83212
Baseline[0,2]	=	79.15458457
Baseline[0,3]	=	105.3456913
Baseline[0,4]	=	49.92998966
Baseline[0,5]	=	60.70695213
Baseline[0,6]	=	66.22207455
Baseline[1,2]	=	38.82850201
Baseline[1,3]	=	111.9512717
Baseline[1,4]	=	48.82121474
Baseline[1,5]	=	59.0537879
Baseline[1,6]	=	18.26778167
Baseline[2,3]	=	73.62784163
Baseline[2,4]	=	36.22617317
Baseline[2,5]	=	28.29191314
Baseline[2,6]	=	43.69229597
Baseline[3,4]	=	80.62542858
Baseline[3,5]	=	57.73338988
Baseline[3,6]	=	111.9270635
Baseline[4,5]	=	23.97696469
Baseline[4,6]	=	38.32823514
Baseline[5,6]	=	55.17263119

for i in range(len(Baseline)):
	for j in range(len(Baseline[0])):
		Baseline[j][i]	=	-1*Baseline[i][j]

if(len(sys.argv)==5):
	if(sys.argv[4]=='cas'):
		RA			=	23.39056
		dec			=	58.8
	elif(sys.argv[4]=='cyg'):
		RA			=	19.9912101667
		dec			=	40.7339155556
	elif(sys.argv[4]=='crab'):
		RA			=	5.57556
		dec			=	22.0145
else:
	RA			=	input('RA(in hrs): ')
	dec			=	input('dec(in deg): ')

file_name	=	sys.argv[1]
time_yr		=	float(file_name[-23:-19])
time_month	=	float(file_name[-19:-17])
time_day	=	float(file_name[-17:-15])
time_hr		=	float(file_name[-14:-12])-5
time_min	=	float(file_name[-12:-10])-30
time_s		=	float(file_name[-10:-8])

if(time_min<0):
	time_min	+=	60
	time_hr		-=	1
if(time_hr<0):
	time_hr		+=	24
	time_day	-=	1
if(time_day<1):
	if(time_month==5 or time_month==7 or time_month==10 or time_month==12):
		time_day	+=	30
	elif(time_month==1 or time_month==2 or time_month==4 or time_month==6 or time_month==8 or time_month==9 or time_month==11):
		time_day	+=	31
	elif(time_month==3):
		if(time_yr%4==0):
			time_day	+=	29
		else:
			time_day	+=	28
	time_month	-=	1
if(time_month<1):
	time_month	+=	12
	time_yr		-=	1


time_yr		=	2018
time_month	=	11
time_day	=	4
time_hr		=	13
time_min	=	30
time_s		=	0

t1			=	str(int(time_yr))+"-"+str(int(time_month))+"-"+str(int(time_day))+"T"+str(int(time_hr))+":"+str(int(time_min))+":"+str(int(time_s))
#print t1
t2			=	Time(t1,format="isot")


d_2000		=	jd-JD2000

das_id1		=	int(sys.argv[2])
das_id2		=	int(sys.argv[3])

baseline	=	Baseline[das_id1-1][das_id2-1]

file		=	open(file_name,"rb")
size		=	getsize(file_name)
i			=	0
gps			=	[]

g			=	np.empty(2)
q			=	0
while(i<size):
	file.read(26)
	q	+=	1
	g[q%2]	=	int(str(file.read(2)).encode('hex'),16)
	if(g[q%2]-g[(q+1)%2]!=0):
		gps.append(g[q%2])
	file.read(4)
	i		+=	32

file.close()

l			=	len(gps)

d_inc		=	[d_2000 for i in range(l)]

UT			=	time_hr+time_min/60+time_s/3600

UT_inc		=	[UT for i in range(l)]

LST			=	np.zeros(l)
ha			=	[0 for i in range(l)]
Alt			=	np.zeros(l)
t			=	np.zeros(l)

for i in range(l):
	d_inc[i]	+=	float(gps[i])*1.15741*10**(-5)
	UT_inc[i]	+=	float(gps[i])/3600
	LST[i]		=	(100.46+0.985647*d_inc[i]+Long+15*UT_inc[i])%360
	ha[i]		=	LST[i]-RA*15
	Alt[i]		=	np.arcsin(np.sin(dec*d2r)*np.sin(Lat*d2r)+np.cos(dec*d2r)*np.cos(ha[i]*d2r)*np.cos(Lat*d2r))
	t[i]		=	1200*wvlngth/pi/baseline/np.cos(Alt[i])

Alt_d	= Alt*(180/pi)

print min(t)


#g_delay_t.write(file_name[:-8]+'_g_delay.txt',format='ascii.no_header',overwrite='true')
plt.figure()
plt.plot(UT_inc,(LST/15))
plt.xlabel('Time(hrs)')
plt.ylabel('LST(hrs)')
plt.figure()
plt.plot(UT_inc,Alt_d)
plt.xlabel('Time(hrs)')
plt.ylabel('Alt(deg)')
plt.show()
