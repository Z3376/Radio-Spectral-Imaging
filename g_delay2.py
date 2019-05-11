#!/usr/bin/python

import sys
from os.path import getsize
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

c			=	299792458
JD2000		=	2451558

Lat			=	13.6029845
Long		=	77.4279978

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

das_id_1	=	int(sys.argv[2])
das_id_2	=	int(sys.argv[3])

file_name	=	sys.argv[1]
time_yr		=	file_name[-23:-19]
time_month	=	file_name[-19:-17]
time_day	=	file_name[-17:-15]
time_hr		=	file_name[-14:-12]
time_min	=	file_name[-12:-10]
time_s		=	file_name[-10:-8]
t1			=	time_yr+"-"+time_month+"-"+time_day+"T"+time_hr+":"+time_min+":"+time_s
t2			=	Time(t1,format="isot")
jd			=	t2.jd

d_2000		=	jd-JD2000

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
	file.read(1028)
	i		+=	1056

file.close()

l			=	len(gps)

d_inc		=	[d_2000 for i in range(l)]

UT			=	float(time_hr)+float(time_min)/60+float(time_s)/3600
UT_inc		=	[UT for i in range(l)]

g_delay		=	np.zeros(l)
baseline	=	Baseline[das_id_1-1][das_id_2-1]
LST			=	[0 for i in range(l)]
ha			=	[0 for i in range(l)]
Alt			=	[0 for i in range(l)]
za			=	[0 for i in range(l)]

for i in range(l):
	d_inc[i]	+=	float(gps[i])*1.15741*10**(-5)
	UT_inc[i]	+=	float(gps[i])/3600
	LST[i]		=	100.46+0.985647*d_inc[i]+Long+15*UT_inc[i]
	ha[i]		=	LST[i]-RA
	Alt[i]		=	np.arcsin(np.sin(dec)*np.sin(Lat)+np.cos(dec)*np.cos(ha[i])*np.cos(Lat))
	za[i]		=	90-Alt[i]
	g_delay[i]	=	baseline*np.sin(za[i])/c

f1	=	open(file_name[:-8]+'_g_delay.txt','w+')

for i in range(l):
	f1.write(str(g_delay[i]*(10**9)))
	f1.write('\n')

f1.close()
#g_delay_t.write(file_name[:-8]+'_g_delay.txt',format='ascii.no_header',overwrite='true')

plt.plot(UT_inc,g_delay*(10**9))
plt.xlabel('Time(hrs)')
plt.ylabel('Geometric Delay(ns)')
plt.show()
