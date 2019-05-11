#!/usr/bin/python

# Cross-correlation averaging time limit calculation due to change in apparent baselines (UV)
# by Harsh_Grover

import sys
from os.path import getsize
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

c			=	299792458
pi			=	3.14159265359
d2r			=	pi/180

Lat			=	13.6029845
Long		=	77.4279978

def read_file(f,a,n1,n2):
	ctr1	=	0
	ctr2	=	0
	ctr3	=	0
	f.seek(0,0)
	while(ctr2<n1):
		while(ctr3<n2):
			b	=	f.read(1)
			if(b=='\n'):
				f.seek(-ctr1-1,1)
				a[ctr2][ctr3]	=	float(f.read(ctr1))
				f.read(1)
				ctr1	=	0
				ctr3	+=	1
			else:
				ctr1	+=	1
		ctr2	+=	1
		ctr3	=	0

wvlngth		=	input('Wavelength:')

Baseline		=	np.zeros((8,8))

fbasel	=	open('Baseline.txt','r')
read_file(fbasel,Baseline,8,8)

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

t1			=	str(int(time_yr))+"-"+str(int(time_month))+"-"+str(int(time_day))+"T"+str(int(time_hr))+":"+str(int(time_min))+":"+str(int(time_s))
#print t1
t2			=	Time(t1,format="isot",location=(str(Long),str(Lat)))

LST			=	t2.sidereal_time('mean')

lst			=	str(LST)

das_id1		=	int(sys.argv[2])
das_id2		=	int(sys.argv[3])

baseline	=	Baseline[das_id1-1][das_id2-1]

i			=	0
ctr_h		=	0
h_flag		=	0
ctr_m		=	0
m_flag		=	0
ctr_s		=	-1
s_flag		=	0
while(i<len(lst)):
	if(lst[i]!='h' and h_flag==0):
		ctr_h	+=	1
	elif(lst[i]=='h'):
		h_flag	=	1
	elif(lst[i]!='m' and m_flag==0):
		ctr_m	+=	1
	elif(lst[i]=='m'):
		m_flag	=	1
	else:
		ctr_s	+=	1
	i	+=	1
lst_hr	=	int(lst[:ctr_h])
lst_min	=	int(lst[ctr_h+1:ctr_h+ctr_m+1])
lst_s	=	float(lst[ctr_h+ctr_m+2:ctr_h+ctr_m+ctr_s+2])

lst_i	=	(lst_hr+float(lst_min)/60+lst_s/3600)*15

ha			=	lst_i-RA*15
Alt			=	np.arcsin(np.sin(dec*d2r)*np.sin(Lat*d2r)+np.cos(dec*d2r)*np.cos(ha*d2r)*np.cos(Lat*d2r))
t			=	1200*wvlngth/pi/baseline/np.cos(Alt)

Alt_d		=	Alt*(180/pi)

print 'Max averaging time allowed:',t,'s'
#print LST
#print Alt_d
