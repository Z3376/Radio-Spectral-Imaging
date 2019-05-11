#!/usr/bin/python
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from astropy.time import Time
from os.path import getsize
from os.path import isfile
from subprocess import call
from astropy.table import Table

c			=	299792458
pi			=	3.14159265359
d2r			=	pi/180
r2d			=	180/pi

#Lat			=	89.9999
Lat			=	13.6029845
Long		=	77.4279978

def magnitude(a):
	mag	=	np.sqrt((a[0]**2)+(a[1]**2)+(a[2]**2))
	return mag
def lst_f(h,m,s):
	lst	=	h+float(m)/60+float(s)/3600
	return lst

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

if(len(sys.argv)==5):
	mode	=	int(sys.argv[4])		#Mode = 1.File with GPS from start of the file ; 2.File with GPS from midnight ; 3.Date, time and source given
	if(mode!=3):
		file_name	=	sys.argv[1]
		das			=	int(sys.argv[2])
	else:
		print('Enter correct mode (Mode = 1.File with GPS from start of the file ; 2.File with GPS from midnight ; 3.Date, time and source given)')
		exit()
	if(sys.argv[3]=='cas'):
		RA			=	23.39056
		dec			=	58.8
	elif(sys.argv[3]=='cyg'):
		RA			=	19.9912101667
		dec			=	40.7339155556
	elif(sys.argv[3]=='crab'):
		RA			=	5.57556
		dec			=	22.0145
	elif(sys.argv[3]=='sun'):
		RA			=	1.6281778
		dec			=	10.16191667
	elif(sys.argv[3]=='vir'):
		RA			=	12.5137287167
		dec			=	12.3911233056
	elif(sys.argv[3]=='orion'):
		RA			=	5.5880555556
		dec			=	-5.3911111111
	elif(sys.argv[3]=='pole'):
		RA			=	2.5
		dec			=	89.000
elif(len(sys.argv)==6):
	mode	=	int(sys.argv[5])		#Mode = 1.File with GPS from start of the file ; 2.File with GPS from midnight ; 3.Date, time and source given
	if(mode==3):
		l	=	int(float(sys.argv[3])*3600)
		ut	=	sys.argv[2]
		date	=	sys.argv[1]
	else:
		print('Enter correct mode (Mode = 1.File with GPS from start of the file ; 2.File with GPS from midnight ; 3.Date, time and source given)')
		exit()
	if(sys.argv[4]=='cas'):
		RA			=	23.39056
		dec			=	58.8
	elif(sys.argv[4]=='cyg'):
		RA			=	19.9912101667
		dec			=	40.7339155556
	elif(sys.argv[4]=='crab'):
		RA			=	5.57556
		dec			=	22.0145
	elif(sys.argv[4]=='sun'):
		RA			=	1.6281778
		dec			=	10.16191667
	elif(sys.argv[4]=='vir'):
		RA			=	12.5137287167
		dec			=	12.3911233056
	elif(sys.argv[4]=='orion'):
		RA			=	5.5880555556
		dec			=	-5.3911111111
	elif(sys.argv[4]=='pole'):
		RA			=	2.5
		dec			=	89.000
else:
	mode		=	input('Mode(GPS from start of the file(1) or midnight(2)) or GPS seconds and start_UT given(3): ')
	if(mode==3):
		date	=	input('Enter Date: ')
		ut		=	input('Start UTC: ')
		l		=	float(input('Enter no. of hours: '))*3600
	else:
		file_name	=	input('File_name: ')
		das			=	int(input('Das_ID: '))
	RA			=	input('RA(in hrs): ')
	dec			=	input('dec(in deg): ')
	

x	=	np.zeros(8,dtype='float')
y	=	np.zeros(8,dtype='float')

x[0]	=	0.0
y[0]	=	0.0
x[1]	=	-64.5
y[1]	=	41.6
x[2]	=	-39.5
y[2]	=	68.0
x[3]	=	27.3
y[3]	=	99.7
x[4]	=	-19.9
y[4]	=	29.5
x[5]	=	-11.4
y[5]	=	58.1
x[6]	=	-62.6
y[6]	=	25.9
x[7]	=	26.0
y[7]	=	53.8

Baseline	=	np.zeros((8,8))
Gamma		=	np.zeros((8,8))

if(not(isfile('Baseline.txt') or not(isfile('Gamma.txt')))):
	call(['./coord2.py'])

fbasel	=	open('Baseline.txt','r')
fgamma	=	open('Gamma.txt','r')

read_file(fbasel,Baseline,8,8)
read_file(fgamma,Gamma,8,8)

if(mode==3):
	gps	=	[i for i in range(l)]

	ut_x	= [float(xi) for xi in ut.split(':')]
	time_hr		=	ut_x[0]
	time_min	=	ut_x[1]
	time_s		=	ut_x[2]
	UT	=	time_hr+time_min/60+time_s/3600

	date_x	=	[int(xi) for xi in date.split('/')]
	time_yr		=	date_x[2]
	time_month	=	date_x[1]
	time_day	=	date_x[0]
	UT	=	time_hr+time_min/60+time_s/3600

	ofname	=	sys.argv[4]+'_'+date.replace('/','')+'_'+ut+'_'+sys.argv[3]

else:
	ofname	=	file_name[:-8]

	time_yr		=	int(file_name[-23:-19])
	time_month	=	int(file_name[-19:-17])
	time_day	=	int(file_name[-17:-15])
	time_hr		=	int(file_name[-14:-12])-5
	time_min	=	int(file_name[-12:-10])-30
	time_s		=	int(file_name[-10:-8])

	while(time_min<0):
		time_min	+=	60
		time_hr		-=	1
	while(time_hr<0):
		time_hr		+=	24
		time_day	-=	1
		while(time_day<1):
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
			while(time_month<1):
				time_month	+=	12
				time_yr		-=	1

	UT			=	time_hr+float(time_min)/60+float(time_s)/3600
	gps			=	[]
	jn			=	0
	while(jn<999):
		j_d		=	jn/10
		d		=	1
		while(j_d):
			j_d	=	j_d/10
			d	+=	1
		s		=	('_'+'0'*(3-d))+str(jn)
		
		if (not(isfile(file_name[:-8]+s+".hdr"))):
			break
		print("Reading file "+file_name[:-8]+s+".hdr")

		jn			+=	1
		file		=	open(file_name[:-8]+s+".hdr","rb")
		size 		=	getsize(file_name[:-8]+s+".hdr")
		i			=	0

		g			=	np.empty(2)
		q			=	0
		while(i<size):
			file.read(7)
			d		=	file.read(1)
			file.read(18)
			g_buff	=	int(str(file.read(2)).encode('hex'),16)
			if(d==str(das)):
				q	+=	1
				g[q%2]	=	g_buff
				if(g[q%2]-g[(q+1)%2]!=0):
					gps.append(g[q%2])
			file.read(4)
			i		+=	32
		file.close()

	l			=	len(gps)
	#print l

	if(mode==2):
		g0			=	gps[0]
		for i in range(l):
			gps[i]	-=	g0

u			=	np.zeros((8*7,l))
v			=	np.zeros((8*7,l))
#w			=	np.zeros((8*7,l))
N			=	np.zeros(3)
N[1]		=	1
E			=	np.zeros(3)
E[0]		=	1
src			=	np.zeros(3)

UP			=	np.zeros(3)
UP[2]		=	1

a			=	np.zeros(3)
u_a			=	np.zeros(3)
v_a			=	np.zeros(3)

lst_i		=	[0 for i in range(l)]
ha			=	[0 for i in range(l)]
Alt			=	[0 for i in range(l)]
Azm			=	[0 for i in range(l)]
UT_i		=	[UT for i in range(l)]

for i in range(l):
	UT_i[i]	+=	float(gps[i])/3600
	if(i==0):
		time_s	+=	gps[i]
	else:
		time_s	+=	gps[i]-gps[i-1]
	if(time_s<0 or time_min<0 or time_hr<0 or time_day<0 or time_month<0):
		break
	while(time_s>=60):
		time_s		-=	60
		time_min	+=	1
		while(time_min>=60):
			time_min	-=	60
			time_hr		+=	1
			while(time_hr>=24):
				time_hr		-=	24
				time_day	+=	1
				while(time_month==2 and time_yr%4==0 and time_day>29):
					time_month	+=	1
					time_day	-=	29
				while(time_month==2 and time_yr%4!=0 and time_day>28):
					time_month	+=	1
					time_day	-=	28
				while((time_month==1 or time_month==3 or time_month==5 or time_month==7 or time_month==8 or time_month==10 or time_month==12) and time_day>31):
					time_month	+=	1
					time_day	-=	31
				while((time_month==4 or time_month==6 or time_month==9 or time_month==11) and time_day>30):
					time_month	+=	1
					time_day	-=	30
				while(time_month>12):
					time_month	-=	12
					time_yr		+=	1
	t1			=	str(int(time_yr))+"-"+str(int(time_month))+"-"+str(int(time_day))+"T"+str(int(time_hr))+":"+str(int(time_min))+":"+str(int(time_s))
	#print t1
	t2			=	Time(t1,format="isot",location=(str(Long),str(Lat)))
	lst			=	str(t2.sidereal_time('mean'))
	j			=	0
	ctr_h		=	0
	h_flag		=	0
	ctr_m		=	0
	m_flag		=	0
	ctr_s		=	-1
	s_flag		=	0
	while(j<len(lst)):
		if(lst[j]!='h' and h_flag==0):
			ctr_h	+=	1
		elif(lst[j]=='h'):
			h_flag	=	1
		elif(lst[j]!='m' and m_flag==0):
			ctr_m	+=	1
		elif(lst[j]=='m'):
			m_flag	=	1
		else:
			ctr_s	+=	1
		j	+=	1
	lst_hr		=	int(lst[:ctr_h])
	lst_min		=	int(lst[ctr_h+1:ctr_h+ctr_m+1])
	lst_s		=	float(lst[ctr_h+ctr_m+2:ctr_h+ctr_m+ctr_s+2])
	lst_i[i]	=	lst_f(lst_hr,lst_min,lst_s)
	ha[i]		=	(lst_i[i]-RA)*15
	Alt[i]		=	(np.arcsin(np.sin(dec*d2r)*np.sin(Lat*d2r)+np.cos(dec*d2r)*np.cos(ha[i]*d2r)*np.cos(Lat*d2r)))*r2d
	#Azm[i]		=	(np.arcsin(np.cos(dec*d2r)*np.sin(ha[i]*d2r)/np.cos(Alt[i]*d2r)))*r2d
	Azm[i]		=	(np.arccos((np.sin(dec*d2r)-np.sin(Alt[i]*d2r)*np.sin(Lat*d2r))/np.cos(Alt[i]*d2r)/np.cos(Lat*d2r)))*r2d
	#Azm[i]		=	(np.arctan2((np.cos(dec*d2r)*np.sin(ha[i]*d2r)),(np.sin(dec*d2r)*np.cos(Lat*d2r)-np.cos(dec*d2r)*np.cos(ha[i]*d2r)*np.sin(Lat*d2r))))*r2d
	if(np.sin(ha[i]*d2r)>0):
		Azm[i]	=	360-Azm[i]

uv	=	np.zeros((240,240))

ctr			=	0
for das_id_1 in range(8):
	for das_id_2 in range(8):
		if(das_id_1!=das_id_2):
			a[0]		=	x[das_id_1]-x[das_id_2]
			a[1]		=	y[das_id_1]-y[das_id_2]
			#print x[das_id_1],':',y[das_id_1],'::',x[das_id_2],':',y[das_id_2]
			for i in range(l):
				src[0]		=	np.cos(Alt[i]*d2r)*np.sin(Azm[i]*d2r)
				src[1]		=	np.cos(Alt[i]*d2r)*np.cos(Azm[i]*d2r)
				src[2]		=	np.sin(Alt[i]*d2r)
				v_a[0]		=	-1*np.sin(Alt[i]*d2r)*np.sin(Azm[i]*d2r)
				v_a[1]		=	-1*np.sin(Alt[i]*d2r)*np.cos(Azm[i]*d2r)
				v_a[2]		=	np.cos(Alt[i]*d2r)
				#u_a			=	np.cross(v_a,src)
				u_a[0]		=	-1*np.cos(Azm[i]*d2r)
				u_a[1]		=	np.sin(Azm[i]*d2r)
				T			=	np.cross(a,src)
				T			=	np.cross(src,T)
				T			=	T/magnitude(T)
				T			=	magnitude(a)*np.sin(Alt[i]*d2r)*T
				u[ctr][i]	=	np.dot(u_a,T)
				v[ctr][i]	=	np.dot(v_a,T)
				uv[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	1
			ctr			+=	1
			print ctr,'/56'

dirty	=	np.fft.fft2(uv)
dirty_c	=	np.fft.fftshift(dirty)
dirty_i	=	np.fft.ifft2(uv)
dirty_i_c	=	np.fft.fftshift(dirty_i)

print('Saving Results')

fw	=	open(os.getcwd()+'/'+ofname+'_uv5.txt',"w+")
fw.write(str(l))
fw.write('\n')
for i in range(l):
	for j in range(ctr):
		fw.write(str(u[j][i]))
		fw.write(' ')
	fw.write('\n')
for i in range(l):
	for j in range(ctr):
		fw.write(str(v[j][i]))
		fw.write(' ')
	fw.write('\n')

fig1=plt.figure()
for i in range(56):
	plt.plot(u[i],v[i],color='k')
plt.xlabel('U')
plt.ylabel('V')
plt.axis([-125, 125, -125, 125])
fig1.savefig(ofname+'_uv5.png')
fig2=plt.figure()
plt.imshow(uv.T,origin='lower')
fig2.savefig(ofname+'_uv5_digi.png')
fig3=plt.figure()
plt.imshow(np.real(dirty),cmap='afmhot',origin='lower')
fig3.savefig(ofname+'_uv5_digi_realfft.png')
fig4=plt.figure()
plt.imshow(np.imag(dirty),cmap='afmhot',origin='lower')
fig4.savefig(ofname+'_uv5_digi_imagfft.png')
fig5=plt.figure()
plt.imshow(np.absolute(dirty),cmap='afmhot',origin='lower')
fig5.savefig(ofname+'_uv5_digi_fft.png')
fig6=plt.figure()
plt.imshow(np.real(dirty_i),cmap='afmhot',origin='lower')
fig6.savefig(ofname+'_uv5_digi_realifft.png')
fig7=plt.figure()
plt.imshow(np.imag(dirty_i),cmap='afmhot',origin='lower')
fig7.savefig(ofname+'_uv5_digi_imagifft.png')
fig8=plt.figure()
plt.imshow(np.absolute(dirty_i),cmap='afmhot',origin='lower')
fig8.savefig(ofname+'_uv5_digi_ifft.png')
fig9=plt.figure()
plt.imshow(np.real(dirty_c),cmap='afmhot',origin='lower')
fig9.savefig(ofname+'_uv5_digi_c_realfft.png')
fig10=plt.figure()
plt.imshow(np.imag(dirty_c),cmap='afmhot',origin='lower')
fig10.savefig(ofname+'_uv5_digi_c_imagfft.png')
fig11=plt.figure()
plt.imshow(np.absolute(dirty_c),cmap='afmhot',origin='lower',clim=[0,800])
plt.colorbar()
fig11.savefig(ofname+'_uv5_digi_c_fft.pdf')
fig12=plt.figure()
plt.imshow(np.real(dirty_i_c),cmap='afmhot',origin='lower')
fig12.savefig(ofname+'_uv5_digi_c_realifft.png')
fig13=plt.figure()
plt.imshow(np.imag(dirty_i_c),cmap='afmhot',origin='lower')
fig13.savefig(ofname+'_uv5_digi_c_imagifft.png')
fig14=plt.figure()
plt.imshow(np.absolute(dirty_i_c),cmap='afmhot',origin='lower',clim=[0,0.005])
plt.colorbar()
fig14.savefig(ofname+'_uv5_digi_c_ifft.pdf')
plt.show()
