#!/usr/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.time import Time
from os.path import getsize
from os.path import isfile
from subprocess import call

c			=	299792458
pi			=	np.pi
d2r			=	pi/180
r2d			=	180/pi

Lat			=	13.6029845
Long		=	77.4279978

def magnitude(a):
	mag	=	np.sqrt((a[0]**2)+(a[1]**2)+(a[2]**2))
	return mag
def lst_f(h,m,s):
	lst	=	h+float(m)/60+float(s)/3600
	return lst

def longlat2ECEF(lat, lon, alt, rad):
	f	=	float(1.0/298.257223563)
	ls	=	np.arctan(((1-f)**2)*np.tan(lat))
	x	=	rad*np.cos(ls)*np.cos(lon)+alt*np.cos(lat)*np.cos(lon)
	y	=	rad*np.cos(ls)*np.sin(lon)+alt*np.cos(lat)*np.sin(lon)
	z	=	rad*np.sin(ls)+alt*np.sin(lat)
	return (x, y, z)

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

file_name	=	sys.argv[1]
das			=	int(sys.argv[2])


if(len(sys.argv)==5):
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
		dec			=	89.00
	mode		=	sys.argv[4]
elif(len(sys.argv)==4):
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
		dec			=	89.00
	mode		=	input('Mode(GPS from start of the file(1) or midnight(2)): ')
else:
	RA			=	input('RA(in hrs): ')
	dec			=	input('dec(in deg): ')
	mode		=	input('Mode(GPS from start of the file(1) or midnight(2)): ')

x	=	np.zeros(8,dtype='float')
y	=	np.zeros(8,dtype='float')
z	=	np.zeros(8,dtype='float')

x[0]		=	0.0
y[0]		=	0.0
x[1]		=	-64.5
y[1]		=	41.6
x[2]		=	-39.5
y[2]		=	68.0
x[3]		=	27.3
y[3]		=	99.7
x[4]		=	-19.9
y[4]		=	29.5
x[5]		=	-11.4
y[5]		=	58.1
x[6]		=	-62.6
y[6]		=	25.9
x[7]		=	26.0
y[7]		=	53.8

Baseline	=	np.zeros((8,8))
Gamma		=	np.zeros((8,8))

if(not(isfile('Baseline.txt') or not(isfile('Gamma.txt')))):
	call(['./coord2.py'])

fbasel	=	open('Baseline.txt','r')
fgamma	=	open('Gamma.txt','r')

read_file(fbasel,Baseline,8,8)
read_file(fgamma,Gamma,8,8)

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

#print time_hr,time_min, time_s

UT			=	time_hr+float(time_min)/60+float(time_s)/3600
#print UT
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
	print "Reading file "+file_name[:-8]+s+".hdr"

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

#Lat_t		=	np.zeros(8)
#Long_t		=	np.zeros((8,l))
#Lat_t[0]	=	13.602427
#Long_t[0][:]=	77.427781
#Lat_t[1]	=	13.602803
#Long_t[1][:]=	77.427185
#Lat_t[2]	=	13.603042
#Long_t[2][:]=	77.427416
#Lat_t[3]	=	13.603328
#Long_t[3][:]=	77.428033
#Lat_t[4]	=	13.602694
#Long_t[4][:]=	77.427597
#Lat_t[5]	=	13.602952
#Long_t[5][:]=	77.427676
#Lat_t[6]	=	13.602661
#Long_t[6][:]=	77.427203
#Lat_t[7]	=	13.602913
#Long_t[7][:]=	77.428021

r			=	6378137.0

#for i in range(8):
#	(x[i],y[i],z[i])=longlat2ECEF(Lat_t[i],Long_t[i][0], 0, r)

#plt.scatter(x,y)
#plt.show()
#exit()

#for i in range(8):
#	for j in range(l):
#		Long_t[i][j]	+=	gps[j]*float(360)/86164

#x	=	np.zeros((8,l),dtype='float')
#y	=	np.zeros((8,l),dtype='float')
#z	=	np.zeros((8,l),dtype='float')

#for i in range(8):
#	for j in range(l):
#		(x[i][j],y[i][j],z[i][j])=longlat2ECEF(Lat_t[i]*d2r,Long_t[i][j]*d2r,100,r)

#x_s	=	np.zeros(l)
#y_s	=	np.zeros(l)
#z_s	=	np.zeros(l)

#for i in range(l):
#	x_s[i]	=	x[0][i]
#	y_s[i]	=	y[0][i]
#	z_s[i]	=	z[0][i]

#for i in range(8):
#	for j in range(l):
#		x[i][j]	-=	x_s[j]
#		y[i][j]	-=	y_s[j]
#		z[i][j]	-=	z_s[j]

#for i in range(8):
#	plt.scatter(x[i][1],y[i][1])
#	plt.scatter(x[i][50],y[i][50])
#plt.show()
#exit()

u			=	np.zeros((8*7,l))
v			=	np.zeros((8*7,l))
w			=	np.zeros((8*7,l))
#N			=	[0,1,0]
#E			=	[1,0,0]

N			=	np.zeros(3)
N[1]		=	1
E			=	np.zeros(3)
E[0]		=	1
UP			=	np.zeros(3)
UP[2]		=	1

a			=	np.zeros(3)
u_a			=	np.zeros((l,3))
v_a			=	np.zeros((l,3))
w_a			=	np.zeros(3)
src			=	np.zeros(3)

lst_i		=	np.zeros(l)
ha			=	np.zeros(l)
Alt			=	np.zeros(l)
Azm			=	np.zeros(l)
UT_i		=	[UT for i in range(l)]

#u_a[0]		=	np.cos(RA*d2r)
#u_a[1]		=	-1*np.sin(RA*d2r)
#v_a[0]		=	-1*np.sin(RA*d2r)
#v_a[1]		=	-1*np.cos(RA*d2r)
#v_a[2]		=	np.cos(dec*d2r)
#w_a[0]		=	np.cos(dec*d2r)*np.sin(RA*d2r)
#w_a[1]		=	np.cos(dec*d2r)*np.cos(RA*d2r)
#w_a[2]		=	np.sin(dec*d2r)

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

#l_t	=	[i for i in range(240)]

ctr			=	0
for das_id_1 in range(8):
	for das_id_2 in range(8):
		if(das_id_1!=das_id_2):
			#print x[das_id_1],':',y[das_id_1],'::',x[das_id_2],':',y[das_id_2]
			a[0]		=	x[das_id_1]-x[das_id_2]
			a[1]		=	y[das_id_1]-y[das_id_2]
			a[2]		=	z[das_id_1]-z[das_id_2]
			for i in range(l):
				#u_a[0]		=	np.cos(Azm[i]*d2r)
				#u_a[1]		=	-1*np.sin(Azm[i]*d2r)
				#v_a[0]		=	-1*np.sin(Azm[i]*d2r)
				#v_a[1]		=	-1*np.cos(Azm[i]*d2r)
				#v_a[2]		=	np.cos(Alt[i]*d2r)
				#w_a[0]		=	np.cos(Alt[i]*d2r)*np.sin(Azm[i]*d2r)
				#w_a[1]		=	np.cos(Alt[i]*d2r)*np.cos(Azm[i]*d2r)
				#w_a[2]		=	np.sin(Alt[i]*d2r)
				src[0]		=	np.cos(Alt[i]*d2r)*np.sin(Azm[i]*d2r)
				src[1]		=	np.cos(Alt[i]*d2r)*np.cos(Azm[i]*d2r)
				src[2]		=	np.sin(Alt[i]*d2r)
				src			=	src/magnitude(src)
				buff		=	np.cross(src,N)
				buff		=	buff/magnitude(buff)
				w_a			=	src
				v_a[i]		=	np.cross(src,buff)
				v_a[i]		=	v_a[i]/magnitude(v_a[i])
				u_a[i]		=	np.cross(src,v_a[i])
				u_a[i]		=	u_a[i]/magnitude(u_a[i])
				#print np.dot(v_a[i],u_a[i]),np.dot(src,u_a[i]),np.dot(src,v_a[i])
				#exit()
				u[ctr][i]	=	np.dot(a,u_a[i])
				v[ctr][i]	=	np.dot(a,v_a[i])
				w[ctr][i]	=	np.dot(a,w_a)
				#uv[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	1
			ctr			+=	1
			print ctr,'/56'

dirty	=	np.fft.fft2(uv)
dirty_c	=	np.fft.fftshift(dirty)
dirty_i	=	np.fft.ifft2(uv)
dirty_i_c	=	np.fft.fftshift(dirty_i)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#for i in range(l):
	#ax.plot(u_a[i][0],u_a[i][1],u_a[i][2])
	#ax.plot(v_a[i][0],v_a[i][1],v_a[i][2])
#plt.show()
#exit()
uv=uv.T

fig1=plt.figure()
for i in range(56):
	plt.plot(u[i],v[i],color='k')
plt.xlabel('U')
plt.ylabel('V')
fig1.savefig(file_name[:-8]+'_uv4.png')
fig2=plt.figure()
for i in range(56):
	plt.plot(v[i],w[i],color='k')
plt.xlabel('V')
plt.ylabel('W')
fig2.savefig(file_name[:-8]+'_vw.png')
fig3=plt.figure()
for i in range(56):
	plt.plot(w[i],u[i],color='k')
plt.xlabel('W')
plt.ylabel('U')
fig3.savefig(file_name[:-8]+'_wu.png')
#############
plt.show()
exit()
#############
#fig3=plt.figure()
#for i in range(56):
#	plt.plot(-1*u[i],-1*v[i],color='k')
#plt.xlabel('-U')
#plt.ylabel('-V')
#fig3.savefig(file_name[:-8]+'_uv-4.png')
fig4=plt.figure()
plt.imshow(uv,origin='lower')
fig4.savefig(file_name[:-8]+'_uv4_digi.png')
#fig3=plt.figure()
#plt.imshow(np.real(dirty),cmap='afmhot',origin='lower')
#fig3.savefig(file_name[:-8]+'_uv4_digi_realfft.png')
#fig4=plt.figure()
#plt.imshow(np.imag(dirty),cmap='afmhot',origin='lower')
#fig4.savefig(file_name[:-8]+'_uv4_digi_imagfft.png')
#fig5=plt.figure()
#plt.imshow(np.absolute(dirty),cmap='afmhot',origin='lower')
#fig5.savefig(file_name[:-8]+'_uv4_digi_fft.png')
#fig6=plt.figure()
#plt.imshow(np.real(dirty_i),cmap='afmhot',origin='lower')
#fig6.savefig(file_name[:-8]+'_uv4_digi_realifft.png')
#fig7=plt.figure()
#plt.imshow(np.imag(dirty_i),cmap='afmhot',origin='lower')
#fig7.savefig(file_name[:-8]+'_uv4_digi_imagifft.png')
#fig8=plt.figure()
#plt.imshow(np.absolute(dirty_i),cmap='afmhot',origin='lower')
#fig8.savefig(file_name[:-8]+'_uv4_digi_ifft.png')
fig5=plt.figure()
plt.subplot(211)
plt.imshow(np.real(dirty_c),cmap='afmhot',origin='lower')
plt.xlabel('digi_c_realfft')
plt.subplot(221)
plt.imshow(np.imag(dirty_c),cmap='afmhot',origin='lower')
plt.xlabel('digi_c_imagfft')
fig5.savefig(file_name[:-8]+'_uv4_digi_c_fft_comp.png')
fig51=plt.figure()
plt.imshow(np.absolute(dirty_c),cmap='afmhot',origin='lower',clim=[0,800])
plt.colorbar()
plt.xlabel('digi_c_fft')
fig51.savefig(file_name[:-8]+'_uv4_digi_c_fft.pdf')
fig6=plt.figure()
plt.subplot(211)
plt.imshow(np.real(dirty_i_c),cmap='afmhot',origin='lower')
plt.xlabel('digi_c_realifft')
plt.subplot(221)
plt.imshow(np.imag(dirty_i_c),cmap='afmhot',origin='lower')
plt.xlabel('digi_c_imagifft')
fig6.savefig(file_name[:-8]+'_uv4_digi_c_ifft_comp.png')
fig61=plt.figure()
plt.imshow(np.absolute(dirty_i_c),cmap='afmhot',origin='lower',clim=[0,0.005])
plt.colorbar()
plt.xlabel('digi_c_ifft')
fig61.savefig(file_name[:-8]+'_uv4_digi_c_ifft.pdf')
plt.show()