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

Baseline	=	np.zeros((8,8))
Gamma		=	np.zeros((8,8))

lamb=1.5

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

UT			=	time_hr+float(time_min)/60+float(time_s)/3600
gps			=	[]

#######

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

if(mode==2):
	g0			=	gps[0]
	for i in range(l):
		gps[i]	-=	g0

u			=	np.zeros((8*7,l))
v			=	np.zeros((8*7,l))
w			=	np.zeros((8*7,l))

lst_i		=	np.zeros(l)
ha			=	np.zeros(l)
Alt			=	np.zeros(l)
Azm			=	np.zeros(l)
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

X	=	np.zeros((8*7,l))
Y	=	np.zeros((8*7,l))
Z	=	np.zeros((8*7,l))
ctr			=	0
for das_id_1 in range(8):
	for das_id_2 in range(8):
		if(das_id_1!=das_id_2):
			b		=	Baseline[das_id_1][das_id_2]
			gamma	=	Gamma[das_id_1][das_id_2]
			H=pi-gamma*d2r
			for i in range(l):
				X[ctr][i]	=	b/lamb*(-np.sin(Lat*d2r)*np.cos(2*pi-gamma*d2r))
				Y[ctr][i]	=	b/lamb*np.sin(Azm[i]*d2r)
				Z[ctr][i]	=	b/lamb*(-np.cos(Lat*d2r)*np.cos(2*pi-gamma*d2r))
				u[ctr][i]	=	(np.sin(ha[i]*d2r)*X[ctr][i])+(np.cos(ha[i]*d2r)*Y[ctr][i])
				v[ctr][i]	=	(-1*np.sin(Alt[i]*d2r)*np.cos(ha[i]*d2r))*X[ctr][i]+(np.sin(Alt[i]*d2r)*np.sin(ha[i]*d2r)*Y[ctr][i])+(np.cos(Alt[i]*d2r)*Z[ctr][i])
				w[ctr][i]	=	(np.cos(Alt[i]*d2r)*np.cos(ha[i]*d2r)*X[ctr][i])-(np.cos(Alt[i]*d2r)*np.sin(ha[i]*d2r)*Y[ctr][i])+(np.sin(Alt[i]*d2r)*Z[ctr][i])
			uv[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	1
			ctr			+=	1
			print ctr,'/56'

dirty	=	np.fft.fft2(uv)
dirty_c	=	np.fft.fftshift(dirty)
dirty_i	=	np.fft.ifft2(uv)
dirty_i_c	=	np.fft.fftshift(dirty_i)

#fig = plt.figure()
#x = fig.gca(projection='3d')
#for i in range(56):
#	ax.plot(u[i], v[i], w[i])
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
