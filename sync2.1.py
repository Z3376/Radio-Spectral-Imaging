#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pyfftw
from astropy.table import Table

def trans(a):
    if(np.abs(a[0]-a[1])==1):
        return 1
    else:
        return 0

def amplitude(a):
	a							=	(a.imag**2+a.real**2)**0.5
	return a


if(sys.argv[1]=='-d'):
	file_name_1	=	'raw_ch03_B_P_20170714_042859_000.mbr'
	file_name_2	=	'raw_ch06_B_P_20170714_042859_000.mbr'
else:
	file_name_1	=	raw_input('File_1:')
	file_name_2	=	raw_input('File_2:')

f1	=	open(file_name_1,'rb')
f2	=	open(file_name_2,'rb')

s1	=	os.path.getsize(file_name_1)
s2	=	os.path.getsize(file_name_2)

f1.read(s1-1030)
g_l1	=	int(str(f1.read(2)).encode('hex'),16)
p_l1	=	int(str(f1.read(4)).encode('hex'),16)

f2.read(s2-1030)
g_l2	=	int(str(f2.read(2)).encode('hex'),16)
p_l2	=	int(str(f2.read(4)).encode('hex'),16)

flag	=	1
p_num1	=	np.zeros(2)
p_num2	=	np.zeros(2)

f1.seek(0)
f2.seek(0)
f1.read(26)
f2.read(26)

while(flag):
	gps1	=	int(str(f1.read(2)).encode('hex'),16)
	gps2	=	int(str(f2.read(2)).encode('hex'),16)
	if(gps1==gps2):
		p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
		p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
		o_t		=	g_l1-gps1
		n_p		=	p_l1-p_num1[0]
		flag	=	0
	elif(gps1>gps2):
		while(gps2!=gps1):
			f2.read(1054)
			gps2	=	int(str(f2.read(2)).encode('hex'),16)
		p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
		p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
		o_t		=	g_l2-gps1
		n_p		=	p_l2-p_num2[0]
		flag	=	0
	else:
		while(gps1!=gps2):
			f1.read(1054)
			gps1	=	int(str(f1.read(2)).encode('hex'),16)
		p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
		p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
		o_t		=	g_l1-gps1
		n_p		=	p_l1-p_num1[0]
		flag	=	0

a_count	=	int(n_p/o_t)
#n	=	int(1056*(1-(a_count-int(a_count))))
#print a_count, n

o_t			=	4		##############!!!!!!!!!!!!################

buff_1      =   np.zeros(1024,dtype='complex')
buff_2      =   np.zeros(1024,dtype='complex')
x_fftw_1    =   pyfftw.empty_aligned(512,dtype='complex')
y_fftw_1    =   pyfftw.empty_aligned(512,dtype='complex')
x_fftw_2    =   pyfftw.empty_aligned(512,dtype='complex')
y_fftw_2    =   pyfftw.empty_aligned(512,dtype='complex')
x_fftw_1    =   pyfftw.FFTW(x_fftw_1,x_fftw_1)
y_fftw_1    =   pyfftw.FFTW(y_fftw_1,y_fftw_1)
x_fftw_2    =   pyfftw.FFTW(x_fftw_2,x_fftw_2)
y_fftw_2    =   pyfftw.FFTW(y_fftw_2,y_fftw_2)
XX          =   np.zeros((o_t,256),dtype='complex')
YY          =   np.zeros((o_t,256),dtype='complex')
XY          =   np.zeros((o_t,256),dtype='complex')
YX          =   np.zeros((o_t,256),dtype='complex')

l_t     =   [i for i in range(256)]
for j in range(128):
    l_t[128+j]  -=  256

q           =   0

print 'Starting Corr'
print p_num1[q%2],p_num2[q%2]

for i in range(o_t):
	for j in range(1024*a_count):
		buff_1[j%1024]	=	l_t[int(str(f1.read(1)).encode('hex'),16)]
		buff_2[j%1024]	=	l_t[int(str(f2.read(1)).encode('hex'),16)]
		if(j%1024==1023):
			X_1			=	x_fftw_1(buff_1[::2])
			Y_1			=	y_fftw_1(buff_1[1::2])
			X_2			=	x_fftw_2(buff_2[::2])
			Y_2			=	y_fftw_2(buff_2[1::2])
			XX[i]		+=	X_1[:256]*np.conj(X_2[:256])
			YY[i]		+=	Y_1[:256]*np.conj(Y_2[:256])
			XY[i]		+=	X_1[:256]*np.conj(Y_2[:256])
			YX[i]		+=	Y_1[:256]*np.conj(X_2[:256])
			f1.read(28)
			f2.read(28)
			q			+=	1
			p_num1[q%2]	=	int(str(f1.read(4)).encode('hex'),16)
			p_num2[q%2]	=	int(str(f2.read(4)).encode('hex'),16)
			if(not(trans(p_num1)) or not(trans(p_num2))):
				print 'slip at:',p_num1[q%2],p_num2[q%2]
				if(not(trans(p_num1))):
					e	=	np.abs(p_num1[0]-p_num1[1])
					j	+=	e*1024
					print e
					f2.read(int(1056*e-4))
					p_num2[q%2]	=	int(str(f2.read(4)).encode('hex'),16)
				elif(not(trans(p_num2))):
					e	=	np.abs(p_num2[0]-p_num2[1])
					j	+=	e*1024
					print e
					f1.read(int(1056*e-4))
					p_num1[q%2]	=	int(str(f1.read(4)).encode('hex'),16)
				else:
					e1	=	np.abs(p_num1[0]-p_num1[1])
					e2	=	np.abs(p_num2[0]-p_num2[1])
					if(e1>e2):
						j	+=	e1*1024
						print e1
						f2.read(int(1056*(e1-e2)-4))
						p_num2[q%2]	=	int(str(f2.read(4)).encode('hex'),16)
					elif(e1<e2):
						j	+=	e2*1024
						print e2
						f1.read(int(1056*(e2-e1)-4))
						p_num1[q%2]	=	int(str(f1.read(4)).encode('hex'),16)
					else:
						j	+=	e1*1024
						print e1
	XX[i][0]	=	XX[i][1]
	XY[i][0]	=	XY[i][1]
	YY[i][0]	=	YY[i][1]
	YX[i][0]	=	YX[i][1]
	print i

print 'Corr Done'
print 'Saving Results'


XX_t	=	Table(XX)
XX_t.write(file_name_1[:-4]+'=X='+file_name_2[:-4]+'.txt',format='ascii.no_header',overwrite='true')
fw		=	open(file_name_1[:-4]+'=X='+file_name_2[:-4]+'.txt',mode='a')
fw.seek(0,os.SEEK_END)
XY_t	=	Table(XY)
XY_t.write(fw,format='ascii.no_header')
YY_t	=	Table(YY)
YY_t.write(fw,format='ascii.no_header')
YX_t	=	Table(YX)
YX_t.write(fw,format='ascii.no_header')
fw.close()

fig	=	plt.figure(1)
plt.subplot(221)
plt.plot(XX)
plt.ylabel('XX')
plt.subplot(222)
plt.plot(XY)
plt.ylabel('XY')
plt.subplot(223)
plt.plot(YX)
plt.ylabel('YX')
plt.subplot(224)
plt.plot(YY)
plt.ylabel('YY')
fig.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'.png')

fig2=	plt.figure(2)
plt.subplot(221)
plt.imshow(amplitude(XX),,interpolation='nearest',aspect='auto')
plt.ylabel('XX')
plt.subplot(222)
plt.imshow(amplitude(XY),,interpolation='nearest',aspect='auto')
plt.ylabel('XY')
plt.subplot(223)
plt.imshow(amplitude(YX),,interpolation='nearest',aspect='auto')
plt.ylabel('YX')
plt.subplot(224)
plt.imshow(amplitude(YY),,interpolation='nearest',aspect='auto')
plt.ylabel('YY')
fig2.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_im.png')
plt.show()
exit()
