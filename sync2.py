#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import getsize
import pyfftw

def trans(a):
    if(np.abs(a[0]-a[1])==1):
        return 1
    else:
        return 0

if(sys.argv[1]=='-d'):
	file_name_1	=	'raw_ch03_B_P_20170714_042859_000.mbr'
	file_name_2	=	'raw_ch06_B_P_20170714_042859_000.mbr'
else:
	file_name_1	=	raw_input('File_1:')
	file_name_2	=	raw_input('File_2:')

f1	=	open(file_name_1,'rb')
f2	=	open(file_name_2,'rb')

s1	=	getsize(file_name_1)
s2	=	getsize(file_name_2)

f1.read(26)
f2.read(26)
g_s1	=	int(str(f1.read(2)).encode('hex'),16)
g_s2	=	int(str(f2.read(2)).encode('hex'),16)
p_s1	=	int(str(f1.read(4)).encode('hex'),16)
p_s2	=	int(str(f2.read(4)).encode('hex'),16)

f1.read(s1-1062)
g_l1	=	int(str(f1.read(2)).encode('hex'),16)
p_l1	=	int(str(f1.read(4)).encode('hex'),16)

f2.read(s2-1062)
g_l2	=	int(str(f2.read(2)).encode('hex'),16)
p_l2	=	int(str(f2.read(4)).encode('hex'),16)

t_l1	=	g_l1-g_s1
t_l2	=	g_l2-g_s2
n_p1	=	p_l1-p_s1
n_p2	=	p_l2-p_s2

if(t_l1<=t_l2):
	t_l	=	t_l1
else:
	t_l	=	t_l2
if(n_p1<=n_p2):
	n_p	=	n_p1
else:
	n_p	=	n_p2

a_count	=	n_p/t_l

flag	=	1
n_flag	=	0
gps1	=	0
gps2	=	0
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
		flag	=	0
	elif(gps1>gps2):
		while(gps2!=gps1):
			f2.read(1054)
			gps2	=	int(str(f2.read(2)).encode('hex'),16)
		p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
		p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
		flag	=	0
	else:
		while(gps1!=gps2):
			f1.read(1054)
			gps1	=	int(str(f1.read(2)).encode('hex'),16)
		p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
		p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
		flag	=	0

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
XX          =   np.zeros((t_l,256),dtype='complex')
YY          =   np.zeros((t_l,256),dtype='complex')
XY          =   np.zeros((t_l,256),dtype='complex')
YX          =   np.zeros((t_l,256),dtype='complex')

l_t     =   [i for i in range(256)]
for j in range(128):
    l_t[128+j]  -=  256

q           =   0

print 'Starting Corr'
print p_num1[q%2],p_num2[q%2]

#for i in range(t_l):
for i in range(4):
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

fig=plt.figure(1)
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
plt.show()
exit()
