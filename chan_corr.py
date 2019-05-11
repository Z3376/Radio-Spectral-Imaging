#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

def amplitude(a):
	a							=	(a.imag**2+a.real**2)**0.5
	return a

def mean(a):
	mean_a						=	np.zeros(256)
	for i in range(256):
		for j in range(len(a)):
			mean_a[i]			+=	amplitude(a[j][i])
		mean_a[i]				=	mean_a[i]/len(a)
	return mean_a

def reject_matrix(XX_m):
	R		=	np.zeros(256)
	for i in range(256):
		if(np.abs(amplitude(XX_m[i])-np.mean(XX_m))>np.std(XX_m)*n_s):
			R[i]	=	1
	return R

def reject(R,corr):
	for j in range(256):
		if(R[j]==1):
			for i in range(i_n):
				corr[i][j]=0

def corr_read(s1,fopen,XX_1):
	j	=	0
	i	=	0
	while(i<int(s1)):
		p		=	fopen.read(1)
		ctr1	=	-1
		ctr2	=	-1
		e_flag	=	0
		while(p!='(' or p=='0' and i<int(s1)):
			if(p=='0'):
				#print p
				XX_1.append(complex(0,0))
				i	+=	2
				fopen.read(2)
			p	=	fopen.read(1)
			i	+=	1
		if(not(i<int(s1))):
			break
		while((p!='+' or e_flag) and (p!='-' or ctr1==0 or e_flag)):
			p	=	fopen.read(1)
			ctr1	+=	1
			if(not(e_flag-ctr1)):
				e_flag	=	0
			if(p=='e'):
				e_flag	=	ctr1+2
		while(p!='j'):
			p	=	fopen.read(1)
			ctr2	+=	1
		fopen.seek(-ctr1-ctr2-2,1)
		r	=	float(fopen.read(ctr1))
		fopen.read(1)
		im	=	float(fopen.read(ctr2))
		#print complex(r,im)
		XX_1.append(complex(r,im))
		fopen.read(3)
		i	+=	ctr1+ctr2+5
		j	+=	1

file_name	=	'raw_ch04_crab_20170925_174908_000=X=raw_ch07_crab_20170925_174908_000.txt'
#file_name	=	raw_input('File: ')
n_s		=	2
fopen	=	open(file_name,'r')
size	=	os.path.getsize(file_name)
s1		=	float(size)/4
ARR		=	[]

corr_read(size,fopen,ARR)

l	= int(float(len(ARR)))

n	=	int(l/4)
i_n	=	int(n/256)

ARR_256	=	np.zeros((4,i_n,256),dtype='complex')
XX		=	np.zeros((i_n,256),dtype='complex')
XY		=	np.zeros((i_n,256),dtype='complex')
YX		=	np.zeros((i_n,256),dtype='complex')
YY		=	np.zeros((i_n,256),dtype='complex')

for k in range(4):
	for i in range(i_n):
		for j in range(256):
			ARR_256[k][i][j]	=	ARR[(k*n)+(i*256)+j]

XX	=	ARR_256[0]
XY	=	ARR_256[1]
YY	=	ARR_256[2]
YX	=	ARR_256[3]

XX_m	=	mean(XX)
XY_m	=	mean(XY)
YY_m	=	mean(YY)
YX_m	=	mean(YX)

R_XX	=	reject_matrix(XX_m)
R_XY	=	reject_matrix(XY_m)
R_YY	=	reject_matrix(YY_m)
R_YX	=	reject_matrix(YX_m)

reject(R_XX,XX)
reject(R_XY,XY)
reject(R_YY,YY)
reject(R_YX,YX)

XX_s	=	np.zeros((256,i_n),dtype='complex')
XY_s	=	np.zeros((256,i_n),dtype='complex')
YY_s	=	np.zeros((256,i_n),dtype='complex')
YX_s	=	np.zeros((256,i_n),dtype='complex')

XX_m	=	mean(XX)
XY_m	=	mean(XY)
YY_m	=	mean(YY)
YX_m	=	mean(YX)

for i in range(256):
	for j in range(i_n):
		XX_s[i][j]	=	XX[j][i]
		XY_s[i][j]	=	XY[j][i]
		YY_s[i][j]	=	YY[j][i]
		YX_s[i][j]	=	YX[j][i]

XX_m	=	XX_m-np.mean(XX_m)
XY_m	=	XY_m-np.mean(XY_m)
YY_m	=	YY_m-np.mean(YY_m)
YX_m	=	YX_m-np.mean(YX_m)

varXX_s	=	np.zeros(256,dtype='complex')
varXY_s	=	np.zeros(256,dtype='complex')
varYY_s	=	np.zeros(256,dtype='complex')
varYX_s	=	np.zeros(256,dtype='complex')

for i in range(256):
	for j in range(i_n):
		varXX_s[i]	+=	(XX_s[i][j]-np.mean(XX_s[i]))**2
		varXY_s[i]	+=	(XY_s[i][j]-np.mean(XY_s[i]))**2
		varYY_s[i]	+=	(YY_s[i][j]-np.mean(YY_s[i]))**2
		varYX_s[i]	+=	(YX_s[i][j]-np.mean(YX_s[i]))**2
	varXX_s[i]	=	np.sqrt((varXX_s[i])/i_n)
	varXY_s[i]	=	np.sqrt((varXY_s[i])/i_n)
	varYY_s[i]	=	np.sqrt((varYY_s[i])/i_n)
	varYX_s[i]	=	np.sqrt((varYX_s[i])/i_n)

chan_corrXX	=	np.zeros((i_n,256,256),dtype='complex')
chan_corrXY	=	np.zeros((i_n,256,256),dtype='complex')
chan_corrYY	=	np.zeros((i_n,256,256),dtype='complex')
chan_corrYX	=	np.zeros((i_n,256,256),dtype='complex')

chan_corrXX_m	=	np.zeros((256,256),dtype='complex')
chan_corrXY_m	=	np.zeros((256,256),dtype='complex')
chan_corrYY_m	=	np.zeros((256,256),dtype='complex')
chan_corrYX_m	=	np.zeros((256,256),dtype='complex')

for k in range(i_n):
	for i in range(256):
		for j in range(256):
			chan_corrXX[k][i][j]	=	XX_s[i][k]*np.conj(XX_s[j][k])/((varXX_s[i]*varXX_s[j])**(1/2))
			chan_corrXY[k][i][j]	=	XY_s[i][k]*np.conj(XY_s[j][k])/((varXY_s[i]*varXY_s[j])**(1/2))
			chan_corrYY[k][i][j]	=	YY_s[i][k]*np.conj(YY_s[j][k])/((varYY_s[i]*varYY_s[j])**(1/2))
			chan_corrYX[k][i][j]	=	YX_s[i][k]*np.conj(YX_s[j][k])/((varYX_s[i]*varYX_s[j])**(1/2))

#print len(varXX_s),varXX_s

for i in range(256):
	for j in range(i_n):
		chan_corrXX_m[i]	+=	chan_corrXX[j][i]
		chan_corrXY_m[i]	+=	chan_corrXY[j][i]
		chan_corrYY_m[i]	+=	chan_corrYY[j][i]
		chan_corrYX_m[i]	+=	chan_corrYX[j][i]
	chan_corrXX_m[i]	=	chan_corrXX_m[i]/i_n
	chan_corrXY_m[i]	=	chan_corrXY_m[i]/i_n
	chan_corrYY_m[i]	=	chan_corrYY_m[i]/i_n
	chan_corrYX_m[i]	=	chan_corrYX_m[i]/i_n

fig1	=	plt.figure(1)
plt.subplot(221)
plt.imshow(amplitude(chan_corrXX_m))
plt.xlabel('XX')
plt.subplot(222)
plt.imshow(amplitude(chan_corrXY_m))
plt.xlabel('XY')
plt.subplot(223)
plt.imshow(amplitude(chan_corrYY_m))
plt.xlabel('YY')
plt.subplot(224)
plt.imshow(amplitude(chan_corrYX_m))
plt.xlabel('YX')
fig1.savefig(file_name[:-4]+'_chan_corr.png')
plt.figure()
plt.plot(XX_m)
plt.plot(XY_m)
plt.plot(YY_m)
plt.plot(YX_m)
plt.figure()
plt.plot(varXX_s)
plt.plot(varXY_s)
plt.plot(varYY_s)
plt.plot(varYX_s)
plt.show()
