#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os

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
			if(not(i<int(s1))):
				break
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

file_name	=	'ch00_CASA_20170602_234118_000_cross_corr.txt'
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

fig1	=	plt.figure(1)
plt.subplot(221)
plt.imshow(amplitude(XX.T),interpolation='bilinear',aspect='auto')
plt.ylabel('XX')
plt.subplot(222)
plt.imshow(amplitude(XY.T),interpolation='bilinear',aspect='auto')
plt.ylabel('XY')
plt.subplot(223)
plt.imshow(amplitude(YY.T),interpolation='bilinear',aspect='auto')
plt.ylabel('YY')
plt.subplot(224)
plt.imshow(amplitude(YX.T),interpolation='bilinear',aspect='auto')
plt.ylabel('YX')
fig1.savefig(file_name[:-4]+'_im(clean).png')

fig11	=	plt.figure(11)
plt.subplot(221)
plt.plot(amplitude(XX_m))
plt.plot([np.mean(XX_m)-n_s*np.std(XX_m) for i in range(len(XX_m))])
plt.plot([np.mean(XX_m)+n_s*np.std(XX_m) for i in range(len(XX_m))])
plt.ylabel('XX')
plt.subplot(222)
plt.plot(amplitude(XY_m))
plt.plot([np.mean(XY_m)-n_s*np.std(XY_m) for i in range(len(XY_m))])
plt.plot([np.mean(XY_m)+n_s*np.std(XY_m) for i in range(len(XY_m))])
plt.ylabel('XY')
plt.subplot(223)
plt.plot(amplitude(YY_m))
plt.plot([np.mean(YY_m)-n_s*np.std(YY_m) for i in range(len(YY_m))])
plt.plot([np.mean(YY_m)+n_s*np.std(YY_m) for i in range(len(YY_m))])
plt.ylabel('YY')
plt.subplot(224)
plt.plot(amplitude(YX_m))
plt.plot([np.mean(YX_m)-n_s*np.std(YX_m) for i in range(len(YX_m))])
plt.plot([np.mean(YX_m)+n_s*np.std(YX_m) for i in range(len(YX_m))])
plt.ylabel('YX')
fig11.savefig(file_name[:-4]+'1(clean).png')

plt.figure()
plt.subplot(221)
plt.plot(R_XX)
plt.subplot(222)
plt.plot(R_XY)
plt.subplot(223)
plt.plot(R_YY)
plt.subplot(224)
plt.plot(R_YX)

plt.show()
