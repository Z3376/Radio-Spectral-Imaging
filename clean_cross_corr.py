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

def reject_matrix(XX_m,n_s):
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

file_name	=	'ch00_CYG-A_20161215_101054_000_cross_corr.txt'
#file_name	=	raw_input('File: ')
n_s		=	2
fopen	=	open(file_name,'r')
size	=	os.path.getsize(file_name)
s1		=	float(size)
ARR		=	[]

corr_read(size,fopen,ARR)

l	= int(float(len(ARR)))

n	=	int(l/4)
i_n	=	int(n/256)

XY		=	np.zeros((i_n,256),dtype='complex')

for i in range(i_n):
		for j in range(256):
			XY[i][j]	=	ARR[(i*256)+j]

XY_m	=	mean(XY)

R_XY	=	reject_matrix(XY_m,n_s)

reject(R_XY,XY)

fig1	=	plt.figure(1)
plt.imshow(amplitude(XY.T),interpolation='bilinear',aspect='auto')
plt.ylabel('XY')
fig1.savefig(file_name[:-4]+'_im(clean).png')

fig11	=	plt.figure(11)
plt.plot(amplitude(XY_m))
plt.plot([np.mean(XY_m)-n_s*np.std(XY_m) for i in range(len(XY_m))])
plt.plot([np.mean(XY_m)+n_s*np.std(XY_m) for i in range(len(XY_m))])
plt.ylabel('XY')
fig11.savefig(file_name[:-4]+'1(clean).png')

plt.figure()
plt.plot(R_XY)

plt.show()
