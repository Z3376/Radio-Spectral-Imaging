#!/usr/bin/python

import sys
from os.path import getsize
import numpy as np
import matplotlib.pyplot as plt
import pyfftw

def amplitude(a):
	a							=	(a.imag**2+a.real**2)**0.5
	return a

def mean(a):
	mean_a						=	np.zeros(256, dtype='complex')
	for i in range(256):
		for j in range(len(a)):
			mean_a[i]			+=	a[j][i]
		mean_a[i]				=	mean_a[i]/len(a)
	return mean_a

def rms(a,mean_a):
	rms_a						=	np.zeros(256, dtype='complex')
	for i in range(256):
		for j in range(len(a)):
			rms_a[i]			+=	(a[j][i]-mean_a[i])**2
		rms_a[i]				=	np.sqrt(rms_a[i]/len(a))
	return rms_a

def efficiency(a,exp):
	mean_a		=	mean(a)
	rms_a		=	rms(a,mean_a)
	SNR			=	mean_a/rms_a
	efficiency	=	SNR/exp
	return efficiency

def reject(R1,corr):
	n_c			=	len(R1)
	for j in range(n_c):
		if(R1[j]==1):
			for i in range(len(corr)):
				corr[i][j]=0


file_name	=	sys.argv[1]
if(len(sys.argv)==3):
	n		=	sys.argv[2]
else:
	n		=	3	

f			=	open(file_name,'rb')

s			=	getsize(file_name)

t_p			=	int(s/1056)

e_c			=	1000
a_t			=	t_p/e_c

buff		=	np.zeros(1024,dtype='complex')
x_fftw		=	pyfftw.empty_aligned(512,dtype='complex')
y_fftw		=	pyfftw.empty_aligned(512,dtype='complex')
X_A			=	np.zeros((a_t,256),dtype='complex')
Y_A			=	np.zeros((a_t,256),dtype='complex')

x_fftw		=	pyfftw.FFTW(x_fftw,x_fftw)						#FFTW_matrix
y_fftw		=	pyfftw.FFTW(y_fftw,y_fftw)

l_t     =   [i for i in range(256)]
for j in range(128):
    l_t[128+j]  -=  256

f.read(32)
b_c			=	0
for j in range(a_t):
	while(b_c<(j+1)*e_c*1024):
		buff[b_c%1024]			=	l_t[int(f.read(1).encode('hex'),16)]
		b_c						+=	1
		if(b_c%1024==0):
			X1					=	x_fftw(buff[::2])
			Y1					=	y_fftw(buff[1::2])
			X_A[j]				+=	X1[:256]*np.conj(X1[:256])
			Y_A[j]				+=	Y1[:256]*np.conj(Y1[:256])
			f.read(32)
	X_A[j][0]					=	X_A[j][1]
	Y_A[j][0]					=	Y_A[j][1]
	X_A[j]						=	X_A[j]/float(e_c)
	Y_A[j]						=	Y_A[j]/float(e_c)
	print j

f.close()

print 'Calculating Efficiency'

tbwp		=	e_c
SNR_exp		=	float(np.sqrt(tbwp))

efficiency_x	=	efficiency(X_A,SNR_exp)
efficiency_y	=	efficiency(Y_A,SNR_exp)

R_X		=	np.zeros(len(efficiency_x))
R_Y		=	np.zeros(len(efficiency_y))

Mean_E	=	np.mean(efficiency_x)
Sigma	=	np.std(efficiency_x)
for i in range(len(efficiency_x)):
	if(np.abs(efficiency_x[i]-Mean_E)>2*Sigma):
		R_X[i]	=	1

Mean_E	=	np.mean(efficiency_y)
Sigma	=	np.std(efficiency_y)
for i in range(len(efficiency_y)):
	if(np.abs(efficiency_y[i]-Mean_E)>2*Sigma):
		R_Y[i]	=	1

reject(R_X,X_A)
reject(R_Y,Y_A)

X_A_m	=	mean(X_A)
Y_A_m	=	mean(Y_A)

X_A_m	=	X_A_m-np.mean(X_A_m)
Y_A_m	=	Y_A_m-np.mean(Y_A_m)

X_A_s	=	np.zeros((256,a_t),dtype='complex')
Y_A_s	=	np.zeros((256,a_t),dtype='complex')

for i in range(256):
	for j in range(len(X_A)):
		X_A_s[i][j]	=	X_A[j][i]
		Y_A_s[i][j]	=	Y_A[j][i]


varX_s	=	np.zeros(256,dtype='complex')
varY_s	=	np.zeros(256,dtype='complex')

for i in range(256):
	for j in range(a_t):
		varX_s[i]	+=	(X_A_s[i][j]-np.mean(X_A_s[i]))**2
		varY_s[i]	+=	(Y_A_s[i][j]-np.mean(Y_A_s[i]))**2
	varX_s[i]	=	np.sqrt((varX_s[i])/a_t)
	varY_s[i]	=	np.sqrt((varY_s[i])/a_t)
	
chan_corr_X	=	np.zeros((a_t,256,256),dtype='complex')
chan_corr_Y	=	np.zeros((a_t,256,256),dtype='complex')
chan_corrX_m	=	np.zeros((256,256),dtype='complex')
chan_corrY_m	=	np.zeros((256,256),dtype='complex')

for k in range(a_t):
	for i in range(256):
		for j in range(256):
			chan_corr_X[k][i][j]	=	X_A_s[i][k]*np.conj(X_A_s[j][k])/(varX_s[i]*varX_s[j]**(1/2))
			chan_corr_Y[k][i][j]	=	Y_A_s[i][k]*np.conj(Y_A_s[j][k])/(varY_s[i]*varY_s[j]**(1/2))

for i in range(256):
	for j in range(a_t):
		chan_corrX_m[i]	+=	chan_corr_X[j][i]
		chan_corrY_m[i]	+=	chan_corr_Y[j][i]
	chan_corrX_m[i]	=	chan_corrX_m[i]/a_t
	chan_corrY_m[i]	=	chan_corrY_m[i]/a_t

fig	=	plt.figure()
plt.subplot(211)
plt.imshow(amplitude(chan_corrX_m))
plt.subplot(212)
plt.imshow(amplitude(chan_corrY_m))
fig.savefig(file_name[:-8]+'_chan_corr_2.png')

plt.show()

