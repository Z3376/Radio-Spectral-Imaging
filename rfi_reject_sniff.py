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

def reject_matrix(eff,a):
	#Mean_E		=	np.mean(eff)
	#Sigma		=	np.std(eff)
	n_c			=	len(eff)
	R			=	np.zeros(n_c)
	for i in range(256):
		if(np.abs(eff[i]-np.mean(eff))>a*np.std(eff)):
			R[i]	=	1
	return R

file_name	=	sys.argv[1]
if(len(sys.argv)==3):
	n		=	sys.argv[2]
else:
	n		=	2	

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

das	=	7
f.read(7)
d	=	int(f.read(1))
f.read(24)
b_c			=	0
for j in range(a_t):
	while(b_c<(j+1)*e_c*1024):
		if(d==das):
			buff[b_c%1024]			=	l_t[int(f.read(1).encode('hex'),16)]
			b_c						+=	1
			if(b_c%1024==0):
				X1					=	x_fftw(buff[::2])
				Y1					=	y_fftw(buff[1::2])
				X_A[j]				+=	X1[:256]*np.conj(X1[:256])
				Y_A[j]				+=	Y1[:256]*np.conj(Y1[:256])
				f.read(7)
				d	=	int(f.read(1))
				f.read(24)
		else:
			b_c	+=	1024
			f.read(1031)
			d	=	int(f.read(1))
			f.read(24)
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


ftxt		=	open(file_name[:-8]+'_reject_matrix.txt','w+')

for i in range(len(R_X)):
	ftxt.write(str(R_X[i]))
	ftxt.write('	')
	ftxt.write(str(R_Y[i]))
	ftxt.write('\n')

ftxt.close()

fig	=	plt.figure()
plt.subplot(211)
plt.plot(efficiency_x)
#plt.plot([np.mean(efficiency_x)-n*np.std(efficiency_x) for i in range(len(efficiency_x))])
#plt.plot([np.mean(efficiency_x)+n*np.std(efficiency_x) for i in range(len(efficiency_x))])
plt.ylabel('Efficiency_X')
plt.subplot(212)
plt.plot(efficiency_y)
#plt.plot([np.mean(efficiency_y)-n*np.std(efficiency_y) for i in range(len(efficiency_y))])
#plt.plot([np.mean(efficiency_y)+n*np.std(efficiency_y) for i in range(len(efficiency_y))])
plt.ylabel('Efficiency_Y')
fig.savefig(file_name[:-8]+'_reject_matrix.png')

plt.show()

