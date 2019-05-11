#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import isfile
import pyfftw
from astropy.table import Table
import time
from subprocess import call

t0	=	time.time()

########-----------------------------########
#			Function Definations
########-----------------------------########

def trans(a):
	if(np.abs(a[0]-a[1])==1):
		return 1
	else:
		return 0

def amplitude(a):
	a							=	(a.imag**2+a.real**2)**0.5
	return a

def mean(a):
	mean_a						=	np.zeros(256, dtype='complex')
	for i in range(256):
		for j in range(len(a)):
			mean_a[i]			+=	amplitude(a[j][i])
		mean_a[i]				=	mean_a[i]/len(a)
	return mean_a

def rms(a,mean_a):
	rms_a						=	np.zeros(256, dtype='complex')
	for i in range(256):
		for j in range(len(a)):
			rms_a[i]			+=	(amplitude(a[j][i])-mean_a[i])**2
		rms_a[i]				=	np.sqrt(rms_a[i]/len(a))
	return rms_a

def efficiency(a,exp):
	efficiency	=	np.zeros((len(a),256))
	for i in range(len(a)):
		mean_a		=	mean(a[i])
		rms_a		=	rms(a[i],mean_a)
		SNR			=	mean_a/rms_a
		efficiency[i]	=	SNR/exp
	return mean(efficiency)

def reject_matrix(eff,a):
	Mean_E		=	np.mean(eff)
	Sigma		=	np.std(eff)
	n_c			=	len(eff)
	R			=	np.zeros(n_c)
	for i in range(n_c):
		if(np.abs(eff[i]-Mean_E)>a*Sigma):
			R[i]	=	1
	return R

def reject(R1,R2,corr):
	n_c			=	len(R1)
	for j in range(n_c):
		if(R1[j]==1 or R2[j]==1):
			for i in range(len(corr)):
				corr[i][j]=0

def hilbert_delay(o_t,Z,z,padding,s_f,o_f,delay):
	for k in range(o_t):
		Z[k]						=	np.fft.fftshift(np.fft.ifft(np.concatenate((z[k],padding))))

		max_x						=	0
		for x in range(len(Z[k])):
			if(Z[k][x]==max(Z[k])):
				max_x				=	x

		n							=	max_x-o_f*512/2
		delay[k]					=	float(n)/float(o_f)/float(s_f)
		print 'Calculated Delay',k,':',delay[k]

def clear(a):
	for i in range(len(a)):
		a[i]	=	0

########-----------------------------########
#				Input Arguments
########-----------------------------########

if(len(sys.argv)==2):
	if(sys.argv[1]=='-d'):
		file_name_1	=	'raw_ch04_Cass_a_191_20180617_203306_000.mbr'
		i_t			=	4
		n_s			=	2
		mode		=	1
elif(len(sys.argv)==5):
	file_name_1	=	sys.argv[1]
	i_t			=	int(sys.argv[2])
	n_s			=	int(sys.argv[3])
	mode		=	int(sys.argv[4])
else:
	file_name_1	=	raw_input('File_1:')
	i_t			=	int(input('Integration Time(in s.): '))
	n_s			=	int(input('Confidence parameter: '))
	mode		=	int(input('Mode(GPS): '))

########-------------------------------########
#				Synchronisation
########-------------------------------########

f1	=	open(file_name_1,'rb')

if(mode==1):
	z_flag	=	1
	z_ctr1	=	-1
	while(z_flag):
		f1.read(28)
		p_s	=	int(f1.read(4).encode('hex'),16)
		f1.read(1024)
		z_ctr1	+=	1
		if(p_s==0):
			z_flag	=	0
else:
	z_ctr1	=	0

g_flag	=	0

f1.seek(-1030,2)
g_l1	=	int(f1.read(2).encode('hex'),16)

g		=	0
c		=	1
while(g!=(g_l1-1)):
	f1.seek(-(1056*c)-1030,2)
	g	=	int(f1.read(2).encode('hex'),16)
	c	+=	1
p_l1	=	int(f1.read(4).encode('hex'),16)

f1.seek(0)
for i in range(z_ctr1):
	f1.read(1056)
f1.read(26)

g1		=	int(f1.read(2).encode('hex'),16)
p_num1	=	np.zeros(2,dtype='int')

gps1	=	g1
while(gps1!=(g1+1)):
	f1.read(1054)
	gps1	=	int(f1.read(2).encode('hex'),16)

p_num1	=	int(f1.read(4).encode('hex'),16)

########----------------------------------########
#			Initialising Corr Variables
########----------------------------------########

o_t			=	g_l1-gps1
n_p			=	p_l1-p_num1+1
a_count		=	int(n_p/o_t)
#o_t		=	18		##############!!!!!!!!!!!!################
o_f			=	16								#oversampling_factor
s_f			=	33000000						#sampling_frequency
b_w			=	s_f/2							#bandwidth
padding		=	np.zeros(256+(512*(o_f-1)))
i_n			=	int(o_t/i_t)
a_count		=	int(i_t*a_count)
e_c			=	1000
a_t			=	int(a_count/e_c)
tbwp		=	e_c								#time-bandwidth-product
buff_1		=	np.zeros(1024,dtype='complex')
x_fftw_1	=	pyfftw.empty_aligned(512,dtype='complex')
y_fftw_1	=	pyfftw.empty_aligned(512,dtype='complex')
x_fftw_1	=	pyfftw.FFTW(x_fftw_1,x_fftw_1)
y_fftw_1	=	pyfftw.FFTW(y_fftw_1,y_fftw_1)
XY			=	np.zeros((i_n,256),dtype='complex')
X_A1_m		=	np.zeros((i_n,256),dtype='complex')
Y_A1_m		=	np.zeros((i_n,256),dtype='complex')
X_A1		=	np.zeros((i_n,a_t,256),dtype='complex')
Y_A1		=	np.zeros((i_n,a_t,256),dtype='complex')
Z_XY		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_X1		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_Y1		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_XY_C		=	np.zeros((i_n,o_f*512), dtype='complex')
delay_XY	=	np.zeros(i_n)

l_t     =   [i for i in range(256)]
for j in range(128):
    l_t[128+j]  -=  256

########-------------------------------------########
#			Calculating Cross Correlation
########-------------------------------------########

print 'Starting Corr'

j	=	0
for i in range(i_n):
	for k in range(a_t):
		while(j<(i+1)*(k+1)*1024*e_c):
			buff_1[j%1024]	=	l_t[int(str(f1.read(1)).encode('hex'),16)]
			j				+=	1
			if(j%1024==0):
				X_1			=	x_fftw_1(buff_1[::2])
				Y_1			=	y_fftw_1(buff_1[1::2])
				X_A1[i][k]	+=	X_1[:256]*np.conj(X_1[:256])
				Y_A1[i][k]	+=	Y_1[:256]*np.conj(Y_1[:256])
				XY[i]		+=	X_1[:256]*np.conj(Y_1[:256])
				f1.read(32)
		X_A1[i][k][0]	=	X_A1[i][k][1]
		Y_A1[i][k][0]	=	Y_A1[i][k][1]
		X_A1[i][k]		=	X_A1[i][k]/float(e_c)
		Y_A1[i][k]		=	Y_A1[i][k]/float(e_c)
	XY[i][0]	=	XY[i][1]
	XY[i]		=	XY[i]/float(a_count)
	print i

f1.close()

########------------------------------########
#			Calculating Efficiency
########------------------------------########

print 'Calculating Efficiency'

SNR_exp			=	float(np.sqrt(tbwp))

efficiency_x1	=	efficiency(X_A1,SNR_exp)
efficiency_y1	=	efficiency(Y_A1,SNR_exp)

########-----------------------------########
#			Clearing Bad Channels
########-----------------------------########

print 'Clearing Bad Channels'

R_X1		=	reject_matrix(efficiency_x1,n_s)
R_Y1		=	reject_matrix(efficiency_y1,n_s)

reject(R_X1,R_Y1,XY)

########-----------------------------########
#			  Hilbert Transform
########-----------------------------########

print 'Performing Hilbert Transform'

hilbert_delay(i_n,Z_XY,XY,padding,s_f,o_f,delay_XY)

for i in range(i_n):
	X_A1_m[i]	=	mean(X_A1[i])
	Y_A1_m[i]	=	mean(Y_A1[i])

Z_X1	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(X_A1_m),padding))))
Z_Y1	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(Y_A1_m),padding))))
Z_XY_C	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(XY),padding))))

########-----------------------########
#			Saving Results
########-----------------------########

print 'Saving Results'

XY_t	=	Table(XY)
XY_t.write(file_name_1[:-4]+'_cross_corr(6).txt',format='ascii.no_header',overwrite='true')

print 'Total time elapsed:',time.time()-t0

########-------------------------########
#			Plotting Results
########-------------------------########

print 'Plotting Results'

#fig		=	plt.figure(1)
#plt.plot(amplitude(XY))
#plt.ylabel('XY')
#fig.savefig(file_name_1[:-4]+'_cross_corr.png')

fig11	=	plt.figure(11)
plt.plot(amplitude(mean(XY)))
plt.ylabel('XY')
fig11.savefig(file_name_1[:-4]+'_cross_corr(6).png')

#fig2	=	plt.figure(2)
#plt.imshow(amplitude(XY),interpolation='nearest',aspect='auto')
#plt.ylabel('XY')
#fig2.savefig(file_name_1[:-4]+'dynamic_spectra.png')

fig21	=	plt.figure(21)
plt.imshow(amplitude(XY),interpolation='bilinear',aspect='auto')
plt.ylabel('XY')
fig21.savefig(file_name_1[:-4]+'dynamic_spectra(6).png')

fig3	=	plt.figure(3)
plt.subplot(211)
plt.plot(amplitude(X_A1_m))
plt.ylabel('X_A1')
plt.subplot(212)
plt.plot(amplitude(Y_A1_m))
plt.ylabel('Y_A1')
fig3.savefig(file_name_1[:-4]+'_auto_corr(6).png')

#fig31	=	plt.figure(31)
#plt.subplot(211)
#plt.imshow(amplitude(X_A1_m),interpolation='nearest',aspect='auto')
#plt.ylabel('X1')
#plt.subplot(212)
#plt.imshow(amplitude(Y_A1_m),interpolation='nearest',aspect='auto')
#plt.ylabel('Y1')
#fig31.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr_im.png')

fig32	=	plt.figure(32)
plt.subplot(211)
plt.imshow(amplitude(X_A1_m),interpolation='bilinear',aspect='auto')
plt.ylabel('X1')
plt.subplot(212)
plt.imshow(amplitude(Y_A1_m),interpolation='bilinear',aspect='auto')
plt.ylabel('Y1')
fig32.savefig(file_name_1[:-4]+'_power_spectra(auto_corr)(6).png')

fig4	=	plt.figure(4)
plt.subplot(211)
plt.plot(efficiency_x1)
plt.plot([np.mean(efficiency_x1)-n_s*np.std(efficiency_x1) for i in range(len(efficiency_x1))])
plt.plot([np.mean(efficiency_x1)+n_s*np.std(efficiency_x1) for i in range(len(efficiency_x1))])
plt.ylabel('Efficiency_X1')
plt.subplot(212)
plt.plot(efficiency_y1)
plt.plot([np.mean(efficiency_y1)-n_s*np.std(efficiency_y1) for i in range(len(efficiency_y1))])
plt.plot([np.mean(efficiency_y1)+n_s*np.std(efficiency_y1) for i in range(len(efficiency_y1))])
plt.ylabel('Efficiency_Y1')
fig4.savefig(file_name_1[:-4]+'_efficiency(6).png')

fig5	=	plt.figure(5)
plt.subplot(211)
plt.plot(Z_X1)
plt.xlabel('Hilbert Transform: X')
plt.subplot(211)
plt.plot(Z_Y1)
plt.xlabel('Hilbert Transform: Y')
fig5.savefig(file_name_1[:-4]+'_hilbert(6).png')

fig51	=	plt.figure(51)
plt.plot(amplitude(Z_XY_C))
plt.xlabel('Hilbert Transform: XY')
fig51.savefig(file_name_1[:-4]+'_hilbert2(6).png')

fig6	=	plt.figure(6)
plt.plot(delay_XY*(10**9))
plt.ylabel('Delay(ns)[XY]')
plt.xlabel('Time(s)')
fig6.savefig(file_name_1[:-4]+'_delay(6).png')

fig12	=	plt.figure(12)
plt.plot(amplitude(XY[5]))
plt.ylabel('XY')
fig12.savefig(file_name_1[:-4]+'_cross_corr[5](6).png')

