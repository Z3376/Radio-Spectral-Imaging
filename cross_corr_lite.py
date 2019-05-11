#!/usr/bin/python

from os.path import getsize
from os.path import isfile
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
import pyfftw
from astropy.table import Table
from progress.bar import Bar

t0	=	time.time()

########-----------------------------########
#			Function Definations
########-----------------------------########

def amplitude(a):
	a							=	(a.imag**2+a.real**2)**0.5
	return a

def mean(a):
	mean_a						=	np.zeros(256, dtype='complex')
	for i in range(256):
		for j in range(a_count):
			mean_a[i]			+=	a[j][i]
		mean_a[i]				=	mean_a[i]/a_count
	return mean_a

def rms(a,mean_a):
	rms_a						=	np.zeros(256, dtype='complex')
	for i in range(256):
		for j in range(a_count):
			rms_a[i]			+=	(a[j][i]-mean_a[i])**2
		rms_a[i]				=	np.sqrt(rms_a[i]/a_count)
	return rms_a

########-----------------------------########
#				Input Arguments
########-----------------------------########

if(len(sys.argv)==2 and sys.argv[1]=='-d'):
	file_name					=	'ch00_CYG-A_20161215_101054_000.mbr'
	avg							=	64
	das_id						=	3
elif(len(sys.argv)==4):
	file_name					=	sys.argv[1]
	avg							=	int(sys.argv[2])
	das_id						=	sys.argv[3]
elif(len(sys.argv)==3):
	file_name					=	sys.argv[1]
	avg							=	int(sys.argv[2])
	das_id						=	raw_input('Das_id: ')
elif(len(sys.argv)==2):
	file_name					=	sys.argv[1]
	avg							=	int(raw_input('No. of Packets to be averaged: '))
	das_id						=	raw_input('Das_id: ')
elif(len(sys.argv)==1):
	file_name					=	raw_input('.mbr File: ')
	avg							=	int(raw_input('No. of Packets to be averaged: '))
	das_id						=	raw_input('Das_id: ')
else:
	print('Too many arguments. Enter onyl .mbr File name and No. of Packets to be averaged.')
	exit()

########-----------------------------########
#			Initialising Variables
########-----------------------------########

size							=	getsize(file_name)
t_p								=	size/1056										#total_packets
b_c								=	0												#byte_counter
#s_c								=	0											#sample_counter
o_f								=	64												#oversampling_factor
s_f								=	33000000										#sampling_frequency
b_w								=	s_f/2											#bandwidth
tbwp							=	avg												#time-bandwidth-product
buff							=	np.zeros(1024,dtype='complex')
x_fftw							=	pyfftw.empty_aligned(512,dtype='complex')
y_fftw							=	pyfftw.empty_aligned(512,dtype='complex')
padding							=	np.zeros(256+(512*(o_f-1)))

x_fftw							=	pyfftw.FFTW(x_fftw,x_fftw)						#FFTW_matrix
y_fftw							=	pyfftw.FFTW(y_fftw,y_fftw)

l_t								=	[i for i in range(256)]							#look-up_table
for j in range(128):
	l_t[j+128]					-=	256

########-----------------------------########
#		Calculating Cross Correlation
########-----------------------------########
print t_p
print('Reading File')

fopen							=	open(file_name,'rb')							#reading_mbr_file_starts
fopen.read(7)

p_d		=	np.zeros(t_p)

t_p_d	=	0
for j in range(t_p):
	if(fopen.read(1)==str(das_id)):
		t_p_d	+=	1
		p_d[j]	=	1
	fopen.read(1055)

a_count							=	t_p_d/avg
z								=	np.zeros((a_count,256), dtype='complex')
x								=	np.zeros((a_count,256), dtype='complex')
y								=	np.zeros((a_count,256), dtype='complex')
Z								=	np.zeros((a_count,o_f*512), dtype='complex')
delay							=	np.zeros(a_count)							

fopen.seek(0,0)
fopen.read(32)

i	=	0
q	=	0
bar	=	Bar('Calculating Cross Correlation',max=a_count)
for j in range(a_count):
	while(i<(j+1)*avg*1024):
		buff[i%1024]						=	l_t[int(str(fopen.read(1)).encode('hex'),16)]
		i	+=	1
		q	+=	1
		if(i%1024==0 and p_d[int(q/1024)-1]==1):
			X						=	x_fftw(buff[0::2])
			Y						=	y_fftw(buff[1::2])
			z[j]					+=	X[:256]*np.conj(Y[:256])
			x[j]					+=	X[:256]*np.conj(X[:256])
			y[j]					+=	Y[:256]*np.conj(Y[:256])
			fopen.read(32)
		elif(i%1024==0):
			fopen.read(32)
			i	-=	1024
	z[j][0]						=	z[j][1]
	x[j][0]						=	x[j][1]
	y[j][0]						=	y[j][1]
	z[j]						=	z[j]/float(avg)/np.sqrt(x[j]*y[j])
	bar.next()
bar.finish()
fopen.close()

########-----------------------------########
#			Calculating Efficiency
########-----------------------------########

print('Calculating Efficiency')

mean_x_o						=	mean(x)
mean_y_o						=	mean(y)

rms_x_o							=	rms(x,mean_x_o)
rms_y_o							=	rms(y,mean_y_o)

SNR_x_obs						=	mean_x_o/rms_x_o
SNR_y_obs						=	mean_y_o/rms_y_o
SNR_exp							=	float(np.sqrt(tbwp))

efficiency_x					=	SNR_x_obs/SNR_exp
efficiency_y					=	SNR_y_obs/SNR_exp

########-----------------------------########
#			Clearing Bad Channels
########-----------------------------########

print('Clearing Bad Channels')

for i in range(256):
	if(np.abs(efficiency_x[i]-np.mean(efficiency_x))>2*np.std(efficiency_x) or np.abs(efficiency_y[i]-np.mean(efficiency_y))>2*np.std(efficiency_y)):
		for j in range(a_count):
			z[j][i]				=	0

print 'Saving Results'

z_t	=	Table(z)
z_t.write(file_name[:-4]+'_cross_corr.txt',format='ascii.no_header',overwrite='true')

print "Total Time Elapsed:", time.time()-t0

########-----------------------------########
#			  Plotting Results
########-----------------------------########

fig1=plt.figure(1)
plt.plot(efficiency_x)
plt.plot(efficiency_y)
plt.xlabel('Frequency Channels')
plt.ylabel('Efficiency')
fig1.savefig(file_name[:-4]+'_eff.png')
fig6=plt.figure(6)
plt.imshow(amplitude(z.T),interpolation='bilinear',aspect='auto')
plt.colorbar()
plt.xlabel('Time(s)')
plt.ylabel('Frequency Channels')
fig6.savefig(file_name[:-4]+'_cross_corr.png')
fig8=plt.figure(8)
plt.plot(amplitude(z[1]))
plt.plot(z[1].real)
plt.plot(z[1].imag)
fig8.savefig(file_name[:-4]+'_cross_corr[5].png')
fig2=plt.figure(2)
plt.imshow(z.real.T,interpolation='bilinear',aspect='auto')
plt.colorbar()
plt.xlabel('Time(s)')
plt.ylabel('Frequency Channels')
fig3.savefig(file_name[:-4]+'real_cross_corr.png')
fig3=plt.figure(3)
plt.imshow(z.imag.T,interpolation='bilinear',aspect='auto')
plt.colorbar()
plt.xlabel('Time(s)')
plt.ylabel('Frequency Channels')
fig3.savefig(file_name[:-4]+'imag_cross_corr.png')

exit()
