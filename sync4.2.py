#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
#import os
import pyfftw
from astropy.table import Table
import time

t0	=	time.time()

def trans(a):
    if(np.abs(a[0]-a[1])==1):
        return 1
    else:
        return 0

def amplitude(a):
	a							=	np.sqrt(a.imag**2+a.real**2)
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

if(sys.argv[1]=='-d'):
	file_name_1	=	'raw_ch07_B_P_20170714_042859_019.mbr'
	file_name_2	=	'raw_ch06_B_P_20170714_042859_019.mbr'
	i_t			=	2
	n_s			=	1
	mode		=	2
else:
	file_name_1	=	raw_input('File_1:')
	file_name_2	=	raw_input('File_2:')
	i_t			=	int(input('Integration Time(in s.): '))
	n_s			=	int(input('Confidence parameter: '))
	mode		=	int(input('Mode(GPS): '))

f1	=	open(file_name_1,'rb')
f2	=	open(file_name_2,'rb')

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
		
	z_flag	=	1
	z_ctr2	=	-1
	while(z_flag):
		p_s	=	int(f2.read(4).encode('hex'),16)
		f2.read(1024)
		z_ctr2	+=	1
		if(p_s==0):
			z_flag	=	0
else:
	z_ctr1	=	0
	z_ctr2	=	0
f1.seek(-1030,2)
g_l1	=	int(str(f1.read(2)).encode('hex'),16)

f2.seek(-1030,2)
g_l2	=	int(str(f2.read(2)).encode('hex'),16)

g		=	0
c		=	1
while(g!=(g_l1-1)):
	f1.seek(-(1056*c)-1030,2)
	g	=	int(str(f1.read(2)).encode('hex'),16)
	c	+=	1
p_l1	=	int(str(f1.read(4)).encode('hex'),16)

g		=	0
c		=	1
while(g!=(g_l2-1)):
	f2.seek(-(1056*c)-1030,2)
	g	=	int(str(f2.read(2)).encode('hex'),16)
	c	+=	1
p_l2	=	int(str(f2.read(4)).encode('hex'),16)

f1.seek(0)
f2.seek(0)

for i in range(z_ctr1):
	f1.read(1056)
for i in range(z_ctr2):
	f2.read(1056)

f1.read(26)
f2.read(26)

g1		=	int(str(f1.read(2)).encode('hex'),16)
g2		=	int(str(f2.read(2)).encode('hex'),16)

gps1	=	0
while(gps1!=(g1+1)):
	f1.read(1054)
	gps1	=	int(str(f1.read(2)).encode('hex'),16)

gps2	=	0
while(gps2!=(g2+1)):
	f2.read(1054)
	gps2	=	int(str(f2.read(2)).encode('hex'),16)

if(gps1==gps2):
	p_num1	=	int(str(f1.read(4)).encode('hex'),16)
	p_num2	=	int(str(f2.read(4)).encode('hex'),16)
	if(g_l1<g_l2):
		o_t		=	g_l1-gps1
		n_p		=	p_l1-p_num1+1
	else:
		o_t		=	g_l2-gps1
		n_p		=	p_l2-p_num2+1

elif(gps1>gps2):
	while(gps2!=gps1):
		f2.read(1054)
		gps2	=	int(str(f2.read(2)).encode('hex'),16)
	p_num1	=	int(str(f1.read(4)).encode('hex'),16)
	p_num2	=	int(str(f2.read(4)).encode('hex'),16)
	o_t		=	g_l2-gps1
	n_p		=	p_l2-p_num2+1

else:
	while(gps1!=gps2):
		f1.read(1054)
		gps1	=	int(str(f1.read(2)).encode('hex'),16)
	p_num1	=	int(str(f1.read(4)).encode('hex'),16)
	p_num2	=	int(str(f2.read(4)).encode('hex'),16)
	o_t		=	g_l1-gps1
	n_p		=	p_l1-p_num1+1

off		=	p_num1-p_num2

a_count	=	int(n_p/o_t)

#o_t			=	4		##############!!!!!!!!!!!!################

tbwp		=	a_count							#time-bandwidth-product
o_f			=	16								#oversampling_factor
s_f			=	33000000						#sampling_frequency
b_w			=	s_f/2							#bandwidth
padding		=	np.zeros(256+(512*(o_f-1)))

i_n			=	int(o_t/i_t)
a_count		=	int(i_t*a_count)

e_c			=	1000
a_t			=	int(a_count/e_c)


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
XX          =   np.zeros((i_n,256),dtype='complex')
YY          =   np.zeros((i_n,256),dtype='complex')
XY          =   np.zeros((i_n,256),dtype='complex')
YX          =   np.zeros((i_n,256),dtype='complex')
X_A1_m		=   np.zeros((i_n,256),dtype='complex')
Y_A1_m		=   np.zeros((i_n,256),dtype='complex')
X_A2_m		=   np.zeros((i_n,256),dtype='complex')
Y_A2_m		=   np.zeros((i_n,256),dtype='complex')
X_A1		=   np.zeros((i_n,a_t,256),dtype='complex')
Y_A1		=   np.zeros((i_n,a_t,256),dtype='complex')
X_A2		=   np.zeros((i_n,a_t,256),dtype='complex')
Y_A2		=   np.zeros((i_n,a_t,256),dtype='complex')
Z_XX		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_XY		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_YY		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_YX		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_X1		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_X2		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_Y1		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_Y2		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_XX_C		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_XY_C		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_YY_C		=	np.zeros((i_n,o_f*512), dtype='complex')
Z_YX_C		=	np.zeros((i_n,o_f*512), dtype='complex')
delay_XX	=	np.zeros(i_n)
delay_XY	=	np.zeros(i_n)
delay_YY	=	np.zeros(i_n)
delay_YX	=	np.zeros(i_n)

l_t     =   [i for i in range(256)]
for j in range(128):
    l_t[128+j]  -=  256

print 'Starting Corr'

j	=	0
for i in range(i_n):
	for k in range(a_t):
		while(j<(i+1)*(k+1)*1024*e_c):
			buff_1[j%1024]	=	l_t[int(str(f1.read(1)).encode('hex'),16)]
			buff_2[j%1024]	=	l_t[int(str(f2.read(1)).encode('hex'),16)]
			j				+=	1
			if(j%1024==0):
				X_1			=	x_fftw_1(buff_1[::2])
				Y_1			=	y_fftw_1(buff_1[1::2])
				X_2			=	x_fftw_2(buff_2[::2])
				Y_2			=	y_fftw_2(buff_2[1::2])
				X_A1[i][k]	+=	X_1[:256]*np.conj(X_1[:256])
				Y_A1[i][k]	+=	Y_1[:256]*np.conj(Y_1[:256])
				X_A2[i][k]	+=	X_2[:256]*np.conj(X_2[:256])
				Y_A2[i][k]	+=	Y_2[:256]*np.conj(Y_2[:256])
				XX[i]		+=	X_1[:256]*np.conj(X_2[:256])
				YY[i]		+=	Y_1[:256]*np.conj(Y_2[:256])
				XY[i]		+=	X_1[:256]*np.conj(Y_2[:256])
				YX[i]		+=	Y_1[:256]*np.conj(X_2[:256])
				f1.read(28)
				f2.read(28)
				p_num1		=	int(str(f1.read(4)).encode('hex'),16)
				p_num2		=	int(str(f2.read(4)).encode('hex'),16)
				if(p_num1-p_num2!=off):
					print 'slip at:',p_num1,p_num2
				while(p_num1-p_num2!=off):
					if(p_num1-p_num2>off):
						j	+=	1024
						f2.read(1052)
						p_num2	=	int(str(f2.read(4)).encode('hex'),16)				
					else:
						j	+=	1024
						f1.read(1052)
						p_num1	=	int(str(f1.read(4)).encode('hex'),16)
		X_A1[i][k][0]	=	X_A1[i][k][1]
		X_A2[i][k][0]	=	X_A2[i][k][1]
		Y_A1[i][k][0]	=	Y_A1[i][k][1]
		Y_A2[i][k][0]	=	Y_A2[i][k][1]
		X_A1[i][k]		=	X_A1[i][k]/float(e_c)
		X_A2[i][k]		=	X_A2[i][k]/float(e_c)
		Y_A1[i][k]		=	Y_A1[i][k]/float(e_c)
		Y_A2[i][k]		=	Y_A2[i][k]/float(e_c)
	XX[i][0]	=	XX[i][1]
	XY[i][0]	=	XY[i][1]
	YY[i][0]	=	YY[i][1]
	YX[i][0]	=	YX[i][1]
	XX[i]		=	XX[i]/float(a_count)
	XY[i]		=	XY[i]/float(a_count)
	YY[i]		=	YY[i]/float(a_count)
	YX[i]		=	YX[i]/float(a_count)
	print i

f1.close()
f2.close()

print 'Calculating Efficiency'

SNR_exp			=	float(np.sqrt(tbwp))

efficiency_x1	=	efficiency(X_A1,SNR_exp)
efficiency_x2	=	efficiency(X_A2,SNR_exp)
efficiency_y1	=	efficiency(Y_A1,SNR_exp)
efficiency_y2	=	efficiency(Y_A2,SNR_exp)

print 'Rejecting Bad Channels'

R_X1		=	reject_matrix(efficiency_x1,n_s)
R_X2		=	reject_matrix(efficiency_x2,n_s)
R_Y1		=	reject_matrix(efficiency_y1,n_s)
R_Y2		=	reject_matrix(efficiency_y2,n_s)

reject(R_X1,R_X2,XX)
reject(R_X1,R_Y2,XY)
reject(R_Y1,R_Y2,YY)
reject(R_Y1,R_X2,YX)

print 'Performing Hilbert Transform'

hilbert_delay(i_n,Z_XX,XX,padding,s_f,o_f,delay_XX)
hilbert_delay(i_n,Z_XY,XY,padding,s_f,o_f,delay_XY)
hilbert_delay(i_n,Z_YY,YY,padding,s_f,o_f,delay_YY)
hilbert_delay(i_n,Z_YX,YX,padding,s_f,o_f,delay_YX)

for i in range(i_n):
	X_A1_m[i]	=	mean(X_A1[i])
	X_A2_m[i]	=	mean(X_A2[i])
	Y_A1_m[i]	=	mean(Y_A1[i])
	Y_A2_m[i]	=	mean(Y_A2[i])


Z_X1	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(X_A1_m),padding))))
Z_X2	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(X_A2_m),padding))))
Z_Y1	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(Y_A1_m),padding))))
Z_Y2	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(Y_A2_m),padding))))

Z_XX_C	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(XX),padding))))
Z_XY_C	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(XY),padding))))
Z_YY_C	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(YY),padding))))
Z_YX_C	=	np.fft.fftshift(np.fft.ifft(np.concatenate((mean(YX),padding))))

print 'Saving Results'

XX_t	=	Table(XX)
XX_t.write(file_name_1[:-4]+'=X='+file_name_2[:-4]+'(4.2).txt',format='ascii.no_header',overwrite='true')
fw		=	open(file_name_1[:-4]+'=X='+file_name_2[:-4]+'(4.2).txt',mode='a')
fw.seek(0,2)
XY_t	=	Table(XY)
XY_t.write(fw,format='ascii.no_header')
YY_t	=	Table(YY)
YY_t.write(fw,format='ascii.no_header')
YX_t	=	Table(YX)
YX_t.write(fw,format='ascii.no_header')
fw.close()

print 'Total time elapsed:',time.time()-t0

print 'Plotting Results'

fig		=	plt.figure(1)
plt.subplot(221)
plt.plot(amplitude(XX))
plt.ylabel('XX')
plt.subplot(222)
plt.plot(amplitude(XY))
plt.ylabel('XY')
plt.subplot(223)
plt.plot(amplitude(YX))
plt.ylabel('YX')
plt.subplot(224)
plt.plot(amplitude(YY))
plt.ylabel('YY')
fig.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'(4.2).png')

fig11	=	plt.figure(11)
plt.subplot(221)
plt.plot(amplitude(mean(XX)))
plt.ylabel('XX')
plt.subplot(222)
plt.plot(amplitude(mean(XY)))
plt.ylabel('XY')
plt.subplot(223)
plt.plot(amplitude(mean(YX)))
plt.ylabel('YX')
plt.subplot(224)
plt.plot(amplitude(mean(YY)))
plt.ylabel('YY')
fig11.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'1(4.2).png')

#fig2	=	plt.figure(2)
#plt.subplot(221)
#plt.imshow(amplitude(XX),interpolation='nearest',aspect='auto')
#plt.ylabel('XX')
#plt.subplot(222)
#plt.imshow(amplitude(XY),interpolation='nearest',aspect='auto')
#plt.ylabel('XY')
#plt.subplot(223)
#plt.imshow(amplitude(YX),interpolation='nearest',aspect='auto')
#plt.ylabel('YX')
#plt.subplot(224)
#plt.imshow(amplitude(YY),interpolation='nearest',aspect='auto')
#plt.ylabel('YY')
#fig2.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_im(4.2).png')

fig21	=	plt.figure(21)
plt.subplot(221)
plt.imshow(amplitude(XX),interpolation='bilinear',aspect='auto')
plt.ylabel('XX')
plt.subplot(222)
plt.imshow(amplitude(XY),interpolation='bilinear',aspect='auto')
plt.ylabel('XY')
plt.subplot(223)
plt.imshow(amplitude(YX),interpolation='bilinear',aspect='auto')
plt.ylabel('YX')
plt.subplot(224)
plt.imshow(amplitude(YY),interpolation='bilinear',aspect='auto')
plt.ylabel('YY')
fig21.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_im2(4.2).png')

fig3	=	plt.figure(3)
plt.subplot(221)
plt.plot(amplitude(X_A1_m))
plt.ylabel('X_A1')
plt.subplot(222)
plt.plot(amplitude(X_A2_m))
plt.ylabel('X_A2')
plt.subplot(223)
plt.plot(amplitude(Y_A1_m))
plt.ylabel('Y_A1')
plt.subplot(224)
plt.plot(amplitude(Y_A2_m))
plt.ylabel('Y_A2')
fig3.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr(4.2).png')

#fig31	=	plt.figure(31)
#plt.subplot(221)
#plt.imshow(amplitude(X_A1_m),interpolation='nearest',aspect='auto')
#plt.ylabel('X1')
#plt.subplot(222)
#plt.imshow(amplitude(X_A2_m),interpolation='nearest',aspect='auto')
#plt.ylabel('X2')
#plt.subplot(223)
#plt.imshow(amplitude(Y_A1_m),interpolation='nearest',aspect='auto')
#plt.ylabel('Y1')
#plt.subplot(224)
#plt.imshow(amplitude(Y_A2_m),interpolation='nearest',aspect='auto')
#plt.ylabel('Y2')
#fig31.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr_im(4.2).png')

fig32	=	plt.figure(32)
plt.subplot(221)
plt.imshow(amplitude(X_A1_m),interpolation='bilinear',aspect='auto')
plt.ylabel('X1')
plt.subplot(222)
plt.imshow(amplitude(X_A2_m),interpolation='bilinear',aspect='auto')
plt.ylabel('X2')
plt.subplot(223)
plt.imshow(amplitude(Y_A1_m),interpolation='bilinear',aspect='auto')
plt.ylabel('Y1')
plt.subplot(224)
plt.imshow(amplitude(Y_A2_m),interpolation='bilinear',aspect='auto')
plt.ylabel('Y2')
fig32.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr_im2(4.2).png')

fig4	=	plt.figure(4)
plt.subplot(221)
plt.plot(efficiency_x1)
plt.plot([np.mean(efficiency_x1)-n_s*np.std(efficiency_x1) for i in range(len(efficiency_x1))])
plt.plot([np.mean(efficiency_x1)+n_s*np.std(efficiency_x1) for i in range(len(efficiency_x1))])
plt.ylabel('Efficiency_X1')
plt.subplot(222)
plt.plot(efficiency_x2)
plt.plot([np.mean(efficiency_x2)-n_s*np.std(efficiency_x2) for i in range(len(efficiency_x2))])
plt.plot([np.mean(efficiency_x2)+n_s*np.std(efficiency_x2) for i in range(len(efficiency_x2))])
plt.ylabel('Efficiency_X2')
plt.subplot(223)
plt.plot(efficiency_y1)
plt.plot([np.mean(efficiency_y1)-n_s*np.std(efficiency_y1) for i in range(len(efficiency_y1))])
plt.plot([np.mean(efficiency_y1)+n_s*np.std(efficiency_y1) for i in range(len(efficiency_y1))])
plt.ylabel('Efficiency_Y1')
plt.subplot(224)
plt.plot(efficiency_y2)
plt.plot([np.mean(efficiency_y2)-n_s*np.std(efficiency_y2) for i in range(len(efficiency_y2))])
plt.plot([np.mean(efficiency_y2)+n_s*np.std(efficiency_y2) for i in range(len(efficiency_y2))])
plt.ylabel('Efficiency_Y2')
fig4.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_efficiency(4.2).png')


fig5	=	plt.figure(5)
plt.subplot(221)
plt.plot(Z_X1,label='X1')
plt.plot(Z_X2,label='X2')
plt.plot(Z_XX_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: XX')
plt.subplot(222)
plt.plot(Z_X1,label='X1')
plt.plot(Z_Y2,label='Y2')
plt.plot(Z_XY_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: XY')
plt.subplot(223)
plt.plot(Z_Y1,label='Y1')
plt.plot(Z_Y2,label='Y2')
plt.plot(Z_YY_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: YY')
plt.subplot(224)
plt.plot(Z_Y1,label='Y1')
plt.plot(Z_X2,label='X2')
plt.plot(Z_YX_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: YX')
fig5.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_hilbert(4.2).png')

fig51	=	plt.figure(51)
plt.subplot(221)
plt.plot(amplitude(Z_XX_C))
plt.xlabel('Hilbert Transform: XX')
plt.subplot(222)
plt.plot(amplitude(Z_XY_C))
plt.xlabel('Hilbert Transform: XY')
plt.subplot(223)
plt.plot(amplitude(Z_YY_C))
plt.xlabel('Hilbert Transform: YY')
plt.subplot(224)
plt.plot(amplitude(Z_YX_C))
plt.xlabel('Hilbert Transform: YX')
fig51.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_hilbert2(4.2).png')

fig6	=	plt.figure(6)
plt.subplot(221)
plt.plot(delay_XX*(10**9))
plt.ylabel('Delay(ns)[XX]')
plt.xlabel('Time(s)')
plt.subplot(222)
plt.plot(delay_XY*(10**9))
plt.ylabel('Delay(ns)[XY]')
plt.xlabel('Time(s)')
plt.subplot(223)
plt.plot(delay_YY*(10**9))
plt.ylabel('Delay(ns)[YY]')
plt.xlabel('Time(s)')
plt.subplot(224)
plt.plot(delay_YX*(10**9))
plt.ylabel('Delay(ns)[YX]')
plt.xlabel('Time(s)')
fig6.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_delay(4.2).png')

fig12	=	plt.figure(12)
plt.subplot(221)
plt.plot(amplitude(XX[1]))
plt.ylabel('XX')
plt.subplot(222)
plt.plot(amplitude(XY[1]))
plt.ylabel('XY')
plt.subplot(223)
plt.plot(amplitude(YX[1]))
plt.ylabel('YX')
plt.subplot(224)
plt.plot(amplitude(YY[1]))
plt.ylabel('YY')
fig12.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'2(4.2)[1].png')


plt.show()
exit()
