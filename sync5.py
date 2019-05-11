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
				corr[i][j]==0

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

if(sys.argv[1]=='-d'):
	file_name_1	=	'raw_ch03_B_P_20170714_042859_000.mbr'
	file_name_2	=	'raw_ch06_B_P_20170714_042859_000.mbr'
	i_t			=	4
else:
	file_name_1	=	raw_input('File_1:')
	file_name_2	=	raw_input('File_2:')
	i_t			=	input('Integration Time(in s.): ')

########-------------------------------########
#				Synchronisation
########-------------------------------########

f1	=	open(file_name_1,'rb')
f2	=	open(file_name_2,'rb')

g_flag	=	0

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
f1.read(26)
f2.read(26)

g1		=	int(str(f1.read(2)).encode('hex'),16)
g2		=	int(str(f2.read(2)).encode('hex'),16)

p_num1		=	np.zeros(2,dtype='int')
p_num2		=	np.zeros(2,dtype='int')

gps1	=	0
while(gps1!=(g1+1)):
	f1.read(1054)
	gps1	=	int(str(f1.read(2)).encode('hex'),16)

gps2	=	0
while(gps2!=(g2+1)):
	f2.read(1054)
	gps2	=	int(str(f2.read(2)).encode('hex'),16)

if(gps1==gps2):
	p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
	p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
	o_t		=	g_l1-gps1
	n_p		=	p_l1-p_num1+1

elif(gps1>gps2):
	while(gps2!=gps1):
		f2.read(1054)
		gps2	=	int(str(f2.read(2)).encode('hex'),16)
	p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
	p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
	o_t		=	g_l2-gps1
	n_p		=	p_l2-p_num2[0]+1

else:
	g_flag	=	1
	while(gps1!=gps2):
		f1.read(1054)
		gps1	=	int(str(f1.read(2)).encode('hex'),16)
	p_num1[0]	=	int(str(f1.read(4)).encode('hex'),16)
	p_num2[0]	=	int(str(f2.read(4)).encode('hex'),16)
	o_t		=	g_l1-gps1
	n_p		=	p_l1-p_num1[0]+1

########----------------------------------########
#			Initialising Corr Variables
########----------------------------------########

off		=	p_num1[0]-p_num2[0]

a_count	=	int(n_p/o_t)

o_t			=	4		##############!!!!!!!!!!!!################

o_f			=	16								#oversampling_factor
s_f			=	33000000						#sampling_frequency
b_w			=	s_f/2							#bandwidth
padding		=	np.zeros(256+(512*(o_f-1)))

i_n			=	int(o_t/i_t)
a_count		=	int(i_t*a_count)

e_c			=	1000
a_t			=	int(a_count/e_c)

tbwp		=	e_c							#time-bandwidth-product

q1		=	0
q2		=	0

buff_1_x	=	np.zeros(1024,dtype='complex')
buff_2_x	=	np.zeros(1024,dtype='complex')
buff_1_y	=	np.zeros(1024,dtype='complex')
buff_2_y	=	np.zeros(1024,dtype='complex')
x_fftw_1	=	pyfftw.empty_aligned(512,dtype='complex')
y_fftw_1	=	pyfftw.empty_aligned(512,dtype='complex')
x_fftw_2	=	pyfftw.empty_aligned(512,dtype='complex')
y_fftw_2	=	pyfftw.empty_aligned(512,dtype='complex')
x_fftw_1	=	pyfftw.FFTW(x_fftw_1,x_fftw_1)
y_fftw_1	=	pyfftw.FFTW(y_fftw_1,y_fftw_1)
x_fftw_2	=	pyfftw.FFTW(x_fftw_2,x_fftw_2)
y_fftw_2	=	pyfftw.FFTW(y_fftw_2,y_fftw_2)
XX			=	np.zeros((i_n,256),dtype='complex')
YY			=	np.zeros((i_n,256),dtype='complex')
XY			=	np.zeros((i_n,256),dtype='complex')
YX			=	np.zeros((i_n,256),dtype='complex')
X_A1_m		=	np.zeros((i_n,256),dtype='complex')
Y_A1_m		=	np.zeros((i_n,256),dtype='complex')
X_A2_m		=	np.zeros((i_n,256),dtype='complex')
Y_A2_m		=	np.zeros((i_n,256),dtype='complex')
X_A1		=	np.zeros((i_n,a_t,256),dtype='complex')
Y_A1		=	np.zeros((i_n,a_t,256),dtype='complex')
X_A2		=	np.zeros((i_n,a_t,256),dtype='complex')
Y_A2		=	np.zeros((i_n,a_t,256),dtype='complex')
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


p_l_t	=	[i for i in range(512)]
p_l_t[0]=	512

########-----------------------------########
#			Calculating S_F_Delay 
########-----------------------------########

p_flag1	=	0
p_flag2	=	0

if(isfile(file_name_1[:-8]+'_packet_info.txt')):
	p_flag1	=	1
if(isfile(file_name_2[:-8]+'_packet_info.txt')):
	p_flag2	=	1

if(p_flag1==0):
	call(['./header2.py',file_name_1[:-8]])
	call(['./packet.py',file_name_1[:-8]])

if(p_flag2==0):
	call(['./header2.py',file_name_2[:-8]])
	call(['./packet.py',file_name_2[:-8]])

fp1	=	open(file_name_1[:-8]+'_packet_info.txt','r')
fp2	=	open(file_name_2[:-8]+'_packet_info.txt','r')

fp1.seek(0,2)
fp2.seek(0,2)
ctr	=	0

ctr1	=	0
ctr2	=	0
while(ctr1!=4):
	ctr2	+=	1
	fp1.seek(-ctr2,2)
	if(fp1.read(1)=='E'):
		ctr1	+=	1

ctr2	+=	3
fp1.seek(-ctr2,2)
ctr3	=	1
while(fp1.read(1)!='	'):
	ctr2	+=	1
	ctr3	+=	1
	fp1.seek(-ctr2,2)

p_blip1	=	float(fp1.read(ctr3))

ctr1	=	0
ctr2	=	0
while(ctr1!=4):
	ctr2	+=	1
	fp2.seek(-ctr2,2)
	if(fp2.read(1)=='E'):
		ctr1	+=	1

ctr2	+=	3
fp2.seek(-ctr2,2)
ctr3	=	1
while(fp2.read(1)!='	'):
	ctr2	+=	1
	ctr3	+=	1
	fp2.seek(-ctr2,2)

p_blip2	=	float(fp2.read(ctr3))

s_f_d1	=	p_blip1-int(p_blip1)
s_f_d2	=	p_blip2-int(p_blip2)

s_f_d	=	s_f_d1-s_f_d2
s_f_del	=	1056*s_f_d
p_del	=	s_f_del-int(s_f_del)
s_f_del	=	int(s_f_del)

s_f_1	=	0
s_f_2	=	0

ctr		=	1
while(fp1.read(1)!='	'):
	ctr	+=	1
	fp1.seek(-ctr,2)

s_f_1	=	float(fp1.read(ctr))

ctr	=	1
while(fp2.read(1)!='	'):
	ctr	+=	1
	fp2.seek(-ctr,2)

s_f_2	=	float(fp2.read(ctr))

fp1.close()
fp2.close()

s_f_dif	=	s_f_1-s_f_2

s_f_dev	=	np.zeros(o_t)

p_s_f_dev	=	np.zeros(o_t)

for i in range(o_t):
	for j in range(i):
		s_f_dev[o_t-j-1]	+=	s_f_dif

for i in range(o_t):
	p_s_f_dev[i]	=	s_f_dev[i]-int(s_f_dev[i])
	s_f_dev[i]		=	int(s_f_dev[i])

s_f_dev_i	=	np.zeros(i_n)
p_s_f_dev_i	=	np.zeros(i_n)

for i in range(i_n):
	for j in range(i_t):
		s_f_dev_i[i]	+=	s_f_dev[(i*i_t)+j]
	s_f_dev_i[i]	=	s_f_dev_i[i]/i_t

for i in range(i_n):
	for j in range(i_t):
		p_s_f_dev_i[i]	+=	p_s_f_dev[(i*i_t)+j]
	p_s_f_dev_i[i]	=	p_s_f_dev_i[i]/i_t

#print p_del, s_f_del
#print s_f_dev,s_f_dev_i

########-----------------------------------########
#			Calculating Geometric Delay
########-----------------------------------########

das_id_1	=	3
das_id_2	=	6
d_file_name	=	'0'
if(g_flag==0):
	if(not(isfile(file_name_1[:-8]+'_g_delay.txt'))):
		call(['./g_delay3.py',file_name_1[:-3]+'hdr',das_id_1,das_id_2,'cas'])
	d_file_name	=	file_name_1[:-8]+'_g_delay.txt'
else:
	if(not(isfile(file_name_2[:-8]+'_g_delay.txt'))):
		call(['./g_delay3.py',file_name_2[:-3]+'hdr',das_id_1,das_id_2,'cas'])
	d_file_name	=	file_name_2[:-8]+'_g_delay.txt'

d_f	=	open(d_file_name,'r')
g_del	=	np.zeros(o_t)
ctr1	=	0
ctr2	=	0
while(ctr2<o_t+1):
	b	=	d_f.read(1)
	if(b=='\n'):
		if(ctr2>0):
			d_f.seek(-ctr1-1,1)
			g_del[ctr2-1]	=	float(d_f.read(ctr1))
			d_f.read(1)
		ctr1	=	0
		ctr2	+=	1
	else:
		ctr1	+=	1

g_del_i	=	np.zeros(i_n)

for i in range(i_n):
	for j in range(i_t):
		g_del_i[i]	+=	g_del[(i*i_t)+j]
	g_del_i[i]	=	g_del_i[i]/i_t

g_del_p	=	np.zeros(i_n)
p_g_del	=	np.zeros(i_n)
t_p_del	=	np.zeros(i_n)

for i in range(i_n):
	g_del_p[i]	=	g_del_i[i]/(10**9)/(1/float(s_f))
	p_g_del[i]	=	g_del_p[i]-int(g_del_p[i])
	g_del_p[i]	=	int(g_del_p[i])

#print g_del_p,p_g_del
#print g_del_i

########-------------------------------------########
#			Calculating Cross Correlation
########-------------------------------------########

print 'Starting Corr'

flag	=	0
if(s_f_del<0):
	j1_s	=	0
	j2_s	=	np.abs(s_f_del)
	j1_sp	=	0
	j2_sp	=	np.abs(p_del)
else:
	j1_s	=	np.abs(s_f_del)
	j2_s	=	0
	j1_sp	=	np.abs(p_del)
	j2_sp	=	0
j1_o	=	[j1_s for i in range(i_n)]
j2_o	=	[j2_s for i in range(i_n)]
j1_p	=	[j1_sp	for i in range(i_n)]
j2_p	=	[j2_sp	for i in range(i_n)]
for i in range(i_n):
	if(g_del_p[i]<0):
		j2_o[i]	+=	np.abs(g_del_p[i])
		j2_p[i]	+=	np.abs(p_g_del[i])
	else:
		j1_o[i]	+=	np.abs(g_del_p[i])
		j1_p[i]	+=	np.abs(p_g_del[i])
	if(s_f_dev_i[i]>0):
		j2_o[i]	+=	np.abs(s_f_dev_i[i])
		j2_p[i]	+=	np.abs(p_s_f_dev_i[i])
	else:
		j1_o[i]	+=	np.abs(s_f_dev_i[i])
		j1_p[i]	+=	np.abs(p_s_f_dev_i[i])

for i in range(i_n):
	j1_o[i]	+=	int(j1_p[i])
	j1_p[i]	=	j1_p[i]-int(j1_p[i])
	j2_o[i]	+=	int(j2_p[i])
	j2_p[i]	=	j2_p[i]-int(j2_p[i])

f1.seek(int(j1_o[0]),1)
f2.seek(int(j2_o[0]),1)
j1	=	int(j1_o[0])
j2	=	int(j2_o[0])
print j1_o,j2_o
for i in range(i_n):
	if(i>0):
		if(int(j1_o[i]-j1_o[i-1])>0 and int(j2_o[i]-j2_o[i-1])>0):
			f1.seek(int(j1_o[i]-j1_o[i-1]),1)
			f2.seek(int(j2_o[i]-j2_o[i-1]),1)
			j1	+=	int(j1_o[i]-j1_o[i-1])
			j2	+=	int(j2_o[i]-j2_o[i-1])
		elif(int(j1_o[i]-j1_o[i-1])<0 and int(j2_o[i]-j2_o[i-1])>0):
			f2.seek(int(j2_o[i]-j2_o[i-1]),1)
			j2	+=	int(j2_o[i]-j2_o[i-1])
			if(np.abs(int(j1_o[i]-j1_o[i-1]))<j1):
				f1.seek(int(j1_o[i]-j1_o[i-1]),1)
				j1	+=	int(j1_o[i]-j1_o[i-1])
			else:
				f1.seek(-j1-32,1)
				f1.seek(int(j1_o[i]-j1_o[i-1])+j1,1)
				j1	+=	512+2*int(j1_o[i]-j1_o[i-1])
		elif(int(j1_o[i]-j1_o[i-1])>0 and int(j2_o[i]-j2_o[i-1])<0):
			f1.seek(int(j1_o[i]-j1_o[i-1]),1)
			j1	+=	int(j1_o[i]-j1_o[i-1])
			if(np.abs(int(j2_o[i]-j2_o[i-1]))<j2):
				f2.seek(int(j2_o[i]-j2_o[i-1]),1)
				j2	+=	int(j2_o[i]-j2_o[i-1])
			else:
				f2.seek(-j2-32,1)
				f2.seek(int(j2_o[i]-j2_o[i-1])+j2,1)
				j2	+=	512+int(j2_o[i]-j2_o[i-1])
		else:
			if(np.abs(int(j1_o[i]-j1_o[i-1]))<j1):
				f1.seek(int(j1_o[i]-j1_o[i-1]),1)
				j1	+=	int(j1_o[i]-j1_o[i-1])
			else:
				f1.seek(-j1-32,1)
				f1.seek(int(j1_o[i]-j1_o[i-1])+j1,1)
				j1	+=	512+int(j1_o[i]-j1_o[i-1])
			if(np.abs(int(j2_o[i]-j2_o[i-1]))<j2):
				f2.seek(int(j2_o[i]-j2_o[i-1]),1)
				j2	+=	int(j2_o[i]-j2_o[i-1])
			else:
				f2.seek(-j2-32,1)
				f2.seek(int(j2_o[i]-j2_o[i-1])+j2,1)
				j2	+=	512+int(j2_o[i]-j2_o[i-1])
	for k in range(a_t):
		while(j1<(i+1)*(k+1)*512*e_c):
			buff_1_x[np.abs(int(j1-j1_o[i]))%512]	=	l_t[int(str(f1.read(1)).encode('hex'),16)]
			buff_2_x[np.abs(int(j2-j2_o[i]))%512]	=	l_t[int(str(f2.read(1)).encode('hex'),16)]
			buff_1_y[np.abs(int(j1-j1_o[i]))%512]	=	l_t[int(str(f1.read(1)).encode('hex'),16)]
			buff_2_y[np.abs(int(j2-j2_o[i]))%512]	=	l_t[int(str(f2.read(1)).encode('hex'),16)]
			j1				+=	1
			j2				+=	1
			if(j2%512==0):
				q2			+=	1
				f2.read(28)
				p_num2[q2%2]	=	int(str(f2.read(4)).encode('hex'),16)
				flag		+=	1
				if((flag%3)==2):
					flag		=	0
					X_1			=	x_fftw_1(buff_1_x)
					Y_1			=	y_fftw_1(buff_1_y)
					X_2			=	x_fftw_2(buff_2_x)
					Y_2			=	y_fftw_2(buff_2_y)
					X_A1[i][k]	+=	X_1[:256]*np.conj(X_1[:256])
					Y_A1[i][k]	+=	Y_1[:256]*np.conj(Y_1[:256])
					X_A2[i][k]	+=	X_2[:256]*np.conj(X_2[:256])
					Y_A2[i][k]	+=	Y_2[:256]*np.conj(Y_2[:256])
					XX[i]		+=	X_1[:256]*np.conj(X_2[:256])
					YY[i]		+=	Y_1[:256]*np.conj(Y_2[:256])
					XY[i]		+=	X_1[:256]*np.conj(Y_2[:256])
					YX[i]		+=	Y_1[:256]*np.conj(X_2[:256])
				if(not(trans(p_num2))):
					print 'Slip at: ',p_num1[q1%2],p_num2[q2%2]
					clear(buff_1)
					clear(buff_2)
					#print 'Done clearing buffers'
					flag		=	0
					while(p_num1[q1%2]-p_num2[q2%2]!=off):
						if(p_num1[q1%2]-p_num2[q2%2]>off):
							j1	+=	512
							j2	+=	512
							p	=	j2%512
							f2.read(1052-2*p)
							p_num2[q2%2]	=	int(str(f2.read(4)).encode('hex'),16)
							f2.read(p)
						else:
							j1	+=	512
							j2	+=	512
							p	=	p_l_t[j1%512]
							f1.read(1052-2*p)
							p_num1[q1%2]	=	int(str(f1.read(4)).encode('hex'),16)
							f1.read(p)
			if(j1%512==0):
				q1			+=	1				
				f1.read(28)
				p_num1[q1%2]		=	int(str(f1.read(4)).encode('hex'),16)
				flag		+=	1
				if((flag%3)==2):
					flag		=	0
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
				if(not(trans(p_num1))):
					print 'Slip at: ',p_num1[q1%2],p_num2[q2%2]
					clear(buff_1)
					clear(buff_2)
					#print 'Done clearing buffers'
					flag		=	0
					while(p_num1[q1%2]-p_num2[q2%2]!=off):
						if(p_num1[q1%2]-p_num2[q2%2]>off):
							j1	+=	512
							j2	+=	512
							p	=	j2%512
							f2.read(1052-2*p)
							p_num2[q2%2]	=	int(str(f2.read(4)).encode('hex'),16)
							f2.read(p)
						else:
							j1	+=	512
							j2	+=	512
							p	=	j1%512
							f1.read(1052-2*p)
							p_num1[q1%2]	=	int(str(f1.read(4)).encode('hex'),16)
							f1.read(p)			
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
	XX[i]		=	XX[i]/float(a_count)/np.sqrt(mean(X_A1[i])*mean(X_A2[i]))
	XY[i]		=	XY[i]/float(a_count)/np.sqrt(mean(X_A1[i])*mean(Y_A2[i]))
	YY[i]		=	YY[i]/float(a_count)/np.sqrt(mean(Y_A1[i])*mean(Y_A2[i]))
	YX[i]		=	YX[i]/float(a_count)/np.sqrt(mean(Y_A1[i])*mean(X_A2[i]))
	print i

f1.close()
f2.close()

########------------------------------########
#			Calculating Efficiency
########------------------------------########

print 'Calculating Efficiency'

SNR_exp			=	float(np.sqrt(tbwp))

efficiency_x1	=	efficiency(X_A1,SNR_exp)
efficiency_x2	=	efficiency(X_A2,SNR_exp)
efficiency_y1	=	efficiency(Y_A1,SNR_exp)
efficiency_y2	=	efficiency(Y_A2,SNR_exp)

########-----------------------------########
#			Clearing Bad Channels
########-----------------------------########

print 'Clearing Bad Channels'

R_X1		=	reject_matrix(efficiency_x1,3)
R_X2		=	reject_matrix(efficiency_x2,3)
R_Y1		=	reject_matrix(efficiency_y1,3)
R_Y2		=	reject_matrix(efficiency_y2,3)

reject(R_X1,R_X2,XX)
reject(R_X1,R_Y2,XY)
reject(R_Y1,R_Y2,YY)
reject(R_Y1,R_X2,YX)

########-----------------------------########
#			  Hilbert Transform
########-----------------------------########

print 'Performing Hilbert Transform'

hilbert_delay(i_n,Z_XX,XX,padding,s_f,o_f,delay_XX)
hilbert_delay(i_n,Z_XY,XY,padding,s_f,o_f,delay_XY)
hilbert_delay(i_n,Z_YY,YY,padding,s_f,o_f,delay_YY)
hilbert_delay(i_n,Z_YX,YX,padding,s_f,o_f,delay_YX)

for i in range(i_n):
	residual_delay_XX	=	delay_XX[i]-j1_o[i]-j2_o[i]
	residual_delay_XY	=	delay_XY[i]-j1_o[i]-j2_o[i]
	residual_delay_YY	=	delay_YY[i]-j1_o[i]-j2_o[i]
	residual_delay_YX	=	delay_YX[i]-j1_o[i]-j2_o[i]

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

########-----------------------########
#			Saving Results
########-----------------------########

print 'Saving Results'

XX_t	=	Table(XX)
XX_t.write(file_name_1[:-4]+'=X='+file_name_2[:-4]+'.txt',format='ascii',overwrite='true')
fw		=	open(file_name_1[:-4]+'=X='+file_name_2[:-4]+'.dat',mode='a')
fw.seek(0,2)
XY_t	=	Table(XY)
XY_t.write(fw,format='ascii.no_header')
YY_t	=	Table(YY)
YY_t.write(fw,format='ascii.no_header')
YX_t	=	Table(YX)
YX_t.write(fw,format='ascii.no_header')
fw.close()

print 'Total time elapsed:',time.time()-t0

########-------------------------########
#			Plotting Results
########-------------------------########

print 'Plotting Results'

fig		=	plt.figure(1)
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

fig11	=	plt.figure(11)
plt.subplot(221)
plt.plot(mean(XX))
plt.ylabel('XX')
plt.subplot(222)
plt.plot(mean(XY))
plt.ylabel('XY')
plt.subplot(223)
plt.plot(mean(YX))
plt.ylabel('YX')
plt.subplot(224)
plt.plot(mean(YY))
plt.ylabel('YY')
fig11.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'1.png')

fig2	=	plt.figure(2)
plt.subplot(221)
plt.imshow(amplitude(XX),interpolation='nearest',aspect='auto')
plt.ylabel('XX')
plt.subplot(222)
plt.imshow(amplitude(XY),interpolation='nearest',aspect='auto')
plt.ylabel('XY')
plt.subplot(223)
plt.imshow(amplitude(YX),interpolation='nearest',aspect='auto')
plt.ylabel('YX')
plt.subplot(224)
plt.imshow(amplitude(YY),interpolation='nearest',aspect='auto')
plt.ylabel('YY')
fig2.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_im.png')

fig21	=	plt.figure(2)
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
fig21.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_im2.png')

fig3	=	plt.figure(3)
plt.subplot(221)
plt.plot(X_A1_m)
plt.ylabel('X_A1')
plt.subplot(222)
plt.plot(X_A2_m)
plt.ylabel('X_A2')
plt.subplot(223)
plt.plot(Y_A1_m)
plt.ylabel('Y_A1')
plt.subplot(224)
plt.plot(Y_A2_m)
plt.ylabel('Y_A2')
fig3.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr.png')

fig31	=	plt.figure(31)
plt.subplot(221)
plt.imshow(amplitude(XX),interpolation='nearest',aspect='auto')
plt.ylabel('XX')
plt.subplot(222)
plt.imshow(amplitude(XY),interpolation='nearest',aspect='auto')
plt.ylabel('XY')
plt.subplot(223)
plt.imshow(amplitude(YX),interpolation='nearest',aspect='auto')
plt.ylabel('YX')
plt.subplot(224)
plt.imshow(amplitude(YY),interpolation='nearest',aspect='auto')
plt.ylabel('YY')
fig31.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr_im.png')

fig32	=	plt.figure(32)
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
fig32.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_auto_corr_im2.png')

fig4	=	plt.figure(4)
plt.subplot(221)
plt.plot(efficiency_x1)
plt.plot([np.mean(efficiency_x1)-3*np.std(efficiency_x1) for i in range(len(efficiency_x1))])
plt.plot([np.mean(efficiency_x1)+3*np.std(efficiency_x1) for i in range(len(efficiency_x1))])
plt.ylabel('Efficiency_X1')
plt.subplot(222)
plt.plot(efficiency_x2)
plt.plot([np.mean(efficiency_x2)-3*np.std(efficiency_x2) for i in range(len(efficiency_x2))])
plt.plot([np.mean(efficiency_x2)+3*np.std(efficiency_x2) for i in range(len(efficiency_x2))])
plt.ylabel('Efficiency_X2')
plt.subplot(223)
plt.plot(efficiency_y1)
plt.plot([np.mean(efficiency_y1)-3*np.std(efficiency_y1) for i in range(len(efficiency_y1))])
plt.plot([np.mean(efficiency_y1)+3*np.std(efficiency_y1) for i in range(len(efficiency_y1))])
plt.ylabel('Efficiency_Y1')
plt.subplot(224)
plt.plot(efficiency_y2)
plt.plot([np.mean(efficiency_y2)-3*np.std(efficiency_y2) for i in range(len(efficiency_y2))])
plt.plot([np.mean(efficiency_y2)+3*np.std(efficiency_y2) for i in range(len(efficiency_y2))])
plt.ylabel('Efficiency_Y2')
fig4.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_efficiency.png')


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
fig5.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_hilbert.png')

fig51	=	plt.figure(51)
plt.subplot(221)
plt.plot(Z_XX_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: XX')
plt.subplot(222)
plt.plot(Z_XY_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: XY')
plt.subplot(223)
plt.plot(Z_YY_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: YY')
plt.subplot(224)
plt.plot(Z_YX_C,label='Cross Corr')
plt.legend(bbox_to_anchor=(0.,1.02,1.,.102),loc=3,ncol=2,mode="expand",borderaxespad=0.)
plt.xlabel('Hilbert Transform: YX')
fig51.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_hilbert2.png')


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
fig6.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_delay.png')

fig61	=	plt.figure(61)
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
fig61.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'_residual_delay.png')

plt.show()
