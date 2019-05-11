#!/usr/bin/python

# Imaging code
# Author: Harsh_Grover

import sys
from astropy.table import Table
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from progress.bar import Bar
import os

#uvf	=	sys.argv[1]
#ccf	=	sys.argv[2]
uvf	=	'./Plots/ch00_CYG-A_20161215_101054_1_3/ch00_CYG-A_20161215_101054_1_3_5.1(U).txt'
ccf	=	'./Plots/Feb12-5/ch00_CYG-A_20161215_101054_cross_corr.txt'
uvf2=	'./Plots/ch00_CYG-A_20161215_101054_7_4/ch00_CYG-A_20161215_101054_7_4_5.1(U).txt'
ccf2=	'./Plots/Feb12-1/ch00_CYG-A_20161215_101054_000_cross_corr.txt'

z	=	Table.read(ccf,format='ascii.no_header')
z2	=	Table.read(ccf2,format='ascii.no_header')
#print len(z),len(z[1])
u	=	Table.read(uvf[:-7]+'(U).txt',format='ascii.no_header')
#print len(u),len(u[1])
v	=	Table.read(uvf[:-7]+'(V).txt',format='ascii.no_header')
#print len(v),len(v[1])
u2	=	Table.read(uvf2[:-7]+'(U).txt',format='ascii.no_header')
#print len(u),len(u[1])
v2	=	Table.read(uvf2[:-7]+'(V).txt',format='ascii.no_header')

img		=	np.zeros((240,240),dtype='complex')
img0	=	np.zeros((240,240))
img1	=	np.zeros((240,240))

das_a	= map(int,raw_input('Enter active dases: ').split())
lda		=	len(das_a)

l	=	min([len(z),len(z2),len(u[1]),len(u2[1])])

k	=	0
bar	=	Bar('Progress',max=l*lda)
for ctr in range(lda/2):
	for i in range(l):
		img[int(u[ctr][i])+120][int(v[ctr][i])+120]		=	complex(z[i][127])
		img[int(u2[ctr][i])+120][int(v2[ctr][i])+120]		=	complex(z2[i][127])
		#img0[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	np.real(complex(z[i][0]))
		#img1[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	np.real(complex(z[i][255]))
		bar.next()
for ctr in range(lda/2,lda):
	for i in range(l):
		img[int(u[ctr][i])+120][int(v[ctr][i])+120]		=	np.conj(complex(z[i][127]))
		img[int(u2[ctr][i])+120][int(v2[ctr][i])+120]		=	np.conj(complex(z2[i][127]))
		#img0[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	np.real(complex(z[i][0]))
		#img1[int(u[ctr][i])+120][int(v[ctr][i])+120]	=	np.real(complex(z[i][255]))
		bar.next()
bar.finish()

dirty	=	np.fft.fft2(img)
#dirty0	=	np.fft.fft2(img0)
#dirty1	=	np.fft.fft2(img1)
dirty_c		=	np.fft.fftshift(dirty)
#dirty_c0	=	np.fft.fftshift(dirty0)
#dirty_c1	=	np.fft.fftshift(dirty1)
dirty_i		=	np.fft.ifft2(img)
#dirty_i0	=	np.fft.ifft2(img0)
#dirty_i1	=	np.fft.ifft2(img1)
dirty_i_c	=	np.fft.fftshift(dirty_i)
#dirty_i_c0	=	np.fft.fftshift(dirty_i0)
#dirty_i_c1	=	np.fft.fftshift(dirty_i1)

cwd			=	os.getcwd()
str_das_a	=	map(str,das_a)
dirname		=	'_'.join(str_das_a)+'_image'
if not os.path.exists(dirname):
	os.makedirs(dirname)


#fig1=plt.figure()
#plt.imshow(img.T,origin='lower')
#plt.xlabel('U')
#plt.ylabel('V')
#plt.axis([-125, 125, -125, 125])
#fig1.savefig(cwd+'/'+dirname+'/imaging.png')
fig2=plt.figure()
plt.imshow(np.absolute(img).T,origin='lower')
fig2.savefig(cwd+'/'+dirname+'/imaging_digi.png')
fig3=plt.figure()
plt.imshow((dirty.real),cmap='afmhot',origin='lower')
plt.colorbar()
fig3.savefig(cwd+'/'+dirname+'/imaging_digi_realfft.png')
fig4=plt.figure()
plt.imshow((dirty.imag),cmap='afmhot',origin='lower')
plt.colorbar()
fig4.savefig(cwd+'/'+dirname+'/imaging_digi_imagfft.png')
fig5=plt.figure()
plt.imshow(np.absolute(dirty),cmap='afmhot',origin='lower')
plt.colorbar()
fig5.savefig(cwd+'/'+dirname+'/imaging_digi_fft.png')
fig6=plt.figure()
plt.imshow((dirty_i.real),cmap='afmhot',origin='lower')
plt.colorbar()
fig6.savefig(cwd+'/'+dirname+'/imaging_digi_realifft.png')
fig7=plt.figure()
plt.imshow((dirty_i.imag),cmap='afmhot',origin='lower')
plt.colorbar()
fig7.savefig(cwd+'/'+dirname+'/imaging_digi_imagifft.png')
fig8=plt.figure()
plt.imshow(np.absolute(dirty_i),cmap='afmhot',origin='lower')
plt.colorbar()
fig8.savefig(cwd+'/'+dirname+'/imaging_digi_ifft.png')
fig9=plt.figure()
plt.imshow((dirty_c.real),cmap='afmhot',origin='lower')
plt.colorbar()
fig9.savefig(cwd+'/'+dirname+'/imaging_digi_c_realfft.png')
fig10=plt.figure()
plt.imshow((dirty_c.imag),cmap='afmhot',origin='lower')
plt.colorbar()
fig10.savefig(cwd+'/'+dirname+'/imaging_digi_c_imagfft.png')
fig11=plt.figure()
plt.imshow(np.absolute(dirty_c),cmap='afmhot',origin='lower',clim=[0,0.05])
plt.colorbar()
fig11.savefig(cwd+'/'+dirname+'/imaging_digi_c_fft.pdf')
fig12=plt.figure()
plt.imshow((dirty_i_c.real),cmap='afmhot',origin='lower')
plt.colorbar()
fig12.savefig(cwd+'/'+dirname+'/imaging_digi_c_realifft.png')
fig13=plt.figure()
plt.imshow((dirty_i_c.imag),cmap='afmhot',origin='lower')
plt.colorbar()
fig13.savefig(cwd+'/'+dirname+'/imaging_digi_c_imagifft.png')
fig14=plt.figure()
plt.imshow(np.absolute(dirty_i_c),cmap='afmhot',origin='lower',clim=[0,0.0000005])
plt.colorbar()
fig14.savefig(cwd+'/'+dirname+'/imaging_digi_c_ifft.pdf')
plt.show()

