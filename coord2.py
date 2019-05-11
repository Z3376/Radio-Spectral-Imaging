#!/usr/bin/python
import numpy as np

pi			=	3.14159265359
d2r			=	pi/180
r2d			=	180/pi

x	=	np.zeros(8,dtype='float')
y	=	np.zeros(8,dtype='float')
gamma	=	np.zeros((8,8),dtype='float')
baseline	=	np.zeros((8,8),dtype='float')

x[0]	=	0
y[0]	=	0
x[1]	=	-64.5
y[1]	=	41.6
x[2]	=	-39.5
y[2]	=	68
x[3]	=	27.3
y[3]	=	99.7
x[4]	=	-19.9
y[4]	=	29.5
x[5]	=	-11.4
y[5]	=	58.1
x[6]	=	-62.6
y[6]	=	25.9
x[7]	=	26
y[7]	=	53.8

for i in range(8):
	for j in range(8):
		if(i!=j):
			xt	=	x[j]-x[i]
			yt	=	y[j]-y[i]
			baseline[i][j]	=	np.sqrt((xt**2)+(yt**2))
			gamma[i][j]	=	(np.arccos(yt/baseline[i][j]))*r2d
			print 'Baseline[',i,'][',j,']	=	',baseline[i][j]
for i in range(8):
	for j in range(8):
		if(i!=j):
			print 'Gamma[',i,'][',j,']		=	',gamma[i][j]

fbasel	=	open('Baseline.txt','w+')
fgamma	=	open('Gamma.txt','w+')
for i in range(8):
	for j in range(8):
		fbasel.write(str(baseline[i][j]))
		fgamma.write(str(gamma[i][j]))
		fbasel.write('\n')
		fgamma.write('\n')