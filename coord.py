#!/usr/bin/python
import numpy as np

pi			=	3.14159265359
d2r			=	pi/180
r2d			=	180/pi

x	=	np.zeros(8,dtype='float')
y	=	np.zeros(8,dtype='float')
id	=	4
gamma	=	np.zeros((8,8),dtype='float')
gamma[3][0]	=	167
gamma[3][1]	=	124
gamma[3][2]	=	117
gamma[3][4]	=	148
gamma[3][5]	=	139
gamma[3][6]	=	130
gamma[3][7]	=	182
x[id-1]	=	0
y[id-1]	=	0
baseline	=	np.zeros((8,8),dtype='float')
baseline[3][0]	=	105.3456913
baseline[3][1]	=	111.9512717
baseline[3][2]	=	73.62784163
baseline[3][4]	=	80.62542858
baseline[3][5]	=	57.73338988
baseline[3][6]	=	111.9270635
baseline[3][7]	=	60
for i in range(8):
	if(i!=id-1):
		#gamma[id-1][i]	=	input('Angle from north: ')
		#baseline[i]	=	float(input('Baseline: '))
		x[i]	=	baseline[id-1][i]*np.cos(90+gamma[id-1][i]*d2r)
		y[i]	=	baseline[id-1][i]*np.sin(90+gamma[id-1][i]*d2r)


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

