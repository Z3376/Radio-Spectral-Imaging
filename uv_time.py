#!/usr/bin/python

# Cross-correlation averaging time limit calculation due to change in apparent baselines (UV)
# Date created: Sep14,2018
# Author: Harsh_Grover

import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.time import Time
from os.path import getsize
from os.path import isfile
from subprocess import call

c			=	299792458
pi			=	np.pi
d2r			=	pi/180
r2d			=	180/pi

Lat			=	13.6029845
#Lat			=	90
Long		=	77.4279978

nb			=	56

tile_d_u	=	5
tile_d_v	=	5


def magnitude(a):
	mag	=	np.sqrt((a[0]**2)+(a[1]**2)+(a[2]**2))
	return mag

file_name	=	sys.argv[1]
f	=	open(file_name,'r')

print('Reading File')

l	=	int(f.readline())
u=[]
v=[]
for i in range(l):
	u.append(map(float,(f.readline().split(' ')[:-1])))
for i in range(l):
	v.append(map(float,(f.readline().split(' ')[:-1])))

#plt.plot(np.take(u,[2,21],axis=1),np.take(v,[2,21],axis=1),'k')
#plt.show()
#exit()

print('Reading Complete')

di	=	5
dl	=	int(l/di)-1
du	=	np.zeros((dl,nb))
dv	=	np.zeros((dl,nb))
duv	=	np.zeros((dl,nb))
for i in range(dl):
	for j in range(nb):
		du[i][j]	=	np.abs((u[5*(i+1)][j]-u[5*i][j])/5)
		dv[i][j]	=	np.abs((v[5*(i+1)][j]-v[5*i][j])/5)
		duv[i][j]	=	np.sqrt(dv[i][j]**2+du[i][j]**2)

muv	=	duv.max(axis=0)
b=list(muv).index(max(muv))
ml=list(np.take(duv,b,axis=1)).index(max(np.take(duv,b,axis=1)))
#print(b,ml)

sti	=	ml-di
k	=	1
dist	=	0
while(dist<float(min([tile_d_v,tile_d_u]))/2):
	dist	=	magnitude([(u[sti][b]-u[sti+k][b]),(v[sti][b]-v[sti+k][b]),0])
	#print(dist)
	k	+=	1

print 'Max averaging time allowed:',k,'s'

#plt.figure()
#plt.plot(du[:])
#plt.figure()
#plt.plot(dv[:])
plt.figure()
plt.plot(duv[:])
#plt.plot(np.take(duv,[2,21],axis=1))
plt.show()

	
