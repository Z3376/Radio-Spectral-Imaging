#!/usr/bin/python

import os
from subprocess import call

file			=	open('names.txt','w+')
file.write(os.popen('ls').read())
file.close()
file			=	open('names.txt','r')
n				=	[[] for i in range(10000)]
i				=	0
p				=	file.read(1)
while(p):
	if(p!="\n"):
		n[i].append(p)
	else:
		n[i]	=	''.join(n[i])
		i		+=	1
	p			=	file.read(1)

i_t		=	str(3)
n_s		=	str(2)
mode	=	str(2)

for j in range(i):
	if(n[j][-4:]=='.mbr'):
		print 'Reading File: ',n[j]
		call(['./sync6.py',n[j],i_t,n_s,mode])

os.remove('names.txt')
