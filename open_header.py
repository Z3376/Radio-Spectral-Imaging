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

for j in range(i):
	if(n[j][-4:]=='.mbr'):
		if(n[j][-8:-4]=='_000'):
			call(['./header2.py',n[j][:-8]])

os.remove('names.txt')
