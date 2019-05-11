#!/usr/bin/python

# Cross-correlation averaging time limit calculation due to bandwidth decorrelation
# Date created: Sep16,18
# by Harsh_Grover

from subprocess import Popen,PIPE
import sys

file_name	=	sys.argv[1]
src	=	sys.argv[2]
mt	=	100

for i in range(8):
	for j in range(8):
		if(i!=j and i<j):
			p	=	Popen(['./s_alt2.py',file_name,str(i),str(j),src],stdout=PIPE,stderr=PIPE)
			output,err	=	p.communicate()
			t	=	p.returncode
			if(t<mt):
				mt=t

print(mt)	