#!/usr/bin/python

import sys
from os import path
from subprocess import call

file_name	=	sys.argv[1]
if(len(sys.argv)==3):
	j		=	int(sys.argv[2])
else:
	j		=	0

while(j<999):
	j_d		=	j/10
	d		=	1
	while(j_d):
		j_d	=	j_d/10
		d	+=	1
	s		=	('_'+'0'*(3-d))+str(j)
	
	if (not(path.isfile(file_name+s+".mbr"))):
		break
	print "Reading file "+file_name+s+".mbr"

	j		+=	1
	i		=	0
	size 	=	path.getsize(file_name+s+".mbr")

	print size

	mbr		=	open(file_name+s+".mbr","rb")
	hdr		=	open(file_name+s+".hdr","wb")
	while (i < size):
		hdr.write(mbr.read(32))
		i	+=	1056
		mbr.read(1024)
	mbr.close()
	hdr.close()

	print "\nFile Name:		", file_name+s
	print ".mbr File Size:		", float(size)/(1024*1024*1024), "GB"
	print ".hdr File Size:		", float(path.getsize(file_name+s+".hdr"))/1024/1024, "MB"
	print "\n"
	print "     ---Info---"
	call(["python","hdr_info2.py",file_name+s+".hdr"])


