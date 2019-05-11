#!/usr/bin/python

import sys
from os.path import getsize

file_name	=	sys.argv[1]
file		=	open(file_name,'rb')
if	(file_name[-3:]=='hdr'):
	file.read(326)
elif	(file_name[-3:]=='mbr'):
	file.read(10566)
das			=	file.read(2)
source		=	file.read(10)
a1			=	int(str(file.read(1)).encode('hex'),16)/2
a2			=	int(str(file.read(1)).encode('hex'),16)/2
a3			=	int(str(file.read(1)).encode('hex'),16)/2
a4			=	int(str(file.read(1)).encode('hex'),16)/2
lo			=	int(str(file.read(2)).encode('hex'),16)
print 'DAS_ID:		',das
print 'SOURCE:		',source
print 'Attenuator_1:	',float(a1)/2
print 'Attenuator_2:	',float(a2)/2
print 'Attenuator_3:	',float(a3)/2
print 'Attenuator_4:	',float(a4)/2
print 'LO_Freq.:	',lo
print 'Obs_Freq.(RF):	',lo-140
if	(file_name[-3:]=='hdr'):
	p_num		=	getsize(file_name)/32
	print 'No. of Packets:',p_num,"\n"
elif	(file_name[-3:]=='mbr'):
	p_num		=	getsize(file_name)/1056
	print 'No. of Packets:',p_num,"\n"
