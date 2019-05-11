import numpy as np

#file_name	=	raw_input('File: ')
file_name	=	'sweetspots.txt'
#file	=	open(file_name)
with open(file_name,'r') as file:
	lines	=	file.readlines()

#for i in range(len(lines[2])):
#	print ord(lines[2][i]),'--',lines[2][i]
#exit()

l	=	len(lines)
a	=	[[] for i in range(2)]
i	=	2
while(i<l):
	j		=	0
	ctr		=	0
	ctr1	=	0
	ctr2	=	0
	while(j<(len(lines[i]))):
		if(ord(lines[i][j])!=32):
			#print lines[i][j]
			ctr	+=	1
		if(ord(lines[i][j])==32 and ctr!=0):
			a[ctr2].append(float(lines[i][ctr1:ctr1+ctr]))
			ctr1	+=	ctr+1
			ctr2	+=	1
			ctr		=	0
		if(ord(lines[i][j])==32 and ctr==0):
			ctr1	+=	1
		if(ctr2==2):
			i	+=	1
			j	=	0
			ctr		=	0
			ctr1	=	0
			ctr2	=	0
		j	+=	1
	i	+=	1
print a[1]