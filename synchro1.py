f1	=	open(file_name_1,'rb')
f2	=	open(file_name_2,'rb')

z_flag	=	1
z_ctr1	=	-1
while(z_flag):
	f1.read(28)
	p_s	=	int(f1.read(4).encode('hex'),16)
	f1.read(1024)
	z_ctr1	+=	1
	if(p_s==0):
		z_flag	=	0
	
z_flag	=	1
z_ctr2	=	-1
while(z_flag):
	f2.read(28)
	p_s	=	int(f2.read(4).encode('hex'),16)
	f2.read(1024)
	z_ctr2	+=	1
	if(p_s==0):
		z_flag	=	0

g_flag	=	0

f1.seek(-1030,2)
g_l1	=	int(f1.read(2).encode('hex'),16)

f2.seek(-1030,2)
g_l2	=	int(f2.read(2).encode('hex'),16)

g		=	0
c		=	1
while(g!=(g_l1-1)):
	f1.seek(-(1056*c)-1030,2)
	g	=	int(f1.read(2).encode('hex'),16)
	c	+=	1
p_l1	=	int(f1.read(4).encode('hex'),16)

g		=	0
c		=	1
while(g!=(g_l2-1)):
	f2.seek(-(1056*c)-1030,2)
	g	=	int(f2.read(2).encode('hex'),16)
	c	+=	1
p_l2	=	int(f2.read(4).encode('hex'),16)

f1.seek(0)
f2.seek(0)

for i in range(z_ctr1):
	f1.read(1056)
for i in range(z_ctr2):
	f2.read(1056)

f1.read(26)
f2.read(26)

g1		=	int(f1.read(2).encode('hex'),16)
g2		=	int(f2.read(2).encode('hex'),16)

p_num1		=	np.zeros(2,dtype='int')
p_num2		=	np.zeros(2,dtype='int')

gps1	=	g1
while(gps1!=(g1+1)):
	f1.read(1054)
	gps1	=	int(f1.read(2).encode('hex'),16)

gps2	=	g2
while(gps2!=(g2+1)):
	f2.read(1054)
	gps2	=	int(f2.read(2).encode('hex'),16)

if(gps1==gps2):
	p_num1[0]	=	int(f1.read(4).encode('hex'),16)
	p_num2[0]	=	int(f2.read(4).encode('hex'),16)
	if(g_l1<g_l2):
		o_t		=	g_l1-gps1
		n_p		=	p_l1-p_num1[0]+1
	else:
		o_t		=	g_l2-gps1
		n_p		=	p_l2-p_num2[0]+1

elif(gps1>gps2):
	while(gps2!=gps1):
		f2.read(1054)
		gps2	=	int(f2.read(2).encode('hex'),16)
	p_num1[0]	=	int(f1.read(4).encode('hex'),16)
	p_num2[0]	=	int(f2.read(4).encode('hex'),16)
	o_t		=	g_l2-gps1
	n_p		=	p_l2-p_num2[0]+1

else:
	g_flag	=	1
	while(gps1!=gps2):
		f1.read(1054)
		gps1	=	int(f1.read(2).encode('hex'),16)
	p_num1[0]	=	int(f1.read(4).encode('hex'),16)
	p_num2[0]	=	int(f2.read(4).encode('hex'),16)
	o_t		=	g_l1-gps1
	n_p		=	p_l1-p_num1[0]+1