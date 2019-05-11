#!/usr/bin/python

import sys
from os.path import getsize
from os.path import isfile
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

def int2str(a,b):
	c							=	a/10
	d							=	1
	while(c):
		c						/=	10
		d						+=	1
	st							=	'0'*(b-d)+str(a)
	return	st

def double_arr_to_list(a):
	b							=	[]
	for i in range(len(a)):
		for j in range(len(a[i])):
			b.append(a[i][j])
	return	b

def slope(a,b):
	s							=	[0 for i in range(len(a)-1)]
	for i in range(len(a)-1):
		s[i]					=	float(b[i+1]-b[i])/(a[i+1]-a[i])
	return	s

def linear_fit(a,b):
	z 							=	np.polyfit(a,b,1)
	r 							= 	np.poly1d(z)
	a_new						=	np.linspace(a[0],a[-1],len(a))
	b_new						= 	r(a_new)
	return	a_new,b_new

if(len(sys.argv)==2):
	file_name					=	sys.argv[1]
elif(len(sys.argv)==1):
	file_name					=	raw_input(".hdr_File: ")
else:
	print("Too many arguments. Enter only .hdr_File name.")
	exit()

jn								=	int(file_name[-7:-4])
file_name						=	file_name[:-8]
G_BLIP							=	[]
P_BLIP							=	[]
F_G_BLIP						=	[]
F_P_BLIP						=	[]
stop							=	[]
SLOPE							=	[]
stability						=	[]
stability_fpga					=	[]
conv							=	[]
conv_fpga						=	[]
julian_day						=	[]
j								=	0
st								=	'_'+int2str(jn,3)

txt_file						=	open(file_name+"_packet_info.txt","w+")
txt_file.write("\t\t---Sampling_Frequency_Summary---")
txt_file.close()
txt_file						=	open(file_name+"_packet_info.txt","a")

#stop.append(0)

while(jn<999):
		
	print "Reading file "+file_name+st+".hdr"

	size						=	getsize(file_name+st+".hdr")
	i							=	26
	q							=	0
	fpga						=	np.empty(2)
	gps							=	np.empty(2)
	p_num 						=	np.empty(2)
	g_blip						=	[]
	p_blip						=	[]
	f_g_blip 					=	[]
	f_p_blip 					=	[]

	file						=	open(file_name+st+".hdr","rb")
	while(i < size):
		file.read(24)
		fpga[q%2]				=	int(str(file.read(2)).encode('hex'),16)
		gps[q%2]				=	int(str(file.read(2)).encode('hex'),16)
		p_num[q%2]				=	int(str(file.read(4)).encode('hex'),16)
		if (gps[q%2]-gps[(q+1)%2]==1 and p_num[q%2]-p_num[(q+1)%2]==1):
			g_blip.append(gps[q%2])
			p_blip.append(p_num[q%2])
			if (fpga[q%2]!=fpga[(q+1)%2]):
				f_g_blip.append(gps[q%2])
				f_p_blip.append(p_num[q%2])
		i						+=	32
		q						+=	1
	file.close()

	s							=	slope(g_blip,p_blip)

	g_blip_new,p_blip_new 		=	linear_fit(g_blip, p_blip)
	f_g_blip_new,f_p_blip_new	=	linear_fit(f_g_blip, f_p_blip)
	

	if(j==0):
		start					=	g_blip[0]-1
		stop.append(g_blip[0])

	time_s						=	int(file_name[-2:])+stop[j]-start
	if(time_s<60):
		st_time_s				=	int2str(time_s,2)
		s_min					=	0
	else:
		while(time_s>=60):
			s_min				=	int(time_s/60)
			time_s				-=	s_min*60
			st_time_s			=	int2str(time_s,2)
	
	time_min					=	int(file_name[-4:-2])+s_min
	if(time_min<60):
		st_time_min				=	int2str(time_min,2)
		min_hr					=	0
	else:
		while(time_min>=60):
			min_hr				=	int(time_min/60)
			time_min			-=	min_hr*60
			st_time_min			=	int2str(time_min,2)
	
	time_hr						=	int(file_name[-6:-4])+min_hr
	if(time_hr<24):
		st_time_hr				=	int2str(time_hr,2)
		hr_day					=	0
	else:
		while(time_hr>=24):
			hr_day				=	int(time_hr/24)
			time_hr				-=	hr_day*24
			st_time_hr			=	int2str(time_hr,2)

	time_yr						=	file_name[-15:-11]
	time_month					=	file_name[-11:-9]
	time_day					=	int(file_name[-9:-7])+hr_day
	st_time_day					=	int2str(time_day,2)


	t1							=	time_yr+"-"+time_month+"-"+st_time_day+"T"+st_time_hr+":"+st_time_min+":"+st_time_s
#print t1
	t2							=	Time(t1,format="isot")

	julian_day.append(t2.jd)

	stop.append(gps[-1])

	txt_file.write("\nFile Name:\t\t\t\t\t"+file_name+st+".hdr")
	txt_file.write("\nJulian Day:\t\t\t\t\t"+str(julian_day[j]))
	txt_file.write("\nGPS counter for first transition:\t\t"+str(g_blip[0]))
	txt_file.write("\nP_Num for first transition:\t\t\t"+str(p_blip[0]))
	txt_file.write("\nEstimate of P_Num for first transition:\t\t"+str(p_blip_new[0]))
	txt_file.write("\nEstimate of P_Num for first FPGA Blip:\t\t"+str(f_p_blip_new[0]))
	txt_file.write("\nEstimated Packets/Sec:\t\t\t\t"+str(np.mean(slope(g_blip_new,p_blip_new))))
	txt_file.write("\nEstimated Packets/Sec(FPGA Blip):\t\t"+str(np.mean(slope(f_g_blip_new,f_p_blip_new))))
	txt_file.write("\nEstimated Sampling Freq.:\t\t\t"+str(np.mean(slope(g_blip_new,p_blip_new))*512))
	txt_file.write("\nEstimated Sampling Freq.(FPGA Blip):\t\t"+str(np.mean(slope(f_g_blip_new,f_p_blip_new))*512))
	txt_file.write("\n")
	
	G_BLIP.append(g_blip)
	P_BLIP.append(p_blip)
	F_G_BLIP.append(f_g_blip)
	F_P_BLIP.append(f_p_blip)
	stability.append(np.mean(slope(g_blip_new,p_blip_new))*512)
	stability_fpga.append(np.mean(slope(f_g_blip_new,f_p_blip_new))*512)
	SLOPE.append(s)

	G_BLIP_c					=	double_arr_to_list(G_BLIP)
	P_BLIP_c					=	double_arr_to_list(P_BLIP)
	F_G_BLIP_c 					=	double_arr_to_list(F_G_BLIP)
	F_P_BLIP_c 					=	double_arr_to_list(F_P_BLIP)
	SLOPE_c						=	double_arr_to_list(SLOPE)

	G_BLIP_new,P_BLIP_new 		=	linear_fit(G_BLIP_c, P_BLIP_c)
	F_G_BLIP_new,F_P_BLIP_new	=	linear_fit(F_G_BLIP_c, F_P_BLIP_c)
	
	conv.append(np.mean(slope(G_BLIP_new,P_BLIP_new))*512)
	conv_fpga.append(np.mean(slope(F_G_BLIP_new,F_P_BLIP_new))*512)

	j							+=	1
	jn							+=	1
	st							=	'_'+int2str(j,3)
		
	if (not(isfile(file_name+st+".hdr"))):
		txt_file.write("\n")
		txt_file.write("\t\t---For_Whole_Observation---")
		txt_file.write("\nP_Num for first transition:\t\t\t"+str(P_BLIP_c[0]))
		txt_file.write("\nP_Num for first FPGA Blip:\t\t\t"+str(F_P_BLIP_c[0]))
		txt_file.write("\nEstimate of P_Num for first transition:\t\t"+str(P_BLIP_new[0]))
		txt_file.write("\nEstimate of P_Num for first FPGA Blip:\t\t"+str(F_P_BLIP_new[0]))
		txt_file.write("\nEstimated Packets/Sec:\t\t\t\t"+str(np.mean(slope(G_BLIP_new,P_BLIP_new))))
		txt_file.write("\nEstimated Packets/Sec(FPGA Blip):\t\t"+str(np.mean(slope(F_G_BLIP_new,F_P_BLIP_new))))
		txt_file.write("\nEstimated Sampling Freq.:\t\t\t"+str(np.mean(slope(G_BLIP_new,P_BLIP_new))*512))
		txt_file.write("\nEstimated Sampling Freq.(FPGA Blip):\t\t"+str(np.mean(slope(F_G_BLIP_new,F_P_BLIP_new))*512))

		print("\t\t---For_Whole_Observation---")
		print("\nP_Num for first transition:\t\t\t"+str(P_BLIP_c[0]))
		print("\nP_Num for first FPGA Blip:\t\t\t"+str(F_P_BLIP_c[0]))
		print("\nEstimate of P_Num for first transition:\t\t"+str(P_BLIP_new[0]))
		print("\nEstimate of P_Num for first FPGA Blip:\t\t"+str(F_P_BLIP_new[0]))
		print("\nEstimated Packets/Sec:\t\t\t\t"+str(np.mean(slope(G_BLIP_new,P_BLIP_new))))
		print("\nEstimated Packets/Sec(FPGA Blip):\t\t"+str(np.mean(slope(F_G_BLIP_new,F_P_BLIP_new))))
		print("\nEstimated Sampling Freq.:\t\t\t"+str(np.mean(slope(G_BLIP_new,P_BLIP_new))*512))
		print("\nEstimated Sampling Freq.(FPGA Blip):\t\t"+str(np.mean(slope(F_G_BLIP_new,F_P_BLIP_new))*512))

		plt.figure()
		plt.plot(G_BLIP_c,P_BLIP_c,"ro")
		plt.plot(G_BLIP_new,P_BLIP_new)
		plt.xlabel("GPS Counter")
		plt.ylabel("Packet Counter")
		plt.figure()
		plt.plot(F_G_BLIP_c,F_P_BLIP_c,"ro")
		plt.plot(F_G_BLIP_new,F_P_BLIP_new)
		plt.xlabel("F_GPS Counter")
		plt.ylabel("F_Packet Counter")
		plt.figure()
		plt.plot(SLOPE_c,"o")
		plt.xlabel("GPS Counter")
		plt.ylabel("Packets/Sec")
		plt.figure()
		plt.plot(julian_day,stability)
		plt.plot(julian_day,stability_fpga)
		plt.xlabel("Julian Day")
		plt.ylabel("Sampling Frequency")
		plt.figure()
		plt.plot(conv)
		plt.plot(conv_fpga)
		plt.ylabel("Sampling Frequency")
		plt.xlabel("No. of files considered")

		break

txt_file.close()

plt.show()
exit()
