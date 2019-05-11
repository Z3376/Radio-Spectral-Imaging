#!/usr/bin/python

from os.path import getsize
import numpy as np
import matplotlib.pyplot as plt
import sys
import pyfftw

def trans(a):
    if(np.abs(a[0]-a[1])==1):
        return 1
    else:
        return 0

if(len(sys.argv)==2):
    if(sys.argv[1]=='-d'):
        file_name_1 =   'raw_ch07_B_P_20170714_042859_019.mbr'
        file_name_2 =   'raw_ch06_B_P_20170714_042859_019.mbr'
else:
    file_name_1 =   sys.argv[1]
    file_name_2 =   sys.argv[2]

f1      =   open(file_name_1,'rb')
f2      =   open(file_name_2,'rb')

gps_1   =   np.empty(2)
gps_2   =   np.empty(2)
p_num_1 =   np.empty(2)
p_num_2 =   np.empty(2)

l_t     =   [i for i in range(256)]
for j in range(128):
    l_t[128+j]  -=  256

size_1  =   getsize(file_name_1)
size_2  =   getsize(file_name_2)

if(size_1<size_2):
    n_p     =   size_1/1056
    f1.read((n_p-1)*1056)
    f1.read(26)
    g_l     =   int(str(f1.read(2)).encode('hex'),16)
    n_p     =   int(str(f1.read(4)).encode('hex'),16)
else:
    n_p     =   size_2/1056
    f2.read((n_p-1)*1056)
    f2.read(26)
    g_l     =   int(str(f2.read(2)).encode('hex'),16)
    n_p     =   int(str(f2.read(4)).encode('hex'),16)

q           =   0 
status      =   0
n_flag      =   0
a_count     =   n_p/g_l
g_blip      =   []
p_blip      =   []

g_l_i       =   [i for i in range(g_l)]

print 'Calculated n_p'

f1.close()
f2.close()
f1      =   open(file_name_1,'rb')
f2      =   open(file_name_2,'rb')

i       =   0

while(i<n_p):
    f1.read(26)
    f2.read(26)
    gps_1[q%2]      =   int(str(f1.read(2)).encode('hex'),16)
    gps_2[q%2]      =   int(str(f2.read(2)).encode('hex'),16)
    p_num_1[q%2]    =   int(str(f1.read(4)).encode('hex'),16)
    p_num_2[q%2]    =   int(str(f2.read(4)).encode('hex'),16)
    f1.read(1024)
    f2.read(1024)
    print gps_1[q%2],gps_2[q%2],p_num_1[q%2],p_num_2[q%2]
    if(trans(gps_1) and trans(gps_2)):
        status      =   1
        print 'status = 1'
        break
    elif(not(trans(gps_1)) and not(trans(gps_2))):
        continue
    else:
        if(len(p_blip)==4):
            status      =   0
            print 'status = 0'
            break
        if(trans(gps_1)):
            p_blip.append(p_num_1[q%2])
            p_blip.append(p_num_2[q%2])
            print '1'
            continue
        else:
            n_flag      =   1
            p_blip.append(p_num_1[q%2])
            p_blip.append(p_num_2[q%2])
            print '1'
            continue

f1.close()
f2.close()

f1      =   open(file_name_1,'rb')
f2      =   open(file_name_2,'rb')

if(status==1):
    d   =   0
else:
    d_1 =   p_blip[0]-p_blip[2]
    d_2 =   p_blip[1]-p_blip[3]
    if(d_1<=d_2):
        d   =   d_1
    else:
        d   =   d_2

if(n_flag==0):
    d1  =   d
    d2  =   0
else:
    d1  =   0
    d2  =   d
    
print 'status, n_flag set'
print d1,d2
        
buff_1      =   np.zeros(1024,dtype='complex')
buff_2      =   np.zeros(1024,dtype='complex')
x_fftw_1    =   pyfftw.empty_aligned(512,dtype='complex')
y_fftw_1    =   pyfftw.empty_aligned(512,dtype='complex')
x_fftw_2    =   pyfftw.empty_aligned(512,dtype='complex')
y_fftw_2    =   pyfftw.empty_aligned(512,dtype='complex')
x_fftw_1    =   pyfftw.FFTW(x_fftw_1,x_fftw_1)
y_fftw_1    =   pyfftw.FFTW(y_fftw_1,y_fftw_1)
x_fftw_2    =   pyfftw.FFTW(x_fftw_2,x_fftw_2)
y_fftw_2    =   pyfftw.FFTW(y_fftw_2,y_fftw_2)
XX          =   np.zeros((g_l,256),dtype='complex')
YY          =   np.zeros((g_l,256),dtype='complex')
XY          =   np.zeros((g_l,256),dtype='complex')
YX          =   np.zeros((g_l,256),dtype='complex')


f1.read(int(28+(1056*d1)))
f2.read(int(28+(1056*d2)))
q           =   0
p_num_1     =   np.zeros(2)
p_num_2     =   np.zeros(2)
p_num_1[q%2]=   int(str(f1.read(4)).encode('hex'),16)
p_num_2[q%2]=   int(str(f2.read(4)).encode('hex'),16)

print 'Starting Corr'
print p_num_1[q%2],p_num_2[q%2]

for i in range(g_l):
    for j in range(1024*a_count):
        if(i!=0 and (not(trans(p_num_1)) or not(trans(p_num_2)))):
            j               +=  1024
            q               +=  1
            f1.read(1052)
            f2.read(1052)
            p_num_1[q%2]    =   int(str(f1.read(4)).encode('hex'))
            p_num_2[q%2]    =   int(str(f2.read(4)).encode('hex'))
            print 'slip at:',p_num_1[q%2],p_num_2[q%2]
            continue
        else:
            buff_1[j%1024]  =   l_t[int(str(f1.read(1)).encode('hex'),16)]
            buff_2[j%1024]  =   l_t[int(str(f2.read(1)).encode('hex'),16)]
            if(j%1024==1023):    
                X_1         =   x_fftw_1(buff_1[::2])
                Y_1         =   y_fftw_1(buff_1[1::2])
                X_2         =   x_fftw_2(buff_2[::2])
                Y_2         =   y_fftw_2(buff_2[1::2])

                XX[i]       +=   X_1[:256]*np.conj(X_2[:256])
                YY[i]       +=   Y_1[:256]*np.conj(Y_2[:256])
                XY[i]       +=   X_1[:256]*np.conj(Y_2[:256])
                YX[i]       +=   Y_1[:256]*np.conj(X_2[:256])

                f1.read(28)
                f2.read(28)
                q           +=  1
                p_num_1[q%2]=   int(str(f1.read(4)).encode('hex'))
                p_num_2[q%2]=   int(str(f2.read(4)).encode('hex'))


    XX[i][0]    =   XX[i][1]
    XY[i][0]    =   XY[i][1]
    YY[i][0]    =   YY[i][1]
    YX[i][0]    =   YX[i][1]

fig=plt.figure(1)
plt.subplot(221)
plt.imshow(XX)
plt.ylabel('XX')
plt.subplot(222)
plt.imshow(XY)
plt.ylabel('XY')
plt.subplot(223)
plt.imshow(YX)
plt.ylabel('YX')
plt.subplot(224)
plt.imshow(YY)
plt.ylabel('YY')
fig.savefig(file_name_1[:-4]+'=X='+file_name_2[:-4]+'.png')
plt.show()
exit()
