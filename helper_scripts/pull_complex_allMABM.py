#!/usr/bin/env python3
# usage: pull_complex $complex_name
# where $complex_name must exist in ComplexHistogram.dat

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import pandas as pd
import sys
import os





filepath = sys.argv[1]
srank=sys.argv[2]
rank=int(srank)
fname="memloc_rank%d.dat" % (rank)
outfile1=open(fname,"w")

print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()
  
inum={}
cnum1={}
cnum2={}
cnum3={}
cnum4={}
cnum5={}
cnum6={}
cnum7={}
cnum8={}
with open(filepath) as fp:
    icnt = 0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))
        
        if(line.find('iter') != -1):
            arr=line.split(' ')
            justnum=arr[1].split('\n')
            inum[icnt]=justnum[0]
            cnum1[icnt]=0
            cnum2[icnt]=0
            cnum3[icnt]=0
            cnum4[icnt]=0
            cnum5[icnt]=0
            cnum6[icnt]=0
            cnum7[icnt]=0
            cnum8[icnt]=0
#            print(justnum[0])
            icnt +=1
        if(line.find('A: 1. B: 1. M: 2') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) > 6):  #this means it is only B: 1. \n, so length=3.
               cnum1[icnt-1]=arr[0]
#               print("found desired line {}, name length: {} ".format(name, len(name)))
        if(line.find('A: 1. B: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 6):  #this means it is only B: 1. \n, so length=3.
               cnum2[icnt-1]=arr[0]
               #print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('A: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 4):  #this means it is only B: 1. \n, so length=3.
               cnum5[icnt-1]=arr[0]
               #print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('B: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 4):  #this means it is only B: 1. \n, so length=3.
               cnum6[icnt-1]=arr[0]
               #print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('M: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 4):  #this means it is only B: 1. \n, so length=3.
               cnum8[icnt-1]=arr[0]
               #print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('B: 1. A: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 6):  #this means it is only B: 1. \n, so length=3.
               cnum2[icnt-1]=arr[0]
               #print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('A: 1. B: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 6):  #this means it is only B: 1. \n, so length=3.
               cnum2[icnt-1]=arr[0]
               #print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('A: 1. M: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 6):  #line "B: 1. \n" has length=3.
               cnum3[icnt-1]=arr[0]
 #              print("found desired line {}, name length: {} ".format(name, len(name)))
           
               
        if(line.find('B: 1. M: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           if(len(name) < 6):  #line "B: 1. \n" has length=3.
               cnum4[icnt-1]=arr[0]
#               print("found desired line {}, name length: {} ".format(name, len(name)))

        if(line.find('A: 1. B: 1. M: 1.') !=-1):
           arr=line.split('\t')
           name=arr[1].split(' ')
#           print("found line {}, name length: {} ".format(name, len(name)))
           cnum7[icnt-1]=arr[0] #in this case, it has LIP, +B and A
#           print(arr[0])

#print(" 2 columns!")
outfile1.write("iter MABM AB AM BM A B ABM M\n")
for i, iters in enumerate(inum):
    step=inum[i]
    v1=cnum1[i]
    v2=cnum2[i]
    v3=cnum3[i]
    v4=cnum4[i]
    v5=cnum5[i]
    v6=cnum6[i]
    v7=cnum7[i]
    v8=cnum8[i]
    outfile1.write("{} {} {} {} {} {} {} {} {}\n".format(step,v1,v2,v3,v4,v5,v6,v7,v8))



   
