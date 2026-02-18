#!/usr/bin/env python3
# usage: concatenate_files $file1 $file2
# where files have same format and first col is iterations/time.

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import pandas as pd
import sys
import os





filepath = sys.argv[1]
filepath2 = sys.argv[2]
srank=sys.argv[3]
rank=int(srank)
fname="combined_rank%d.dat" % (rank)
outfile1=open(fname,"w")

print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()
  
iter1={}
iter2={}
with open(filepath) as fp:
    icnt = 0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))
        
        if(line.find('MABM') == -1): #this checks for the find being false, then enters if statement
            arr=line.split(' ')
#            print(" file 1: {}".format(arr[0]))
            iter1[icnt]=arr[0]
            icnt +=1
            

lastiter=arr[0]
print("last iteration: {} ".format(lastiter))


with open(filepath2) as fp:
    icnt = 0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))
        
        if(line.find('MABM') == -1):
            
            arr=line.split(' ')
#            print(" file 2: {}".format(arr[0]))     
            iter2[icnt]=arr[0]
            icnt +=1

startiter=iter2[0]
print(" first iter next file: {}".format(startiter))

#read them in again and print up to startiter only

with open(filepath) as fp:
    icnt = 0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))
        if(line.find('MABM') == -1):    
            arr=line.split(' ')
            if(int(arr[0]) < int(startiter) ):
#                print("write out iter from file 1 {}".format(arr[0]))
                outfile1.write(line)
        else:
            outfile1.write(line)


#now read in file2
with open(filepath2) as fp:
    icnt = 0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))
        
        if(line.find('MABM') == -1):
            arr=line.split(' ')
            outfile1.write(line)



   
