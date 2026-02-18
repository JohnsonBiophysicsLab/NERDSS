#!/usr/bin/env python3
# usage: pull_complex $complex_name
# where $complex_name must exist in ComplexHistogram.dat

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import pandas as pd
import sys
import os





filepath = sys.argv[1]
box_x=sys.argv[2]
box_y=sys.argv[3]
box_z=sys.argv[4]


print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()
  
inum={}
atom="ATOM"
ninety=90.0
nint = 2
#half the box length
hbx=float(box_x)*0.5
hby=float(box_y)*0.5
hbz=float(box_z)*0.5

with open(filepath) as fp:
    icnt = 0
    t=0
    for line in fp:
#        print("line {} contents {}".format(icnt, line))
        if(line.find('iter') != -1):  #true
            arr=line.split(' ')
            justnum=arr[1].split('\n')
            iter=justnum[0]
            fname="%d.pdb" % (int(iter))
            outfile1=open(fname,"w")
            outfile1.write("TITLE  PDB TIMESTEP {}\n".format(iter))
            outfile1.write("CRYST1 {:8.3f}{:8.3f}{:8.3f}{:8.2f}{:8.2f}{:8.2f}    P1\n".format(float(box_x),float(box_y),float(box_z),ninety, ninety, ninety))

            outfile1.write("{:6s}{:5d}  {:^3s}{:3s}{:1s} {:4d}{:1s}    {:8.3f}{:8.3f}{:8.3f}     {:1s}    {:>2s}{:2s}\n".format(atom,0,"A","","A",0,"",float(box_x),float(box_y),float(box_z),"0","0","CL"))
            outfile1.write("{:6s}{:5d}  {:^3s}{:3s}{:1s} {:4d}{:1s}    {:8.3f}{:8.3f}{:8.3f}     {:1s}    {:>2s}{:2s}\n".format(atom,1,"B","","B",1,"",float(box_x),float(box_y),float(box_z),"0","0","CL"))
            outfile1.write("{:6s}{:5d}  {:^3s}{:3s}{:1s} {:4d}{:1s}    {:8.3f}{:8.3f}{:8.3f}     {:1s}    {:>2s}{:2s}\n".format(atom,2,"C","","C",2,"",float(box_x),float(box_y),float(box_z),"0","0","CL"))
            #            print(justnum[0])
            icnt +=1
            t=3 #start new file
            ts=3
            res=-1 #increment by ninterfaces
        else:
           arr=line.split(' ')
           print(" length of line: {} ".format(len(arr)))
           if(len(arr) > 3):
               #this is not the molecule number.
               name=arr[0]   
               x=arr[1]
               y=arr[2]
               z=arr[3]
               if( (t-ts) % nint == 0):
                   res=res+1
#               outfile1.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(atom,t,"","",name,name,res,"",float(x),float(y),float(z),0,0,"0","CL"))
               outfile1.write("{:6s}{:5d}  {:^3s}{:3s}{:1s} {:4d}{:1s}    {:8.3f}{:8.3f}{:8.3f}     {:1s}    {:>2s}{:2s}\n".format(atom,t,name,"",name,res,"",float(x)+hbx,float(y)+hby,float(z)+hbz,"0","0","CL"))

               t=t+1
#               outfile1.write("{} {} {} {} \n".format(name,x,y,z))
           else:
               values=line.split(' ')
               nmols=values[0]    
               print("Number of elements {}".format(nmols));
            

   
