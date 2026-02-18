#!/usr/bin/env python3                                                                                                                                            
# usage: remove_molecules.py $crds.inp $just_complex.inp $Original_ReducedMolecul.inp Original_ReducedComplex.inp
# where $crds.inp can only contain coordinates.
#just_complex is from that section of the restart, excluding the last Nempty line
#Original_.inp files contain two columns, original protein num, new renumbered index
#Original_complex.inp is the same, contains original complex num, new renumbered index
#formatting must have the index on first line, and second line starts with '1.000000000000'
#for each new molecule.
#renumber both molecules and complexes, and then match them based on these new numbers.
#WARNING: Be CAREFUL with printing out last element of molecule array and complex array
#Double check this is correct!
#WARNING: SETS N Empty complex=0 for last line. 

import sys
import os





filepath = sys.argv[1]
filepath2 = sys.argv[2]
#srank=sys.argv[2]
#rank=int(srank)
fname="UPDATED_RENUM_cleancrds.dat"
fname2="DELETED_crds.dat"
#% (rank)

outfile=open(fname,"w")
delfile=open(fname2,"w")

print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()

inum={}
delnum={}
orig_to_reduce={}
reduce_to_orig={}

origpro={}
origcomplex={}
filepath3 = sys.argv[3]
n=0
with open(filepath3) as fp3:
    for line in fp3:
        n+=1
        proarr=line.strip()
        proindices=line.split(' ')
        origpro[int(proindices[0])]=int(proindices[1])
        print("orig pro from Omolecule: {}",format(proindices[0]))

c=0
filepath4 = sys.argv[4]
with open(filepath4) as fp4:
    for line in fp4:
        c+=1
        proarr=line.strip()
        proindices=line.split(' ')
        origcomplex[int(proindices[0])]=int(proindices[1])


print "Length ",n
totMol=n
totComplex=c

with open(filepath) as fp:
    i = 0
    molnum=-1
    keepnum=0
    linestart=0;
    check='false'
    delete = 'true'
    nlipid=0
    d=0
    outfile.write("#All Molecules and coordinates\n")
    outfile.write("{} {}\n".format(totMol, totMol))
    for line in fp:
        print("line contents {}".format(line))

        #line1:line7, store.
        #line 8 lists U (state)
        # line 9, next interface
        inum[i]=line.strip();
        arr=line.split(' ');
        if(arr[0]=='1.000000000000'):
            #print out the previous lines, if they are not U
            if(delete == 'false'):
                #print all lines
#                outfile.write("Molecule {} Keep {}\n".format(molnum, keepnum))
                
                
                firstline=inum[linestart].split(' ');
                print "first element ",firstline[0]
                molnum=int(firstline[0])
                orig_to_reduce[molnum]=keepnum
                reduce_to_orig[keepnum]=molnum
                keepnum+=1 #this tracks only printed molecules.
                
                newpro=origpro[int(firstline[0])]
                newcomplex=origcomplex[int(firstline[2])]
                outfile.write("{} {} {} {} {}\n".format(newpro, firstline[1], newcomplex, firstline[3], firstline[4]))
                currline=linestart+1
                for j in range(currline, currline+4):
                    outfile.write("{} \n".format(inum[j]))
                currline+=4
                #now write out Nbnd partner, updating protein indexes if they are bound.
                prolabels=inum[currline].split(' ');
                outfile.write("{}".format(prolabels[0]))
                for j in range(1, int(prolabels[0])+1):
                    outfile.write(" {}".format(origpro[int(prolabels[j])]))
                outfile.write("\n")
                currline+=1
                nifaces=int(inum[currline])
                outfile.write("{} \n".format(inum[currline])) #print niface to file.
                currline+=1
                for j in range(0, nifaces):
                    ifaceline=inum[currline].split(' ')
                    isbound=int(ifaceline[5])
                    outfile.write("{} \n".format(inum[currline])) #print isbound line
                    #next line is crds
                    currline+=1
                    outfile.write("{} \n".format(inum[currline])) #print coords to file
                    currline+=1
                    if(isbound == 1):
                        bndlist=inum[currline].split(' ')
                        num=int(bndlist[0])
                        outfile.write("{} {} {} \n".format(origpro[num], bndlist[1], bndlist[2]))
                        currline+=1
                    #end if bound
                #end all ifaces
                for j in range(currline, i-1):
                    outfile.write("{} \n".format(inum[j]))

                #rest of lines
            else:
                delfile.write("Deleted {}\n".format(molnum))
                if(molnum>-1):
                    delnum[d]=molnum
                    d=d+1
            molnum+=1
            if(molnum==6409):
                print "Last molecule at line ", i
      
            delete = 'false' #reset this for each new molecule.
            linestart=i-1
            #now check if it is a U
            if(arr[1]=='1'):
                #this is a lipid
                
                nlipid+=1
                print('found a lipid')
           

        #now figure out if it is U.
        if(line.find('U') != -1):
            #we want to delete this molecule, and all its information.
            delete='true'
        i+=1 # increment lines
                



#done looping over all lines, print out the last protein you read in
#outfile.write(" LAST ELEMENT\n")
if(delete == 'false'):
    #print all lines
    #                outfile.write("Molecule {} Keep {}\n".format(molnum, keepnum))
    
    
    firstline=inum[linestart].split(' ');
    print "first element ",firstline[0]
    molnum=int(firstline[0])
    orig_to_reduce[molnum]=keepnum
    reduce_to_orig[keepnum]=molnum
    keepnum+=1 #this tracks only printed molecules.
                
    newpro=origpro[int(firstline[0])]
    newcomplex=origcomplex[int(firstline[2])]
    outfile.write("{} {} {} {} {}\n".format(newpro, firstline[1], newcomplex, firstline[3], firstline[4]))
    currline=linestart+1
    for j in range(currline, currline+4):
        outfile.write("{} \n".format(inum[j]))
    currline+=4
    #now write out Nbnd partner, updating protein indexes if they are bound.
    prolabels=inum[currline].split(' ');
    outfile.write("{}".format(prolabels[0]))
    for j in range(1, int(prolabels[0])+1):
        outfile.write(" {}".format(origpro[int(prolabels[j])]))
    outfile.write("\n")
    currline+=1
 #   outfile.write("currline: ".format(inum[currline]))
    nifaces=int(inum[currline])
    outfile.write("{} \n".format(inum[currline])) #print niface to file.
    currline+=1
    for j in range(0, nifaces):
        ifaceline=inum[currline].split(' ')
        isbound=int(ifaceline[5])
        outfile.write("{} \n".format(inum[currline])) #print isbound line
        #next line is crds
        currline+=1
        outfile.write("{} \n".format(inum[currline])) #print coords to file
        currline+=1
        if(isbound == 1):
            bndlist=inum[currline].split(' ')
            num=int(bndlist[0])
            outfile.write("{} {} {} \n".format(origpro[num], bndlist[1], bndlist[2]))
            currline+=1
        #end if bound
    #end all ifaces
    for j in range(currline, i-1):
        outfile.write("{} \n".format(inum[j]))
    
    #rest of lines

else:
    delfile.write("Deleted {}\n".format(molnum))
    delnum[d]=molnum
    d=d+1
molnum+=1

########
#Before printing complex, there is another element: N empty Molecules
#inum.reverse()
outfile.write("{}\n".format(inum[j])) #last element

print('Total molecules read in') 
print(molnum)
molnumsave=molnum
print ('total lipids found ')
print(nlipid)
nlipidsave=nlipid
print ('total deleted ')
totdel=d
print(totdel)
print ('total kept ')
totKept=keepnum  
print(totKept)

fname3="Original_ReducedMolecule.dat"
rfile=open(fname3,"w")
for t, vals in enumerate(reduce_to_orig):
    molnum=reduce_to_orig[t]
    num=int(reduce_to_orig[t])
    rfile.write("{} {}\n".format(molnum, orig_to_reduce[num]))

    
    


        
#now loop over all complexes.
print "read in file: ", filepath2
old_to_reduce={}
reduce_to_old={}
outfile.write("#All Complexes and components\n")
outfile.write("{} {}\n".format(totComplex, totComplex))
with open(filepath2) as fp2:
    i = 0
    complexNum=0
    linestart=0;
    check='false'
    delete = 'true'
    nlipid=0
    d=0
    remove='false'
    keepnum=0
    for line in fp2:
#        print("line {} contents {}".format(icnt, line))

        inum[i]=line.strip();
#        print("line {} contents {}".format(i, line))
  
        if(i==4):
            #this is the number of members and their indices
            arr=line.split(' ');
            num=int(arr[0])
            print "complex number ",complexNum
            print "Nmembers ",num
            if(num >0 ):
                print "member 1 ",arr[1]
            
                pro=int(arr[1])
                print "protein index ",pro
                #print "total delete ",totdel
                for dp in range(0, totdel):
                    #   print "deleted proteins ",delnum[dp]
                    if(pro == int(delnum[dp])):
                        #this protein has been deleted.
                        remove = 'true'
                #done for loop
            else:
                remove = 'true' #remove empty complexes.
                
        if(i==5):
            #this is the last line for this complex, write or not.
            #print all lines
            if(remove=='false'):
#                outfile.write("Complex {}\n".format(complexNum))
                reduce_to_old[keepnum]=complexNum
                old_to_reduce[complexNum]=keepnum
                keepnum+=1
                comlabel=inum[0].split(' ')
                prolabels=inum[4].split(' ')
                newcomplex=origcomplex[int(comlabel[0])]
                outfile.write("{} {} {} {} \n".format(newcomplex, comlabel[1], comlabel[2], comlabel[3]))
                for j in range(1, 4):
                    outfile.write("{} \n".format(inum[j]))
                outfile.write("{}".format(prolabels[0]))
                for j in range(1, int(prolabels[0])+1):
                    outfile.write(" {}".format(origpro[int(prolabels[j])]))
                outfile.write("\n")
                
                outfile.write("{} \n".format(inum[5]))
            else:
                delfile.write("Complex Deleted {}\n".format(complexNum))
            remove='false' #reset
            i=-1 #restart the counter.
            complexNum+=1
        #done
        i=i+1
    #done loop over lines

#Now Do last line
#Comment below out, it is repeating the last complex
#arr=inum[4].split(' ');
#num=int(arr[0])
#print "complex number ",complexNum
#print "Nmembers ",num
#if(num >0 ):
#    print "member 1 ",arr[1]
#    
#    pro=int(arr[1])
#    print "protein index ",pro
#    #print "total delete ",totdel
#    for dp in range(0, totdel):
#        #   print "deleted proteins ",delnum[dp]
#        if(pro == int(delnum[dp])):
#            #this protein has been deleted.
#            remove = 'true'
#            #done for loop
#else:
#    remove = 'true' #remove empty complexes.
remove='true'
if(remove =='false'):
#    outfile.write("Complex {}\n".format(complexNum))
    reduce_to_old[keepnum]=complexNum
    old_to_reduce[complexNum]=keepnum
    keepnum+=1
    comlabel=inum[0].split(' ')
    prolabels=inum[4].split(' ')
    print "last line",inum[0], comlabel[0]
    newcomplex=origcomplex[int(comlabel[0])]
    outfile.write("{} {} {} {} \n".format(newcomplex, comlabel[1], comlabel[2], comlabel[3]))
    for j in range(1, 4):
        outfile.write("{} \n".format(inum[j]))
    outfile.write("{}".format(prolabels[0]))
    for j in range(1, int(prolabels[0])+1):
        outfile.write(" {}".format(origpro[int(prolabels[j])]))
    outfile.write("\n")
    outfile.write("{} \n".format(inum[5]))
else:
    delfile.write("Complex Deleted {}\n".format(complexNum))

####
###print out last element: N empty complexes
#inum.reverse()
outfile.write("0 \n")

print "Total molecules read in ",molnumsave

print "total lipids found ",nlipidsave

print "total deleted ",totdel
print "total deleted: (should be same repeate) ", (molnumsave-totKept)

print "total molecules kept ",totKept
lipidLeft=nlipidsave-totdel
print " Total lipids remaining (should be same repeated)", lipidLeft
print "total complexes kept ", keepnum

totAtomClApSy=1930 #300*4+100*7+10*3

totAtomSys=totAtomClApSy+lipidLeft*2
print "n total site units ",totAtomSys

fname4="Original_thenReducedComplex.dat"
rfile2=open(fname4,"w")
for t, vals in enumerate(reduce_to_old):
    molnum=reduce_to_old[t]
    num=int(reduce_to_old[t])
    rfile2.write("{} {} \n".format(molnum, old_to_reduce[num]))
        
