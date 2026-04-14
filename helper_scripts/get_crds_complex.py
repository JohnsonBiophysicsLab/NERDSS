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

#srank=sys.argv[2]
#rank=int(srank)
fname="crds.inp"
fname2="complex1.inp"
fname3="complex.inp"
#% (rank)

outfile1=open(fname,"w")
outfile2=open(fname2,"w")
outfile3=open(fname3,"w")

print("File path {}: ".format(filepath))
if not os.path.isfile(filepath):
    print("File path {} does not exist. Exiting...".format(filepath))
    sys.exit()




with open(filepath) as fp:
    start=3
    cstart =3
    nline=0
    for line in fp:

        if(line.find('Molecules') != -1):
            #start writing after next line
            start=0
        else:
            start +=1
            
        if(line.find('Complexes') != -1):
            #start writing after next line to other file
            start=3 #stop writing to outfile
            cstart=0
        else:
            cstart +=1

        if(line.find('Observables') != -1):
            cstart=3 #stop writing to outfile
            

        if(start == 2):
            outfile1.write("{}".format(line))
            start = 1 #keep coming back to this loop from now on.

        if(cstart == 2):
            outfile2.write("{}".format(line))
            nline +=1
            cstart = 1 #keep coming back to this loop from now on.

    keepline=nline-1
    print 'n lines with complexes to keep: ', keepline
    print(keepline)

    
newfilepath="complex1.inp"
with open(newfilepath) as fp2:

    nline=0
    for line in fp2:
        if(nline<keepline):
            outfile3.write("{}".format(line))
            nline +=1
