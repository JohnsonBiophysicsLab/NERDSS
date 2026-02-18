#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:37:07 2020

Used to fit the pucadyilData by exp_lag

@author: sikao
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#define the function exp_lag
def exp_lag(x, a, b, c, d):
    return d + np.heaviside(x - c, 0.5) * a * (1.0 - np.exp(-(x - c) / b))

needCreatTotal = True

for itrDirectory in range(1,6):
    print("/Users/sikao/Downloads/Simulations/1001/case1/data/%d/histogram_complexes_time.dat" % (itrDirectory))
    # create copies table for copy numbers of clathrin on membrane
    #iter copies  this line stored in headline
    #0    1       else is in output
    #1    2
    ##########################
    headline = ['time']
    headline.append('copies')
    output = []
    tempLine = []
    clathrinOnMem = 0
    clathrinInSol = 0
    
    # Load the ComplexHistoram file
    f = open("/Users/sikao/Downloads/Simulations/1001/case1/data/%d/histogram_complexes_time.dat" % (itrDirectory))
    line = f.readline()
    if line[0] == 'T': # this is a time line
        # we need to find the position of space, it is at index 5
        
        # then we get the time
        iterStep = float(line[10:])
        
        #update the copies table
        tempLine.append(iterStep) 
        
    else: # this is a line contained the copies and name of a complex, seperated by \t
        # we get the copies  and name first, need to find the position of the \t
        firstSpacePosition = line.find('\t')
        copies = int(line[:firstSpacePosition])
        name = line[firstSpacePosition+1:len(line)-1]
        
        # analyze the components of the complex
        if name[0] == 'c' and name.find('I') != -1: #this is a clathrin struct on membrane
            # get the number of clathrin in one complex, between the : and .
            sizeStruct = int(name[name.find(' ') + 1 : name.find('.')])
            clathrinOnMem += sizeStruct * copies
            
        if name[0] == 'c' and name.find('I') == -1: #this is a clathrin struct in solution
            # get the number of clathrin in one complex, between the : and .
            sizeStruct = int(name[name.find(' ') + 1 : name.find('.')])
            clathrinInSol += sizeStruct * copies
        
    while line:
        line = f.readline()
        if len(line) == 0:
            continue
        if line[0] == 'T': # this is a new iter line, append tempLine to output first
            tempLine.append(clathrinOnMem)
            tempLine.append(clathrinInSol)
            clathrinOnMem = 0
            clathrinInSol = 0
            output.append(tempLine)
            tempLine = []
            
            # then we get the iter step number
            iterStep = float(line[10:])
        
            #update the copies table
            tempLine.append(iterStep)
                    
        else: # this is a line contained the copies and name of a complex, seperated by \t
            # we get the copies  and name first, need to find the position of the \t
            firstSpacePosition = line.find('\t')
            copies = int(line[:firstSpacePosition])
            name = line[firstSpacePosition+1:len(line)-1]
        
            # analyze the components of the complex
            if name[0] == 'c' and name.find('I') != -1: #this is a clathrin struct on membrane
                # get the number of clathrin in one complex, between the : and .
                sizeStruct = int(name[name.find(' ') + 1 : name.find('.')])
                clathrinOnMem += sizeStruct * copies
                
            if name[0] == 'c' and name.find('I') == -1: #this is a clathrin struct in solution
                # get the number of clathrin in one complex, between the : and .
                sizeStruct = int(name[name.find(' ') + 1 : name.find('.')])
                clathrinInSol += sizeStruct * copies
                
                            
    f.close()
    
    clathrinOnMemTable = [headline] + output
    
    tmpResult = np.array(output)
    # write result of each trace to file
    np.savetxt("/Users/sikao/Downloads/Simulations/1001/case1/data/%d/sim_time.dat" % (itrDirectory), tmpResult[:,0])
    np.savetxt("/Users/sikao/Downloads/Simulations/1001/case1/data/%d/sim_clatMem.dat" % (itrDirectory), tmpResult[:,1])
    np.savetxt("/Users/sikao/Downloads/Simulations/1001/case1/data/%d/sim_clatSol.dat" % (itrDirectory), tmpResult[:,2])
    
#    if needCreatTotal == True:
#        tmpTotal = np.array(output)
#        needCreatTotal = False
#    else:
#        # align
#        minLength = min(len(tmpResult), len(tmpTotal))
#        tmpTotal = tmpTotal[:minLength, :]
#        tmpResult = tmpResult[:minLength, :]
#        
#        # append
#        tmpTotal = np.append(tmpTotal, tmpResult, axis=1)
#        
#simDataTime = tmpTotal[:,0] * timeStep / 1e6
#simDataAll = tmpTotal[:,1::2]
#simDataMean = np.mean(simDataAll, axis=1)
#
#plt.plot(simDataTime, simDataMean, 'b-', label='simDataMean')
#
#popt, pcov = curve_fit(exp_lag, simDataTime, simDataMean, bounds=(0, [10e3, 2000.0, 20.0, 50.0]))
#
#x_plot_fit = np.linspace(0, 100, 1000)
#plt.plot(x_plot_fit, exp_lag(x_plot_fit, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt))
#
#plt.xlabel('x')
#plt.ylabel('y')
#plt.legend()
#plt.show()