#! Program Files (x86)/PythonDownload/python
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import statistics

def determinePE(xPos,yPos,zPos):
    
    nAtoms = len(xPos)
    allR = []
    allPE = []
    totalPE = 0
    totalR = 0

    for i in range(nAtoms):
        for j in range(i+1,nAtoms):
            #only execute if two different atoms
            #positions
            x1 = xPos[i]
            x2 = xPos[j]
            y1 = yPos[i]
            y2 = yPos[j]
            z1 = zPos[i]
            z2 = zPos[j]
 
            #difference in pos
            diffX = x2-x1
            diffY = y2-y1
            diffZ = z2-z1

            #find r
            r = math.sqrt((diffX**2)+(diffY**2)+(diffZ**2))
            totalR = totalR + r

            #constants
            Ar_diam = 142.0 #nm
            #from White 1999
            eps = 1.0 #125.0*(8.6173*10**-5) #eV
            sigma = 0.3345 #nm

            rCorrected = r*sigma
            sigOverR = sigma/rCorrected
            tempPE = 4.0*eps*((sigOverR**12)-(sigOverR**6))
            
            #add to total potential
            totalPE = totalPE + tempPE
        
        avgPE = totalPE/(nAtoms-1)
        allPE.append(avgPE)
        avgR = totalR/(nAtoms-1)
        allR.append(avgR)
        totalPE = 0
        totalR = 0
    
    overallAvgPE = statistics.mean(allPE)
    overallAvgR = statistics.mean(allR)

    return overallAvgPE, overallAvgR
    

def main():
    #specify where the text file is located
    filePath = open('/Users/gbonn/Summer_Research_2020/lammps_tut/testCoords.txt','r')
    finalPE = []
    finalR = []
    xData = []
    yData = []
    zData = []
    
    for lineNum, line in enumerate(filePath):
        testLine = lineNum%502
        if(testLine != 0 and testLine!=1):
            lineList = line.split()
            xData.append(float(lineList[1]))
            yData.append(float(lineList[2]))
            zData.append(float(lineList[3]))
            if(testLine == 501):
                tempPE,tempR = determinePE(xData,yData,zData)
                finalPE.append(tempPE)
                finalR.append(tempR)
        
    filePath.close()
    #outFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/testOutput.txt','w')
    #for indx in len(finalR):
     #   outFile.write(finalR[indx] + ' PE:')
     #   outFile.write(finalPE[indx]+'\n')
    #outFile.close()
    
    #plot potential energy at each time step
    plt.plot(finalR,finalPE)
    plt.xlabel('Distance')
    plt.ylabel('Potential Energy (eV)')
    plt.show()

#run functions    
main()


