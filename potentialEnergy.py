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
            if(diffX>5):
                diffX = diffX - 10
            if(diffX<-5):
                diffX = diffX + 10
            diffY = y2-y1
            if(diffY>5):
                diffY = diffY - 10
            if(diffY<-5):
                diffY = diffY + 10
            diffZ = z2-z1
            if(diffZ>5):
                diffZ = diffZ - 10
            if(diffZ<-5):
                diffZ = diffZ + 10

            #find r
            r = math.sqrt((diffX**2)+(diffY**2)+(diffZ**2))
            totalR = totalR + r

            tempPE = 4.0*(((1/r)**12)-((1/r)**6))
            
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
    #for indx in range(len(finalR)):
    #    outFile.write(str(finalR[indx]) + ' PE:')
    #    outFile.write(str(finalPE[indx])+'\n')
    #outFile.close()
    
    #plot potential energy at each time step
    plt.plot(finalR,finalPE)
    plt.xlabel('Distance')
    plt.ylabel('Potential Energy (eV)')
    plt.show()

#run functions    
main()


