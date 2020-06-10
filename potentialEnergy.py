#! Program Files (x86)/PythonDownload/python
import matplotlib.pyplot as plt
import numpy as np
import math

def determinePE(xPos,yPos,zPos):
    
    nAtoms = len(xPos)
    allR = []
    allPE = []

    for i in range(nAtoms-1):
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
            rSq = (diffX**2)+(diffY**2)+(diffZ**2)
            allR.append(rSq)
            tempPE = 4.0*(((1/rSq)**6)-((1/rSq)**3))
            allPE.append(tempPE)
        
    allR = np.sqrt(allR)
    return allR, allPE
    
def binning(allRValues,allPEValues,numBins):
    numBins = float(numBins)
    binSpace = math.ceil((math.sqrt(10.0*math.sqrt(3))/(numBins)))
    rPEPairs = [[]]
    currR = 0
    for indx in range(numBins):
        rPEPairs[indx][0] = currR
        currR = currR + binSpace
    for timeStep in range(len(allRValues)):
        for atomRIndx in range(len(allRValues[timeStep])):
            tempR = allRValues[timeStep][atomRIndx]
            binIndx = -1
            while (tempR < rPEPairs[binIndx]):
                binIndx = binIndx - 1
            tempPE = allPEValues[timeStep][atomRIndx]
            if(len(rPEPairs[binIndx]) == 1):
                rPEPairs[binIndx][1] = tempPE
                rPEPairs[binIndx][2] = 1.0
            else:
                totPrevPE = rPEPairs[binIndx][1]*rPEPairs[binIndx][2]
                rPEPairs[binIndx][2] = rPEPairs[binIndx][2] + 1.0
                rPEPairs[binIndx][1] = (totPrevPE + tempPE)/rPEPairs[binIndx][2]


    return rPEPairs

def main():
    #specify where the text file is located
    filePath = open('/Users/gbonn/Summer_Research_2020/lammps_tut/melt.lmpdump','r')
    #outFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/rAndPE_pairs.txt','w')
    xData = []
    yData = []
    zData = []
    allMinR = []
    allMaxR = []
    allR = []
    allPE = []
    
    for lineNum, line in enumerate(filePath):
        testLine = lineNum%509
        if(testLine not in range(9)):
            lineList = line.split()
            xData.append(float(lineList[2]))
            yData.append(float(lineList[3]))
            zData.append(float(lineList[4]))
            if(testLine == 508):
                #print(xData[-1],yData[-1],zData[-1])
                tempR,tempPE = determinePE(xData,yData,zData)
                allR.append(tempR)
                allPE.append(tempPE)
                xData = []
                yData = []
                zData = []
                #if(lineNum > 508):
                    #outFile.write(str(tempR) + " " + str(tempPE) + "\n")
                if(lineNum >= 500):
                    break
    
    finalPairs = binning(allR,allPE,10)
    print(finalPairs)
    #print(currMin,currMax)
    filePath.close()
    #outFile.close()
    
    #plot potential energy at each time step
    #print(len(allR[0]),len(allPE[0]))
    #plt.plot(allR[0],allPE[0],'ro',ms = 1)
    #plt.xlim(0,10)
    #plt.ylim(-2,15)
    #plt.plot(allR[1],allPE[1],'ro')
    #plt.xlabel('Distance')
    #plt.ylabel('Potential Energy')
    #plt.show()

#run functions    
main()


