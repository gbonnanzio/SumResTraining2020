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
    
def binning(timeStepRValues,timeStepPEValues,allBinned,numBins):
    #takes current r and PE values for one time step
    #adds them in to the overall bin created
    binnedR = []
    binnedPE = []
    minR = 0.9
    maxR = 5.0*math.sqrt(3)
    binSpace = (maxR-minR)/(numBins - 1)
    rRange = np.linspace(minR, maxR, numBins)
    
    #if we are creating the final bins for the first time
    if(len(allBinned) == 0):
        for indx in range(numBins):
            tmpList = [rRange[indx], 0, 0]
            allBinned.append(tmpList)
        
    for timeStep in range(len(timeStepRValues)):
        tempR = timeStepRValues[timeStep]
        tempPE = timeStepPEValues[timeStep]
        binIndx = -1
        while (tempR < allBinned[binIndx][0]):
            binIndx = binIndx - 1
            if(binIndx == -len(allBinned)):
                binIndx = 0
                break
        #if PE hasn't been added to bin yet
        if(allBinned[binIndx][2] == 0): 
            allBinned[binIndx][1] = tempPE
            allBinned[binIndx][2] = 1.0 #added one PE so far
        #if PE has already been added to this bin before 
        else:
            totPrevPE = allBinned[binIndx][1]*allBinned[binIndx][2]
            allBinned[binIndx][2] = allBinned[binIndx][2] + 1.0
            allBinned[binIndx][1] = (totPrevPE + tempPE)/allBinned[binIndx][2]
            
    
    return allBinned, binSpace

def main():
    #specify where the text file is located
    filePath = open('/Users/gbonn/Summer_Research_2020/lammps_tut/melt.lmpdump','r')
    totalPEOutFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/totalPotEnPerTS.txt','w')
    rdfOutFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/radialDistFun.txt','w')
    xData = []
    yData = []
    zData = []
    binsTotal = []
    binnedR = []
    binnedPE = []
    binWidth = 0

    for lineNum, line in enumerate(filePath):
        testLine = lineNum%509
        if(testLine not in range(9)):
            lineList = line.split()
            xData.append(float(lineList[2]))
            yData.append(float(lineList[3]))
            zData.append(float(lineList[4]))
            if(testLine == 508):
                #print(xData[-1],yData[-1],zData[-1])
                #print(len(binsTotal))
                tempR,tempPE = determinePE(xData,yData,zData)
                totalPEOutFile.write(str(sum(tempPE))+'\n')
                binsTotal, binSpace = binning(tempR,tempPE,binsTotal,150)
                xData = []
                yData = []
                zData = []
                #if(lineNum > 508):
                    #outFile.write(str(tempR) + " " + str(tempPE) + "\n")
                if(lineNum >= 75000):
                    break
    
    #print(len(binsTotal))
    #print(binsTotal)
    #reimannApproxTotalPE = 0
    totalPts = 0
    freq = []
    for indx in range(len(binsTotal)):
        if(binsTotal[indx][2] == 0):
            indx = indx - 1
        else:
            binnedR.append(binsTotal[indx][0])
            binnedPE.append(binsTotal[indx][1])
            #reimannApproxTotalPE = reimannApproxTotalPE + binWidth*binsTotal[indx][1]
            totalPts = totalPts + binsTotal[indx][2]
    
    
    for indx in range(len(binsTotal)):
        currFreq = binsTotal[indx][2]/totalPts
        freq.append(currFreq)
        currR = binsTotal[indx][0]
        rdfOutFile.write(str(currR) + " " + str(currFreq) + "\n")
    
    #print(reimannApproxTotalPE)
    #print(sum(freq))

    filePath.close()
    totalPEOutFile.close()
    rdfOutFile.close()
        
    #plot potential energy at each time step
    
    plt.plot(binnedR,binnedPE,lw = 1)
    #plt.plot(binnedR,freq)
    #plt.plot(allR[2],allPE[2],'ro',ms = 1)
    plt.xlim(0.5,4)
    plt.ylim(-1.5,1.5)
    #plt.plot(allR[1],allPE[1],'ro')
    plt.xlabel('Distance (D)')
    plt.ylabel('Potential Energy (kT)')
    plt.legend(['Average Potential Energy'])
    plt.show()

#run functions    
main()


