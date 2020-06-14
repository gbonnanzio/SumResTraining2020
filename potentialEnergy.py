#! Program Files (x86)/PythonDownload/python
import matplotlib.pyplot as plt
import numpy as np
import math
import sys

################## determinePE ##############################################################

def determinePE(xPos,yPos,zPos):
    
    # takes the coordinates of a certain number of atoms 
    # determines the distance between every possible pair and the Lennard Jones potential 
    # and returns the radii with the pair (assumes periodic boundary conditions)

    nAtoms = len(xPos)
    allR = []
    allPE = []
    
    # j is adjusted so that you never compare any two atoms twice 
    for i in range(nAtoms-1):
        for j in range(i+1,nAtoms):
            x1 = xPos[i]
            x2 = xPos[j]
            y1 = yPos[i]
            y2 = yPos[j]
            z1 = zPos[i]
            z2 = zPos[j]

            # implement PBC for box length 10 
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

            # calculate r and potential energy 
            rSq = (diffX**2)+(diffY**2)+(diffZ**2)
            allR.append(rSq)
            tempPE = 4.0*(((1/rSq)**6)-((1/rSq)**3))
            allPE.append(tempPE)
        
    allR = np.sqrt(allR)
    return allR, allPE


######################## binning ##############################################################

def binning(binSpace, radRange, currRValues,currPEValues,allBinned,numBins):
    
    # takes the current r and PE lists for the current time step 
    # as well as a list of all the past binned r and PE values 
    # and the number of bins 
    # finds what bin r is in and averages in PE to running total
    # returns compiled bin and the specified 
   
    # if we are creating the final bins for the first time
    if(len(allBinned) == 0):
        for indx in range(numBins):
            # adding arrays of format [r, PE, count] to overall array
            # count is number of times a PE has been added to that bin
            tmpList = [radRange[indx], 0, 0]
            allBinned.append(tmpList)
    # bin the current (r,PE) pairs   
    for atomNum in range(len(currRValues)):
        tempR = currRValues[atomNum]
        tempPE = currPEValues[atomNum]
        # iterate from back of list starting at the largest r until
        # we find spot for tempR in bin 
        binIndx = int(tempR//binSpace)
        # if initializing this bin
        if(allBinned[binIndx][2] == 0): 
            allBinned[binIndx][1] = tempPE
            allBinned[binIndx][2] = 2.0 #first PE added, accounts for ij and ji 
        else:
            # calculate new average PE
            totPrevPE = allBinned[binIndx][1]*allBinned[binIndx][2]
            allBinned[binIndx][2] = allBinned[binIndx][2] + 2.0
            allBinned[binIndx][1] = (totPrevPE + 2*tempPE)/allBinned[binIndx][2]
            
    return allBinned

######################### radialDistFun ##########################################################

def radialDistFun(binWidth,allRad,allCts,runs):
    
    # takes a binSize, all of the radii from all time steps, the bin count for each
    # radius and the number of time steps performed
    # creates the radial distribution function g(r)
    radDistVal = []
    n = []
    nID = []
    meanDens = (500)/(1000) #500.0/((4.0/3.0)*math.pi*(max(allRad)**3))
    #print(densList)
    #print(meanDens)
    for i in range(len(allRad)):
        n.append(allCts[i]/(500*runs))
        tempNId = (4.0/3.0)*math.pi*(meanDens)*((allRad[i] + 0.5*binWidth)**3 - (allRad[i]- 0.5*binWidth)**3)
        nID.append(tempNId)
        radDistVal.append(n[i]/nID[i])

    return radDistVal



################################### MAIN ##########################################################

def main():

    #specify where files are located
    coordFile = open('melt.lmpdump','r') #contains all the coords of the atoms at every time stamp
    
    # initialize necessary variables
    xData = []
    yData = []
    zData = []
    binsTotal = [] # has length numBins and contains arrays [r, PE, count]
    binnedR = [] # r bins
    binnedPE = [] # pe bins
    binnedPts = [] # how many pairs we added to each bin
    currTimeStep = 0
    totBins = 100
    minR = 0.0
    maxR = 5.0*math.sqrt(3) # for box length 10 with PBC
    binSize = (maxR-minR)/(totBins - 1)
    rRange = np.linspace(minR, maxR, totBins)

    for lineNum, line in enumerate(coordFile):
        testLine = lineNum%509 #dump file has 9 lines of garbage 500 lines of coords
        if(testLine not in range(9)):#if not a garbage line
            lineList = line.split()
            xData.append(float(lineList[2]))
            yData.append(float(lineList[3]))
            zData.append(float(lineList[4]))
            # if last coordinate in time step, analyze time step
            if(testLine == 508):
                sys.stdout.write(str(currTimeStep)+'\n')
                sys.stdout.flush()
                tempR,tempPE = determinePE(xData,yData,zData)
                binsTotal = binning(binSize,rRange,tempR,tempPE,binsTotal,totBins)
                #numDens.append(timeDens)
                # reset xyz arrays to empty
                xData = []
                yData = []
                zData = []
                currTimeStep = currTimeStep + 1 # how many time steps we've done
                # control how many time steps we analyze
                if(currTimeStep > 3000):
                    break
        
    #close files we don't need
    coordFile.close()

    # organize our binned data
    for indx in range(len(binsTotal)):
        # if we didn't initialize bin, skip it
        if(binsTotal[indx][2] == 0):
            indx = indx - 1
        # split up bin into its three parts
        else:
            binnedR.append(binsTotal[indx][0] + 0.5*binSize) #want to plot midpoint of bin
            binnedPE.append(binsTotal[indx][1])
            binnedPts.append(binsTotal[indx][2])
    
    # generate radial distribution function 
    finalG = radialDistFun(binSize,binnedR,binnedPts, currTimeStep)

    # plot the RDF 
    VMDFile = open('RDF.AGR','r')
    data = np.loadtxt(VMDFile)
    VMDFile.close()
    VMDRad = data[:,0]
    gr = data[:,1]
    rdfPlot = plt.figure(1)
    plt.plot(binnedR,finalG,'b') 
    plt.plot(VMDRad,gr,'r')
    plt.xlim(0,5)
    plt.ylim(0, 3)
    plt.xlabel("Distance (d)")
    plt.ylabel("g(r)")
    plt.legend(["Calculated", "VMD"])
    plt.show()

        
    #plot potential energy at each time step
    pePlot = plt.figure(2)    
    plt.plot(binnedR,binnedPE,'b^',ms = 3) #points
    plt.plot(binnedR,binnedPE, 'b', lw = 1) #line
    lammpsFile = open('potentialOuput.data','r')
    PEdata = np.loadtxt(lammpsFile)
    lammpsFile.close()
    allRad = PEdata[:,1]
    PE = PEdata[:,2]
    plt.plot(allRad,PE,'ro',ms = 3) #points
    plt.plot(allRad,PE,'r',lw = 1) #line
    plt.xlabel("Distance (d)")
    plt.ylabel("Potential Energy (kT)")
    plt.legend(["Calculated", "Calculated Trajectory", "LAMMPS", "LAMMPS Trajectory"])
    plt.xlim(0.5,4)
    plt.ylim(-1.5,1.5)
    plt.show()

#run functions    
main()


