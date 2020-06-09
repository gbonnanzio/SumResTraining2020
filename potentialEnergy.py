#! Program Files (x86)/PythonDownload/python
import matplotlib.pyplot as plt
import numpy as np
import math
import os

def determinePE(coordFile):
    #coordFile is a txt file of the coords at a single time step
    data = np.loadtxt(fname = coordFile)

    xPos = data[:,2]
    yPos = data[:,3]
    zPos = data[:,4]

    nAtoms = data.shape[0]
    totalPotential = 0

    for i in range(nAtoms):
        for j in range(i+1,nAtoms):
            if i != j:
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

                #constants
                Ar_diam = 142.0 #nm
                #from White 1999
                eps = 125.0*(8.6173*10**-5) #eV
                sigma = 0.3345 #nm

                rCorrected = r*sigma
                sigOverR = sigma/rCorrected
                tempPot = 4.0*eps*((sigOverR**12)-(sigOverR**6))

                #add to total potential
                totalPotential = totalPotential + tempPot


    return totalPotential
    

def main():
    #specify where text files are located
    filePath = '/Users/gbonn/Summer Research 2020/lammps_tut/coordTime'
    #using 6 text files
    textFiles = ['0000000.txt','0200000.txt','0400000.txt','0600000.txt','0800000.txt','1000000.txt']
    timeStep = []
    potentialEnergy = []
    #calculate and keep track of potential energy for each time step
    for fileName in textFiles:
        tempFile = filePath + fileName
        tempPE = determinePE(tempFile)
        timeStep.append(float(fileName[:7]))
        potentialEnergy.append(tempPE)
        #stringPE = str(determinePE(tempFile)) 
        #print("At time step " + fileName[:7] + " potential energy equals " + stringPE + " eV.")
    
    #plot potential energy at each time step
    plt.plot(timeStep,potentialEnergy)
    plt.xlabel('Time Step')
    plt.ylabel('Potential Energy (eV)')
    plt.show()

#run functions    
main()


