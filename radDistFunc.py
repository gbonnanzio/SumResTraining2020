import numpy as np
import matplotlib.pyplot as plt
import math
def radialDistFun(binSize):
    radialFreqFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/radialDistFun.txt','r')
    data = np.loadtxt(radialFreqFile)
    allRad = data[:,0]
    allCts = data[:,1]
    #print(len(allRad))
    #print(len(allCts))
    #binSize = (5.0*math.sqrt(10))/50.0	
    radDistVal = []
    n = []
    nID = []
    for i in range(len(allRad)):
        n.append(allCts[i]/(500.0))
        nID.append(((4.0*math.pi)/3.0)*(math.pow((allRad[i]+binSize),3)-math.pow(allRad[i],3)))
        radDistVal.append(n[i]/nID[i])

    plt.plot(allRad,radDistVal)
    plt.xlim[0,4]
    plt.show()

