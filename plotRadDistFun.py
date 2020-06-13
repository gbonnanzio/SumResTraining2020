import numpy as np
import matplotlib.pyplot as plt
import math

def GrVIM():
    file = open('/Users/gbonn/Summer_Research_2020/lammps_tut/RDF.AGR','r')
    data = np.loadtxt(radialFreqFile)
    allRad = data[:,0]
    gr = data[:,1]
    plt.plot(allRad,gr)
    plt.xlim[0,4]
    plt.show()
    plt.xlabel("Distance (D)")
    lt.ylabel("g(r)")
    plt.legend(["g(r)"])
    
def main():
    GrVIM()
    
main()

