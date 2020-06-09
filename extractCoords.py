#! /Users/Program Files (x86)/PythonDownload

coordFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/allCoords.txt','r')
outFile = open('/Users/gbonn/Summer_Research_2020/lammps_tut/testOutput.txt','w')
for lineNum, line in enumerate(coordFile):
    if(lineNum == 50100):
        break
    else:
        outFile.write(line)

coordFile.close()
outFile.close()
