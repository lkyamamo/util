import sys
import numpy as np


def LoadXYZ(xyz):
        atoms = {}
        bonds = {}
        with open(xyz,'r') as file1:


                lines = file1.readline()        #1
                lines = file1.readline()        #2
                lines = file1.readline()        #3
                natoms = int(lines.strip().split()[0])
                lines = file1.readline()        #4
                lines = file1.readline()        #5
                la = float(lines.strip().split()[1]) - float(lines.strip().split()[0])
                lines = file1.readline()        #6
                lb = float(lines.strip().split()[1]) - float(lines.strip().split()[0])
                lines = file1.readline()        #7
                lc = float(lines.strip().split()[1]) - float(lines.strip().split()[0])
                lines = file1.readline()        #8

                #counter = 0
                for lines in file1:
                        #line = lines[i]
                        line = lines.strip().split()
                        #print(line)
                        atomID = int(line[0])
                        x,y,z,q = float(line[3]), float(line[4]), float(line[5]), float(line[6])
                        bondID = int(line[1])
                        itype = int(line[2])
                        if bondID in bonds:
                                bonds[bondID].append(atomID)
                        else:
                                bonds[bondID] = [atomID]

                        #print(atomID, bondID, itype, x, y, z, q)

                        atoms[atomID] = [itype,x,y,z,q]
                        #print("%s   %14.6f  %14.6f   %14.6f   %14.6f   %d" %(atype,x,y,z,q,atomID))
        return natoms,atoms,bonds,la,lb,lc

def TotalMoment(atoms, bonds, la, lb, lc):
        totalMoment = [0.0,0.0,0.0]

        dipole = []

        counter = 0
        for iid in bonds.keys():
                counter += 1

                atomIDs = sorted(bonds[iid])

                atomID1, atomID2, atomID3 = atomIDs[0], atomIDs[1], atomIDs[2]


                dx12 = atoms[atomID1][1] - atoms[atomID2][1]
                if dx12 >= la/2.0:
                        dx12 -= la
                elif dx12 <= -la/2.0:
                        dx12 += la


                dx31 = atoms[atomID3][1] - atoms[atomID1][1]
                if dx31 >= la/2.0:
                        dx31 -= la
                elif dx31 <= -la/2.0:
                        dx31 += la


                dy12 = atoms[atomID1][2] - atoms[atomID2][2]
                if dy12 >= lb/2.0:
                        dy12 -= lb
                elif dy12 <= -lb/2.0:
                        dy12 += lb


                dy31 = atoms[atomID3][2] - atoms[atomID1][2]
                if dy31 >= lb/2.0:
                        dy31 -= lb
                elif dy31 <= -lb/2.0:
                        dy31 += lb


                dz12 = atoms[atomID1][3] - atoms[atomID2][3]
                if dz12 >= lc/2.0:
                        dz12 -= lc
                elif dz12 <= -lc/2.0:
                        dz12 += lc


                dz31 = atoms[atomID3][3] - atoms[atomID1][3]
                if dz31 >= lc/2.0:
                        dz31 -= lc
                elif dz31 <= -lc/2.0:
                        dz31 += lc

                dipoleMoment = [0,0,0]

                q = 0.32983	
	
                dipoleMoment[0] += (dx12)*q
                dipoleMoment[0] -= (dx31)*q
 
                dipoleMoment[1] += (dy12)*q
                dipoleMoment[1] -= (dy31)*q

                dipoleMoment[2] += (dz12)*q
                dipoleMoment[2] -= (dz31)*q

                dipole_mag = np.linalg.norm(np.array(dipoleMoment))
                dipole.append(dipole_mag)

                #print(dipoleMoment[0], dipoleMoment[1], dipoleMoment[2])

                totalMoment[0] += dipoleMoment[0]
                totalMoment[1] += dipoleMoment[1]
                totalMoment[2] += dipoleMoment[2]


        print('---------------------')
        print(np.mean(np.array(dipole)))


        return totalMoment

start = int(sys.argv[1])
end = int(sys.argv[2])
incr = int(sys.argv[3])
counter = 0

file10 = open(sys.argv[4],"w")

for i in range(start,end,incr):
        xyz = "298.%d" %i
        moment = []

        natoms, atoms, bonds, la, lb, lc = LoadXYZ(xyz)

        moment = TotalMoment(atoms,bonds, la, lb, lc)

        counter += 1
        print(la, lb, lc)
        print("%d %14.6f  %14.6f  %14.6f" %(counter,moment[0],moment[1],moment[2]))
        file10.write("%14.6f  %14.6f  %14.6f\n" %(moment[0],moment[1],moment[2]))


file10.close()

