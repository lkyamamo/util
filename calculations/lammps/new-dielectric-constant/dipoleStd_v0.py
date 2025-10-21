import sys

file1 = sys.argv[1]
T = float(sys.argv[2])
la = float(sys.argv[3])
lb = float(sys.argv[4])
lc = float(sys.argv[5])

mavg = [0.0,0.0,0.0]
msq = 0.0
counter = 0.0
minitial = [0.0,0.0,0.0]

mxx, myy, mzz = 0, 0, 0

with open(file1,"r") as f1:
	for line in f1.readlines():

		line = line.strip().split()
		
		m1,m2,m3 = float(line[0]), float(line[1]), float(line[2])
		if counter == 0:
			minitial[0], minitial[1], minitial[2] = m1, m2, m3
		#m1,m2,m2 = m1 - minitial[0], m2 - minitial[1], m3 - minitial[2]
		msq += (m1*m1 + m2*m2 + m3*m3)
		mxx += m1* m1
		myy += m2*m2
		mzz += m3*m3

		mavg[0] += m1
		mavg[1] += m2
		mavg[2] += m3

		counter += 1.0


msq /= float(counter)

mxx /= float(counter)
myy /= float(counter)
mzz /= float(counter)

mavg[0] /= float(counter)
mavg[1] /= float(counter)
mavg[2] /= float(counter)

print("%14.6f %14.6f %14.6f" %(mavg[0], mavg[1],mavg[2]))
"""
with open(file1,"r") as f1:
	for line in f1.readlines():
		line = line.strip().split()
		m1,m2,m3 = float(line[0]), float(line[1]), float(line[2])
		m1,m2,m3 = m1-mavg[0],m2-mavg[1],m3-mavg[2]
		msq += (m1*m1) + (m2*m2) + (m3*m3)
"""

prefactor = (2.56*(10**7))/(1.38*8.85*T*la*lb*lc)

mdev = msq - (mavg[0]*mavg[0] + mavg[1]*mavg[1] + mavg[2]*mavg[2])
mx_dev = mxx - mavg[0]**2
my_dev = myy - mavg[1]**2
mz_dev = mzz - mavg[2]**2
mdev = mdev / 3
print("x_dev = %12.6f, y_dev = %12.6f, z_dev = %12.6f,deviation = %14.6f" %(mx_dev, my_dev, mz_dev, mdev))
print("eps_x = %12.6f, eps_y = %12.6f, eps_z = %12.6f, eps_total = %12.6f" %(prefactor * mx_dev, prefactor * my_dev, prefactor * mz_dev, prefactor * mdev))
