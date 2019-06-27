#Imports
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

#Load File
mat = scipy.io.loadmat('IMU_slow.mat')

#Print file
#print mat["gyroReadings"]
gyroReadingsMatrix = np.matrix(mat['gyroReadings'])
accelReadingsMatrix = np.matrix(mat['accelReadings'])
magReadingsMatrix = np.matrix(mat['magReadings'])
trueAngleMatrix = np.matrix(mat['Euler_ang_true'])

gyro = []
true = []
times = []
for i in range(0,1000):
	print (accelReadingsMatrix.item((i,1)))
