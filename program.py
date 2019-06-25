#Imports
import FaBo9Axis_MPU9250
import time
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import os

#Creating Object
sensor = FaBo9Axis_MPU9250.MPU9250()

#Empty list to hold full set of raw data
dataSet = []

#Empty lists to hold raw data
gyroYaw = []
gyroPitch = []
gyroRoll = []

accX = []
accY = []
accZ = []

magX = []
magY = []
magZ = []

#Empty lists to hold Calculated Values
gyroEulerAngle = np.matrix([[0],[0],[0]])
gyroEulerYaw = []
gyroEulerPitch = []
gyroEulerRoll = []

accEulerAngle = np.matrix([[0],[0],[0]])
accEulerYaw = []
accEulerPitch = []
accEulerRoll = []

magEulerAngle = np.matrix([[0],[0],[0]])
magEulerYaw = []
magEulerPitch = []
magEulerRoll = []

#Variables for Calibration Sequence
calAccEuler = np.matrix([[0],[0],[0]])
calMagEuler = np.matrix([[0],[0],[0]])
magPsiValues = []

calGyroYaw = []
calGyroPitch = []
calGyroRoll = []

gyroRollCorrection = 0
gyroPitchCorrection = 0
gyroYawCorrection = 0

#Time-Related Variables:
dt = .0145 #delta T for integration
times = [] #x-axis for plotting
runTime = int(input("Time for trial:   ")) #user desired run-time

#Creating directories if user wants to save data
saveIndicator = str(input("Do you want to save data from this run? y/n    "))
if (saveIndicator== "y" or saveIndicator == "Y"):
	fileName = str(input("Enter File Name:   "))
	os.mkdir(fileName)
	os.mkdir(fileName+'/rawOutput')
	os.mkdir(fileName+'/angleEstimations')
	os.mkdir(fileName+'/comparisons')
	
#Functions for calculating Euler Angles	
def accCalcEuler(sensorValues):
    phi = math.atan((sensorValues[5])/(math.sqrt((sensorValues[4])**2 + (sensorValues[6])**2)))
    phi = math.degrees(phi)
    theta = math.atan((-1*sensorValues[4])/(math.sqrt((sensorValues[5])**2 + (sensorValues[6])**2)))
    theta = math.degrees(theta)
    psi = 0
    eulerEstimate = np.matrix([[-1*phi],[-1*theta],[-1*psi]]) 
    return eulerEstimate    
    
def gyroCalcEuler(sensorValues, prevEulerAngle):
    currentAngularVelocities = np.matrix([[math.radians(sensorValues[1])],[math.radians(sensorValues[2])],[math.radians(sensorValues[3])]])
    phi = prevEulerAngle.item((0,0))
    theta = prevEulerAngle.item((1,0))
    M = np.matrix([[1, (math.tan(theta)*math.sin(phi)), (math.tan(theta)*math.cos(phi))], [0, math.cos(phi), -1*math.sin(phi)], [0, (math.sin(phi)/math.cos(theta)), (math.cos(phi)/math.cos(theta))]])
    eulerAngle = prevEulerAngle + M*currentAngularVelocities*dt
    return eulerAngle
    
def magCalcEuler(sensorValues, accEulerAngle):
    phi = 0
    theta = 0
    accPhi = -1 * math.radians(accEulerAngle.item((0,0)))
    accTheta = -1 * math.radians(accEulerAngle.item((1,0)))
    psi = math.atan2((sensorValues[9]*math.sin(accPhi) - sensorValues[8]*math.cos(accPhi)),((sensorValues[7]*math.cos(accTheta))+(sensorValues[8]*math.sin(accTheta)*math.sin(accPhi))+(sensorValues[9]*math.sin(accTheta)*math.cos(accPhi))))
    psi = math.degrees(psi)
    eulerEstimate = np.matrix([[phi],[theta],[psi]])
    return eulerEstimate
    
#Function to collect data and store in temporary list    
def getData(currentTime):
	accel = sensor.readAccel()             #Reads the acceleration list? from the sensor
	gyro = sensor.readGyro()               #Reads the gyro list? from the sensor
	mag = sensor.readMagnet()              #Reads the magnetometer list? from the sensor
	times.append(currentTime)
	gyroXCorr = -1 * gyro['x'] - gyroRollCorrection
	gyroRoll.append(gyroXCorr)
	gyroYCorr = gyro['y'] - gyroPitchCorrection
	gyroPitch.append(gyroYCorr)
	gyroZCorr = -1 * gyro['z'] -gyroYawCorrection
	gyroYaw.append(gyroZCorr)
	accelXCorr = -1 * (10 * accel['x'] - .0660602258469)
	accX.append(accelXCorr)
	accelYCorr = 10 * accel['y'] - .0519196988708
	accY.append(accelYCorr)
	accelZCorr = 10 * accel['z'] + (9.8005-9.78941028858)
	accZ.append(accelZCorr)
	magXCorr = mag['x']
	magX.append(magXCorr)
	magYCorr = mag['y']
	magY.append(magYCorr)
	magZCorr = mag['z']
	magZ.append(magZCorr)
	sensorValues = [currentTime, gyroXCorr, gyroYCorr, gyroZCorr, accelXCorr, accelYCorr, accelZCorr, magXCorr, magYCorr, magZCorr]
	dataSet.append(sensorValues)
	return sensorValues
	
#Template function for graphing and saving three-series of data
def graph(figNum,series1,label1,series2,label2,series3,label3,plotTitle,plotYLabel,savePath):
	plt.figure(figNum)
	plt.plot(times, series1, label = label1)
	plt.plot(times, series2, label = label2)
	plt.plot(times, series3, label = label3)
	plt.legend(loc = 'upper right')
	plt.title(plotTitle)
	plt.xlabel("Time(ms)")
	plt.ylabel(plotYLabel)
	plt.savefig("/home/pi/MPU9250/"+fileName+savePath)	
	
#Nine calls to graph function with proper data series
def saveAndPlot():
	np.savetxt("/home/pi/MPU9250/"+fileName+"/rawData.csv", dataSet, delimiter = ',')
	graph(1,gyroYaw,'gyroYaw',gyroPitch,'gyroPitch',gyroRoll,'gyroRoll',"Gyro Output vs. Time","Gyro Output (dps)","/rawOutput/gyroOutput.png")
	graph(2,accX,'accX',accY,'Acc Y',accZ,'Acc Z',"Acc. Output vs. Time","Acc. Output (m/s/s)","/rawOutput/accelerometerOutput.png")
	graph(3,magX,'Mag X',magY,'Mag Y',magZ,'Mag Z',"Magnetometer Output vs. Time","Magnetometer Output (uT)","/rawOutput/magnetometerOutput.png")
	graph(4,gyroEulerYaw,'gyroEulerYaw',gyroEulerPitch,'gyroEulerPitch',gyroEulerRoll,'gyroEulerRoll',"Gyro Calculated Angles","Angle (degrees)","/angleEstimations/gyroEulerAngle.png")
	graph(5,accEulerYaw,'accEulerYaw',accEulerPitch,'accEulerPitch',accEulerRoll,'accEulerRoll',"Accelerometer Calculated Angles","Angle (degrees)","/angleEstimations/accEulerAngle.png")
	graph(6,magEulerYaw,'magEulerYaw',magEulerPitch,'magEulerPitch',magEulerRoll,'magEulerRoll',"Magnetometer Calculated Angles","Angle (degrees)","/angleEstimations/magEulerAngle.png")
	graph(7,gyroEulerYaw,'gyroEulerYaw',accEulerYaw,'accEulerPitch',magEulerYaw,'magEulerRoll',"Yaw Angle Comparisons","Angle (degrees)","/comparisons/yawComparisons.png")
	graph(8,gyroEulerPitch,'gyroEulerPitch',accEulerPitch,'accEulerPitch',magEulerPitch,'magEulerPitch',"Pitch Angle Comparisons","Angle (degrees)","/comparisons/pitchComparisons.png")
	graph(9,gyroEulerRoll,'gyroEulerRoll',accEulerRoll,'accEulerRoll',magEulerRoll,'gmagEulerRoll',"Roll Angle Comparisons","Angle (degrees)","/comparisons/rollComparisons.png")
	
#Calibration Sequence
print 'calibrating'
for i in range(100):
	calSensorValues = getData(0)
	calAccEuler = accCalcEuler(calSensorValues)
	calMagEuler = magCalcEuler(calSensorValues, calAccEuler)
	calGyroRoll.append(calSensorValues[1])
	calGyroPitch.append(calSensorValues[2])
	calGyroYaw.append(calSensorValues[3])
	magPsiValues.append(calMagEuler.item((2,0)))
	time.sleep(.007)
magCorrectionMean = np.mean(magPsiValues)
print magCorrectionMean
gyroRollCorrection = np.mean(calGyroRoll)
gyroPitchCorrection = np.mean(calGyroPitch)
gyroYawCorrection = np.mean(calGyroYaw)
print 'calibration over, starting in 3 seconds'
time.sleep(3)

#Resetting Variables Common to calibration sequence and data-collection
times = []
dataSet = []

gyroYaw = []
gyroPitch = []
gyroRoll = []

accX = []
accY = []
accZ = []

magX = []
magY = []
magZ = []	

#Starting Count for data-collection
startTime = time.time()
currentTime = 0

#data collection sequence
while (currentTime<runTime):
	currentTime = time.time() - startTime
	#print currentTime
	sensorValues = getData(currentTime)
	
	#Call to calculations and storing calculated Values for Plotting
	gyroEulerAngle = gyroCalcEuler(sensorValues, gyroEulerAngle)
	gyroEulerRoll.append(math.degrees(gyroEulerAngle.item((0,0))))
	gyroEulerPitch.append(math.degrees(gyroEulerAngle.item((1,0))))
	gyroEulerYaw.append(math.degrees(gyroEulerAngle.item((2,0))))
                    
                    
	accEulerAngle = accCalcEuler(sensorValues)
	accEulerRoll.append(accEulerAngle.item((0,0)))
	accEulerPitch.append(accEulerAngle.item((1,0)))
	accEulerYaw.append(accEulerAngle.item((2,0)))
        
	magEulerAngle = magCalcEuler(sensorValues, accEulerAngle)
	magEulerRoll.append(magEulerAngle.item((0,0)))
	magEulerPitch.append(magEulerAngle.item((1,0)))
	magEulerYaw.append(2 * (magEulerAngle.item((2,0)) - magCorrectionMean))
	time.sleep(.007)

#Saving plots if requested
if(saveIndicator == "y" or saveIndicator == "Y"):
	print "Saving data. This takes forever."
	saveAndPlot()

print magCorrectionMean
