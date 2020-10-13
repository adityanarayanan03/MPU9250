'''
NOTES:
1.	The Euler Angle Calculation functions are being
	called multiple times with the same values. If
	an efficient code is desired, pass the 
	Kalman Filter calculated values right off the bat.
	
2.	The magnetometer calibration done in getData() must
	be tailored to the environment, and re-done every time
	an axis change is desired.
'''
#Imports
import FaBo9Axis_MPU9250
import time
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import RPi.GPIO as GPIO
import threading

#Creating Objects
sensor = FaBo9Axis_MPU9250.MPU9250()

#Setting up the GPIO ports for stepper motor
GPIO.setwarnings(False)
GPIO.setmode(GPIO.BOARD)

pulsePin = 8
directionPin = 10
enablePin = 12

GPIO.setup(pulsePin, GPIO.OUT)
GPIO.setup(directionPin, GPIO.OUT)
GPIO.setup(enablePin, GPIO.OUT)

GPIO.output(enablePin, False)
GPIO.output(directionPin, False) #Switch to true for other way

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
global gyroEulerAngle
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

accXCorrection = 0
accYCorrection = 0
accZCorrection = 0

#Time-Related Variables:
global dt
dt = 0.0208 #delta T for integration
global times
times = [] #x-axis for plotting
runTime = 0#user desired run-time

#Variables for Filtering
A = np.matrix([[1,0,0,-dt,0,0],[0,1,0,0,-dt,0],[0,0,1,0,0,-dt],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
B = np.matrix([[dt,0,0],[0,dt,0],[0,0,dt],[0,0,0],[0,0,0],[0,0,0]])
C = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0]])
P = np.identity(6)
Q = np.identity(6) * .0005
R = np.matrix([[.1,0,0],[0,.1,0],[0,0,10]]) * 100
state_estimate = np.matrix([[0],[0],[0],[0],[0],[0]])
filteredData = np.matrix([[0],[0],[0]])
filteredRoll = []
filteredPitch = []
filteredYaw = []

#Variables for true data collection
global true_angles
global true_times
global end

end = "dont end"
true_angles = []
true_times = []
	
#Creating directories if user wants to save data
saveIndicator = 'y'
if (saveIndicator== "y" or saveIndicator == "Y"):
	fileName = str(input("Enter File Name:   "))
	os.mkdir(fileName)
	os.mkdir(fileName+'/rawOutput')
	os.mkdir(fileName+'/angleEstimations')
	os.mkdir(fileName+'/comparisons')
	os.mkdir(fileName+'/comparisons/filteredComparisons')
	
#Functions for calculating Euler Angles	
def accCalcEuler(sensorValues):
	phi = math.atan2((sensorValues[5]),(math.sqrt((sensorValues[4])**2 + (sensorValues[6])**2)))
	phi = math.degrees(phi)
	theta = math.atan2((-1*sensorValues[4]),(math.sqrt((sensorValues[5])**2 + (sensorValues[6])**2)))
	theta = math.degrees(theta)
	psi = 0
	eulerEstimate = np.matrix([[phi],[theta],[psi]]) 
	return eulerEstimate

def gyroCalcEuler(sensorValues, prevEulerAngle, dt_local):
	currentAngularVelocities = np.matrix([[math.radians(sensorValues[1])],[math.radians(sensorValues[2])],[math.radians(sensorValues[3])]])
	phi = prevEulerAngle.item((0,0))
	theta = prevEulerAngle.item((1,0))
	M = np.matrix([[1, (math.tan(theta)*math.sin(phi)), (math.tan(theta)*math.cos(phi))], [0, math.cos(phi), -1*math.sin(phi)], [0, (math.sin(phi)/math.cos(theta)), (math.cos(phi)/math.cos(theta))]])
	eulerAngle = prevEulerAngle + M*currentAngularVelocities*dt_local
	return eulerAngle

def magCalcEuler(sensorValues, accEulerAngle, magCorrectionMean, gyroEstimate):
	phi = 0
	theta = 0
	accPhi = math.radians(accEulerAngle.item((0,0)))
	accTheta = math.radians(accEulerAngle.item((1,0)))
	psi = math.atan2((sensorValues[9]*math.sin(accPhi) - sensorValues[8]*math.cos(accPhi)),((sensorValues[7]*math.cos(accTheta))+(sensorValues[8]*math.sin(accTheta)*math.sin(accPhi))+(sensorValues[9]*math.sin(accTheta)*math.cos(accPhi))))
	if (magCorrectionMean == 0):
		psi = math.degrees(psi)
	else:
		psi = (math.degrees(psi) - magCorrectionMean)
	eulerEstimate2 = np.matrix([[phi],[theta],[psi]])
	return eulerEstimate2

#Function to collect data and store in temporary list
def getData(currentTime, gyroRollCorrection, gyroPitchCorrection, gyroYawCorrection):
	#read data from sensor
	accel = sensor.readAccel()				
	gyro = sensor.readGyro()				
	mag = sensor.readMagnet()				
	
	#get Time into a list
	times.append(currentTime)
	
	#get gyro values into proper list w/ corrections
	gyroXCorr = gyro['x'] - gyroRollCorrection
	gyroRoll.append(gyroXCorr)
	gyroYCorr = gyro['y'] - gyroPitchCorrection
	gyroPitch.append(gyroYCorr)
	gyroZCorr = gyro['z'] - gyroYawCorrection
	gyroYaw.append(gyroZCorr)
	
	#Correct and append accelerometer values
	accelXCorr = accel['x'] - accXCorrection
	accX.append(accelXCorr)
	accelYCorr = accel['y'] - accYCorrection
	accY.append(accelYCorr)
	accelZCorr = accel['z'] + (9.8005-9.78941028858)
	accZ.append(accelZCorr)
	
	#Sphere method for calibration of magnetometer:
	#Correction Matrices
	#A_sphere = np.matrix([[1.0, 0.0, .0020],[-.0055, 1.0287, .0136],[.0020, .0136, .9459]])
	b_sphere = np.matrix([[281.8552,-1.3048,-23.3093]])
	rawMagVals = np.matrix([[mag['x'],mag['y'],mag['z']]])
	correctedMag = np.matrix([[0,0,0]])
	
	#calculate
	correctedMag = (rawMagVals - b_sphere)
	
	#Append
	magXCorr = correctedMag.item((0,0))
	magX.append(magXCorr)
	magYCorr = correctedMag.item((0,1))
	magY.append(magYCorr)
	magZCorr = correctedMag.item((0,2))
	magZ.append(magZCorr)
	
	#append corrected values to a list and that list to data set
	sensorValues = [currentTime, gyroXCorr, gyroYCorr, gyroZCorr, accelXCorr, accelYCorr, accelZCorr, magXCorr, magYCorr, magZCorr]
	dataSet.append(sensorValues)
	return sensorValues
	
#Function for Kalman Filtering
def kalmanFilter(sensorValues, magCorrectionMean, gyroYaw):
	#tell python to modify global forms of variables
	global A
	global B
	global C
	global P
	global Q
	global R
	global state_estimate
	
	#Collect Gyro Data
	p = math.radians(sensorValues[1])
	q = math.radians(sensorValues[2])
	r = math.radians(sensorValues[3])
	
	#calculate Angles from measurements
	phi_hat_acc = math.radians(accCalcEuler(sensorValues).item((0,0)))
	theta_hat_acc = math.radians(accCalcEuler(sensorValues).item((1,0)))
		
	phi_hat = state_estimate.item((0,0))
	theta_hat = state_estimate.item((1,0))
	psi_hat = state_estimate.item((2,0))
	
	psi_hat_mag= math.radians(magCalcEuler(sensorValues,accCalcEuler(sensorValues), magCorrectionMean, gyroYaw).item((2,0)))
	#psi_hat_mag = 2 * psi_hat_mag - math.radians(2 * magCorrectionMean)
	
	phi_dot = p+math.sin(phi_hat)*math.tan(theta_hat)*q+math.cos(phi_hat)*math.tan(theta_hat)*r
	theta_dot = math.cos(phi_hat)*q - math.sin(phi_hat)*r
	psi_dot = math.sin(phi_hat)/math.cos(theta_hat)*q + math.cos(phi_hat)/math.cos(theta_hat)*r
	delta_angle = np.matrix([[phi_dot],[theta_dot],[psi_dot]])
	
	#predict actual attitude
	state_estimate = A * state_estimate + B * delta_angle
	P = A*P*np.transpose(A) + Q
	
	#Update readings
	Z = np.matrix([[phi_hat_acc],[theta_hat_acc],[psi_hat_mag]])
	r = Z - C*state_estimate
	S = R + C*P*np.transpose(C)
	K = P*np.transpose(C)*(np.linalg.inv(S)) 
	state_estimate = state_estimate + K*r
	P = (np.identity(6) - K*C) * P;
	state_estimate_degrees = state_estimate * (180/math.pi)
	return state_estimate_degrees
	
#Template function for graphing and saving three-series of data
def graph(figNum,series1,label1,series2,label2,series3,label3,plotTitle,plotYLabel,savePath):
	plt.figure(figNum)
	plt.plot(times, series1, label = label1)
	plt.plot(times, series2, label = label2)
	plt.plot(times, series3, label = label3)
	plt.legend(loc = 'upper right')
	plt.title(plotTitle)
	plt.xlabel("Time(s)")
	plt.ylabel(plotYLabel)
	plt.savefig("/home/pi/MPU9250/"+fileName+savePath)	
	
#Nine calls to graph function with proper data series
def saveAndPlot():
	global true_times
	global true_angles
	np.savetxt("/home/pi/MPU9250/"+fileName+"/rawData.csv", dataSet, delimiter = ',')
	graph(1,gyroYaw,'gyroYaw',gyroPitch,'gyroPitch',gyroRoll,'gyroRoll',"Gyro Output vs. Time","Gyro Output (dps)","/rawOutput/gyroOutput.png")
	graph(2,accX,'accX',accY,'Acc Y',accZ,'Acc Z',"Acc. Output vs. Time","Acc. Output (m/s/s)","/rawOutput/accelerometerOutput.png")
	graph(3,magX,'Mag X',magY,'Mag Y',magZ,'Mag Z',"Magnetometer Output vs. Time","Magnetometer Output (uT)","/rawOutput/magnetometerOutput.png")
	graph(4,gyroEulerYaw,'gyroEulerYaw',gyroEulerPitch,'gyroEulerPitch',gyroEulerRoll,'gyroEulerRoll',"Gyro Calculated Angles","Angle (degrees)","/angleEstimations/gyroEulerAngle.png")
	graph(5,accEulerYaw,'accEulerYaw',accEulerPitch,'accEulerPitch',accEulerRoll,'accEulerRoll',"Accelerometer Calculated Angles","Angle (degrees)","/angleEstimations/accEulerAngle.png")
	#graph(6,magEulerYaw,'magEulerYaw',magEulerPitch,'magEulerPitch',magEulerRoll,'magEulerRoll',"Magnetometer Calculated Angles","Angle (degrees)","/angleEstimations/magEulerAngle.png")
	#graph(7,gyroEulerYaw,'gyroEulerYaw',accEulerYaw,'accEulerYaw',magEulerYaw,'magEulerYaw',"Yaw Angle Comparisons","Angle (degrees)","/comparisons/yawComparisons.png")
	#graph(8,gyroEulerPitch,'gyroEulerPitch',accEulerPitch,'accEulerPitch',magEulerPitch,'magEulerPitch',"Pitch Angle Comparisons","Angle (degrees)","/comparisons/pitchComparisons.png")
	graph(9,gyroEulerRoll,'gyroEulerRoll',accEulerRoll,'accEulerRoll',magEulerRoll,'gmagEulerRoll',"Roll Angle Comparisons","Angle (degrees)","/comparisons/rollComparisons.png")
	graph(10,accEulerRoll,'accEulerRoll',filteredRoll,'filteredRoll',gyroEulerRoll,'gyroRoll','Roll Angle Filter Comparison','Angle (degrees)','/comparisons/filteredComparisons/roll.png')
	#graph(11,accEulerPitch,'accEulerPitch',filteredPitch,'filteredPitch',gyroEulerPitch,'gyroPitch','Pitch Angle Filter Comparison','Angle (degrees)','/comparisons/filteredComparisons/pitch.png')
	#graph(12,magEulerYaw,'magEulerYaw',filteredYaw,'filteredYaw',gyroEulerYaw,'gyroEulerYaw','Yaw Angles Filter Comparison','Angle (degrees)','/comparisons/filteredComparisons/yaw.png')
	plt.figure(13)
	plt.plot(true_times, true_angles, label = 'True Value')
	plt.plot(times, gyroEulerRoll, label = 'Gyro Estimate')
	plt.plot(times, accEulerRoll, label = 'Accelerometer Estimate')
	plt.plot(times, filteredRoll, label = 'Filtered Estimate')
	plt.legend(loc = 'upper right')
	plt.title('Roll estimates Compared to True Value')
	plt.xlabel('Time (sec)')
	plt.ylabel('Angle (degrees)')
	plt.savefig('/home/pi/MPU9250/'+fileName+'/comparisons/rollEstimatesAndTrue.png')
	
	
#Calibration Sequence
print 'calibrating'
for i in range(100):
	calSensorValues = getData(0,0,0,0)
	calAccEuler = accCalcEuler(calSensorValues)
	calMagEuler = magCalcEuler(calSensorValues, calAccEuler, 0, -50)
	magPsiValues.append(calMagEuler.item((2,0)))
	time.sleep(.009)
magCorrectionMean = np.mean(magPsiValues)
print magCorrectionMean
gyroRollCorrection = np.mean(gyroRoll)
gyroPitchCorrection = np.mean(gyroPitch)
gyroYawCorrection = np.mean(gyroYaw)
accXCorrection = np.mean(accX)
accYCorrection = np.mean(accY)

print 'calibration over, starting in 3 seconds'
time.sleep(3)
print 'starting'

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
global startTime
startTime = time.time()
currentTime = 0
errorVal = []

#Stepper Motor Sequence
def stepper_motor():
	global end
	global true_times
	global true_angles
	
	for i in range (25000/6):
		GPIO.output(pulsePin, True)
		time.sleep(.00025)
		GPIO.output(pulsePin, False)
		time.sleep(.00025)
		true_angles.append((360.0*i)/25000.0)
		true_times.append(time.time() - startTime)
        
	GPIO.output(directionPin, True)
    
	for r in range (25000/3):
		GPIO.output(pulsePin, True)
		time.sleep(.00025)
		GPIO.output(pulsePin, False)
		time.sleep(.00025)
		true_angles.append(60 - (360.0*r)/25000.0)
		true_times.append(time.time() - startTime)
        
	GPIO.cleanup()
	end = "end"

def stepper_100():
	global end
	end = "wait"
	time.sleep(100)
	end = "end"

#data collection sequence
def data_collect():
	global gyroEulerAngle
	global true_times
	global true_angles
	dt_gyro = 0.0
	
	while(end!="end"):
		currentTime = time.time() - startTime
		sensorValues = getData(currentTime, gyroRollCorrection, gyroPitchCorrection, gyroYawCorrection)
		print currentTime
		
		if (len(times) >= 2): #does a dt exist yet?
			dt_gyro = times[-1] - times[-2]
		
		#Call to calculations and storing calculated Values for Plotting
		gyroEulerAngle = gyroCalcEuler(sensorValues, gyroEulerAngle, dt_gyro)
		gyroEulerRoll.append(math.degrees(gyroEulerAngle.item((0,0))))
		gyroEulerPitch.append(math.degrees(gyroEulerAngle.item((1,0))))
		gyroEulerYaw.append(math.degrees(gyroEulerAngle.item((2,0))))

		accEulerAngle = accCalcEuler(sensorValues)
		accEulerRoll.append(accEulerAngle.item((0,0)))
		accEulerPitch.append(accEulerAngle.item((1,0)))
		accEulerYaw.append(accEulerAngle.item((2,0)))
		
		magEulerAngle = magCalcEuler(sensorValues, accEulerAngle, magCorrectionMean, gyroEulerAngle.item((2,0)))
		magEulerRoll.append(magEulerAngle.item((0,0)))
		magEulerPitch.append(magEulerAngle.item((1,0)))
		magEulerYaw.append(magEulerAngle.item((2,0)))
			
		#Call to filtering functions
		filteredData = kalmanFilter(sensorValues, magCorrectionMean, math.degrees(gyroEulerAngle.item((2,0))))
		filteredRoll.append(filteredData.item((0,0)))
		filteredPitch.append(filteredData.item((1,0)))
		filteredYaw.append(filteredData.item((2,0)))
		time.sleep(.009)
		
	#save Plots if requested
	if(saveIndicator == "y" or saveIndicator == "Y"):
		print "Saving data. This takes forever."
		saveAndPlot()

#Threading
threading.Thread(target = stepper_100).start()
threading.Thread(target = data_collect).start()
