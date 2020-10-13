#Imports
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import csv
import scipy.io

#Create some global variables
dt = .010
A = np.matrix([[1,0,0,-dt,0,0],[0,1,0,0,-dt,0],[0,0,1,0,0,-dt],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
B = np.matrix([[dt,0,0],[0,dt,0],[0,0,dt],[0,0,0],[0,0,0],[0,0,0]])
C = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0]])
P = np.identity(6)
Q = np.identity(6) * .0008
R = np.matrix([[.1,0,0],[0,.1,0],[0,0,10]])
state_estimate = np.transpose([[0,0,0,0,0,0]])

phi_angles = []
theta_angles = []
psi_angles = []

true_phi = []
true_theta= []
true_psi = []

times = []

roll_error = []
pitch_error = []
yaw_error = []

#get, open, and read the csv file
mat = scipy.io.loadmat('IMU_slow.mat')
gyroReadingsMatrix = np.matrix(mat['gyroReadings'])
accelReadingsMatrix = np.matrix(mat['accelReadings'])
magReadingsMatrix = np.matrix(mat['magReadings'])
trueAngleMatrix = np.matrix(mat['Euler_ang_true'])

#iterate through dataSet
for i in range(0,1000):
	#get list of values into individual values
	currentTime = i
	times.append(currentTime)
	
	#read gyro values
	p = float(gyroReadingsMatrix.item((i,0))) - .349
	q = (float(gyroReadingsMatrix.item((i,1))) - .349)
	r = float(gyroReadingsMatrix.item((i,2))) - .349
	
	#read and normalize acc Values
	a_x = float(accelReadingsMatrix.item((i,0)))
	a_y = float(accelReadingsMatrix.item((i,1)))
	a_z = float(accelReadingsMatrix.item((i,2)))
	accNorm = 1.00#math.sqrt(a_x**2 + a_y**2 + a_z**2)
	a_x = a_x/accNorm
	a_y = a_y/accNorm
	a_z = a_z/accNorm
	
	#read and normalize mag values
	m_x = float(magReadingsMatrix.item((i,0)))
	m_y = float(magReadingsMatrix.item((i,1)))
	m_z = float(magReadingsMatrix.item((i,2)))
	magNorm = 1.00#math.sqrt(m_x**2 + m_y**2 + m_z**2)
	x_mag = m_x/magNorm
	y_mag = m_y/magNorm
	z_mag = m_z/magNorm
	
	#save true values in right places
	true_phi.append(trueAngleMatrix.item((i,2)))
	true_theta.append(trueAngleMatrix.item((i,1)))
	true_psi.append(trueAngleMatrix.item((i,0)))
	
	
	#Get Angles from raw measurements
	phi_hat_acc = math.atan2((a_y),math.sqrt(a_x**2 + a_z**2))
	theta_hat_acc = math.atan2((a_x),math.sqrt(a_y**2+a_z**2))
	
	phi_hat = state_estimate.item((0,0))
	theta_hat = state_estimate.item((1,0))
	psi_hat = state_estimate.item((2,0))
	
	psi_hat_mag=math.atan2((-1*y_mag*math.cos(phi_hat)+z_mag*math.sin(phi_hat)),(x_mag*math.cos(theta_hat)+y_mag*math.sin(theta_hat)*math.sin(phi_hat)+z_mag*math.sin(theta_hat)*math.cos(phi_hat)))
	psi_hat_mag -= .0873
	
	phi_dot = p+math.sin(phi_hat)*math.tan(theta_hat)*q+math.cos(phi_hat)*math.tan(theta_hat)*r
	theta_dot = math.cos(phi_hat)*q - math.sin(phi_hat)*r
	psi_dot = math.sin(phi_hat)/math.cos(theta_hat)*q + math.cos(phi_hat)/math.cos(theta_hat)*r
	delta_angle = np.matrix([[phi_dot],[theta_dot],[psi_dot]])
	
	#predict attitude
	state_estimate = A * state_estimate + B * delta_angle
	P = A*P*np.transpose(A) + Q
	
	#update
	Z = np.matrix([[phi_hat_acc],[theta_hat_acc],[psi_hat_mag]])
	r = Z - C*state_estimate
	S = R + C*P*np.transpose(C)
	K = P*np.transpose(C)*(np.linalg.inv(S)) 
	state_estimate = state_estimate + K*r
	P = (np.identity(6) - (K*C)) * P
	state_estimate_degrees = state_estimate * (180/math.pi)
	phi_angles.append(state_estimate_degrees.item((0,0)))
	theta_angles.append(-1*state_estimate_degrees.item((1,0)))
	psi_angles.append(state_estimate_degrees.item((2,0)))
	roll_error.append(state_estimate_degrees.item((0,0)) - trueAngleMatrix.item((i,2)))
	pitch_error.append(-1 * state_estimate_degrees.item((1,0)) - trueAngleMatrix.item((i,1)))
	yaw_error.append(state_estimate_degrees.item((2,0)) - trueAngleMatrix.item((i,0)))	
	
plt.figure(1)
plt.plot(times, roll_error, label = 'Phi error')
#plt.plot(times, phi_angles, label = 'Filtered Phi')
plt.legend(loc = 'upper right')
plt.title('True phi vs Filtered')
plt.xlabel('Time (readings)')
plt.ylabel('Angle (degrees)')
plt.savefig('rollError')
plt.show()

plt.figure(2)
plt.plot(times, pitch_error, label = 'Theta Error')
#plt.plot(times, theta_angles, label = 'Filtered Theta')
plt.legend(loc = 'upper right')
plt.title('True theta vs Filtered')
plt.xlabel('Time (readings)')
plt.ylabel('Angle (degrees)')
plt.savefig('pitchError')
plt.show()

plt.figure(3)
plt.plot(times, yaw_error, label = 'Psi Error')
#plt.plot(times, psi_angles, label = 'Filtered Psi')
plt.legend(loc = 'upper right')
plt.title('True psi vs Filtered')
plt.xlabel('Time (readings)')
plt.ylabel('Angle (degrees)')
plt.savefig('yawError')
plt.show()
