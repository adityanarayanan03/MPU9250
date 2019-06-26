#Imports
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import csv

#Create some global variables
dt = .0145
A = np.matrix([[1,0,0,-dt,0,0],[0,1,0,0,-dt,0],[0,0,1,0,0,-dt],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
B = np.matrix([[dt,0,0],[0,dt,0],[0,0,dt],[0,0,0],[0,0,0],[0,0,0]])
C = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0]])
P = np.identity(6)
Q = np.identity(6) * .0008
R = np.matrix([[.1,0,0],[0,.1,0],[0,0,10]])
state_estimate = np.matrix([[0],[0],[0],[0],[0],[0]])

phi_angles = []
theta_angles = []
psi_angles = []

times = []

#get, open, and read the csv file
fileName = str(input(":              "))
dataSet = open("/home/pi/MPU9250/"+fileName+"/rawData.csv")
fullData = csv.reader(dataSet)

#iterate through dataSet
for row in fullData:
	#convert all string readings to float values
	for i in range(0,10):
		row[i] = float(row[i])
		
	#get list of values into individual values
	currentTime = row[0]
	times.append(currentTime)
	
	#read gyro values
	p = math.radians(row[1])
	q = math.radians(row[2])
	r = math.radians(row[3])
	
	#read and normalize acc Values
	a_x = row[4]
	a_y = row[5]
	a_z = row[6]
	accNorm = 1#math.sqrt(a_x**2 + a_y**2 + a_z**2)
	a_x = a_x/accNorm
	a_y = a_y/accNorm
	a_z = a_z/accNorm
	
	#read and normalize mag values
	m_x = row[7]
	m_y = row[8]
	m_z = row[9]
	magNorm = 1#math.sqrt(m_x**2 + m_y**2 + m_z**2)
	x_mag = m_x/magNorm
	y_mag = m_y/magNorm
	z_mag = m_z/magNorm
	
	#Get Angles from raw measurements
	phi_hat_acc = math.atan(a_y/ math.sqrt(a_x**2 + a_z**2))
	theta_hat_acc = math.atan((a_x) / math.sqrt(a_y**2+a_z**2))
	
	phi_hat = state_estimate.item((0,0))
	theta_hat = state_estimate.item((1,0))
	psi_hat = state_estimate.item((2,0))
	
	psi_hat_mag=math.atan((-1*y_mag*math.cos(phi_hat)+z_mag*math.sin(phi_hat))/(x_mag*math.cos(theta_hat)+y_mag*math.sin(theta_hat)*math.sin(phi_hat)+z_mag*math.sin(theta_hat)*math.cos(phi_hat)))
	psi_hat_mag = 2 * psi_hat_mag - math.radians(-15.1)
	
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
	P = (np.identity(6) - K*C) * P;
	state_estimate_degrees = state_estimate * (180/math.pi)
	phi_angles.append(-1 * state_estimate_degrees.item((0,0)))
	theta_angles.append(state_estimate_degrees.item((1,0)))
	psi_angles.append(state_estimate_degrees.item((2,0)))	
	
plt.figure(1)
plt.plot(times, psi_angles, label = 'yaw Axis')
plt.plot(times, theta_angles, label = 'pitch Axis')
plt.plot(times, phi_angles, label = 'roll Axis')
plt.legend(loc = 'upper right')
plt.savefig('/home/pi/MPU9250/'+fileName+'/angleEstimations/inverseFunction.png')
