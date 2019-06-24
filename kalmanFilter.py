#Imports
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import csv

#Create some global, one-time variables
dt = .0145
A = np.matrix([[1,0,0,-dt,0,0],[0,1,0,0,-dt,0],[0,0,1,0,0,-dt],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
B = np.matrix([[dt,0,0],[0,dt,0],[0,0,dt],[0,0,0],[0,0,0],[0,0,0]])
C = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0]])
P = np.identity(6)
Q = np.identity(6) * .0008
R = np.matrix([[.1,0,0],[0,.1,0],[0,0,10]])
stateEstimate = np.matrix([[0],[0],[0],[0],[0],[0]])

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
	
	p = row[1]
	q = row[2]
	r = row[3]
	
	a_x = row[4]
	a_y = row[5]
	a_z = row[6]
	accNorm = math.sqrt(a_x^2 + a_y^2 + a_z^2)
	a_x = a_x/accNorm
	a_y = a_y/accNorm
	a_z = a_z/accNorm
	
	m_x = row[7]
	m_y = row[8]
	m_z = row[9]
	magNorm = math.sqrt(m_x^2 + m_y^2 + m_z^2)
	m_x = m_x/magNorm
	m_y = m_y/magNorm
	m_z = m_z/magNorm
	
	
