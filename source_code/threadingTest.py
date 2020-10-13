import threading
import time
import stepperControl as SC
import matplotlib.pyplot as plt

def loop_3():
	global end
	global angleTime
	end = "no"
	angleTime = SC.motor_control()
	end = "end"

def loop_4():
	while (end != "end"):
		print "threading"
	plt.plot(angleTime[1], angleTime[0])
	plt.show()
threading.Thread(target=loop_3).start()
threading.Thread(target=loop_4).start()
