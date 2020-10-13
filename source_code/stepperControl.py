def motor_control():
    import time
    import RPi.GPIO as GPIO

    pulsePin = 8
    directionPin = 10
    enablePin = 12


    GPIO.setmode(GPIO.BOARD)
    GPIO.setup(pulsePin, GPIO.OUT)
    GPIO.setup(directionPin, GPIO.OUT)
    GPIO.setup(enablePin, GPIO.OUT)

    GPIO.output(enablePin, False)
    GPIO.output(directionPin, False)
    angles = []
    times = []
    for i in range (25000/6):
        GPIO.output(pulsePin, True)
        time.sleep(.00025)
        GPIO.output(pulsePin, False)
        time.sleep(.00025)
        angles.append((360.0*i)/25000.0)
        times.append(i)
    GPIO.output(directionPin, True)
    for r in range (25000/3):
        GPIO.output(pulsePin, True)
        time.sleep(.00025)
        GPIO.output(pulsePin, False)
        time.sleep(.00025)
        angles.append(60 - (360.0*r)/25000.0)
        times.append(25000/6 + r)
    GPIO.cleanup()
    return [angles, times]
    #return times
