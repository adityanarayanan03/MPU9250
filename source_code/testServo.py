import RPi.GPIO as GPIO
import time

GPIO.setmode(GPIO.BOARD)

GPIO.setup(16, GPIO.OUT)

p = GPIO.PWM(16, 50)

p.start(2.5)
print 'started'
time.sleep(3)
a = 2.5
try:
    while (a <= (2.5+1.6666)):
        p.ChangeDutyCycle(a)
        a+= .005
        time.sleep(.05)
except KeyboardInterrupt:
    p.stop()
    GPIO.cleanup()
