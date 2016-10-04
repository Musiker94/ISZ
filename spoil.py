import numpy as np
import math as m
import random


eph = open('EPH(20).OUT', 'r')
coordinates = []
count = 1
obs = 0
for line in eph:
    if count == 2:
        nulCor = line
    if count == 3:
        nulVel = line
    if count == 10:
        coordinates.append(line)
        obs += 1
    if count == 8 * obs + 10:
        coordinates.append(line)
        obs += 1
    count += 1

eph.close()

startCor = np.array([float(nulCor[7:29]), float(nulCor[33:55]), float(nulCor[60:82])])
startVel = np.array([float(nulVel[7:29]), float(nulVel[33:55]), float(nulVel[60:82])])

new = []
for count in range(0, 100):
    new.append(float(coordinates[count][7:29]))
    new.append(float(coordinates[count][33:55]))
    new.append(float(coordinates[count][60:82]))

coordinates = None
nulCor = None
nulVel = None

dx = []
for count in range(0, 298, 3):
    dx.append(m.sqrt(new[count] ** 2 + new[count + 1] ** 2 + new[count + 2] ** 2))

new = np.array(new)
dx = np.array(dx)
for count in range(0, 100):
    new[3 * count:3 * count+3] += dx[count] * random.uniform(-1.0, 1.0) * m.radians(1.0 / 360000.0)

dx0 = np.sqrt(np.dot(startCor, startCor)) * random.uniform(-1.0, 1.0) * m.radians(1.0 / 360000.0) * 100
dv0 = np.sqrt(np.dot(startVel, startVel)) * random.uniform(-1.0, 1.0) * m.radians(1.0 / 360000.0) * 100
startCor += dx0
startVel += dv0
start = np.array([startCor, startVel])

# np.savetxt('Observed(20).txt', new)
# np.savetxt('Start(20).txt', start)
