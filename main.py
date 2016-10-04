import os
import numpy as np


def variation(a=10, param=0, component=0):
    start = np.loadtxt('Start.txt')

    start[param, component] += a

    coordinates = '  ' + str(start[0, 0]) + ' ' + str(start[0, 1]) + ' ' + str(start[0, 2])
    velocities = '  ' + str(start[1, 0]) + ' ' + str(start[1, 1]) + ' ' + str(start[1, 2])

    text = []
    f = open('iszm_puc.in', 'r')
    for line in f:
        text.append(line)

    f.close()

    text[3] = coordinates + text[3][-2] + text[3][-1]
    text[4] = velocities + text[4][-2] + text[4][-1]

    f = open('iszm_puc.in', 'w')
    for line in text:
        f.write(line)

    f.close()

    os.system('wine iszm_puc.exe')

    eph = open('EPH.OUT', 'r')
    coordinates = []
    count = 1
    obs = 0
    for line in eph:
        # if count == 2:
        #     nulCor = line
        # if count == 3:
        #     nulVel = line
        if count == 10:
            coordinates.append(line)
            obs += 1
        if count == 8 * obs + 10:
            coordinates.append(line)
            obs += 1
        count += 1

    eph.close()

    # startCor = np.array([float(nulCor[7:29]), float(nulCor[33:55]), float(nulCor[60:82])])
    # startVel = np.array([float(nulVel[7:29]), float(nulVel[33:55]), float(nulVel[60:82])])

    new = []
    for count in range(0, 100):
        new.append(float(coordinates[count][7:29]))
        new.append(float(coordinates[count][33:55]))
        new.append(float(coordinates[count][60:82]))

    # coordinates = None
    # nulCor = None
    # nulVel = None

    return new

a = 10 ** -5
dictionary = {'coordinate': 0,
              'velocity': 1,
              'x': 0,
              'y': 1,
              'z': 2}
model = np.loadtxt('Observed.txt')
start = np.loadtxt('Start.txt')
iteration = 1
sigma = 10.0
while True:
    obs = np.array(variation(0))
    # pro = np.array(variation(a, 0, 0))
    A = np.zeros([300, 6])

    for count in range(0, 3):
        A[:, count] = (np.array(variation(a, 0, count)) - obs) / a

    for count in range(0, 3):
        A[:, count + 3] = (np.array(variation(a, 1, count)) - obs) / a

    c = obs - model

    Q = np.dot(A.T, A)
    d = np.dot(A.T, c)
    y = np.dot(np.linalg.inv(Q), d)

    sigma0 = np.sqrt(np.dot(c, c) / (300 - 6))
    sigma_old = sigma
    sigma = np.array([np.linalg.inv(Q)[0, 0], np.linalg.inv(Q)[1, 1], np.linalg.inv(Q)[2, 2]])
    sigma *= sigma0

    start[0] -= y[0:3]
    start[1] -= y[3:6]

    np.savetxt('Start.txt', start)
    print iteration
    iteration += 1
    eph = open('EPH.OUT', 'r')
    eph.readline()
    print eph.readline()
    print eph.readline()
    eph.close()
    print np.sqrt(np.dot(sigma, sigma))
    if abs(np.sqrt(np.dot(sigma, sigma)) - np.sqrt(np.dot(sigma_old, sigma_old))) < 0.000001:
        break
