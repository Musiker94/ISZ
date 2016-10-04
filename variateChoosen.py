import os
import numpy as np
import matplotlib.pyplot as plt


def variation(a=10, param=0, component=0):
    start = np.loadtxt('Start(20).txt')

    if param != 2:
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

    elif param == 2:
        f_garm = open(r'MODTOOLS/garm360.in', 'r')
        f_garm.seek(147)
        str_garm = f_garm.readline()
        print(str_garm.split(' '))
        garm = float(str_garm.split(' ')[8]) + a
        f_garm.seek(0)
        str_garm = '   2   0 ' + str(garm) + ' 0.000000000000E+00' + str_garm[-1:]

        text = []
        for line in f_garm:
            text.append(line)

        f_garm.close()
        text[component] = str_garm

        f_garm = open(r'MODTOOLS/garm360.in', 'w')
        for line in text:
            f_garm.write(line)

        f_garm.close()

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


def get_variate(param='coordinate', component='x', count_start=10.0**-10.0, stop=10000.0):
    dictionary = {'coordinate': 0,
                  'velocity': 1,
                  'garm': 2,
                  'x': 0,
                  'y': 1,
                  'z': 2,
                  '20': 3}
    variate = []
    while count_start != stop * 10.0:
        variate.append(variation(count_start, dictionary[param], dictionary[component]))
        print count_start
        count_start *= 10

    variate = np.array(variate)

    return variate

parr1 = 'garm'
parr = '20'
var = get_variate(parr1, parr)
obs = np.array(variation(0))

a = 10.0 ** -10.0
for count in range(0, var.shape[0]):
    print a
    var[count] = (var[count] - obs) / a
    a *= 10

ss = 'vari' + parr1 + parr.upper() + '.txt'
np.savetxt(ss, var)
