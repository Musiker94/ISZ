__author__ = 'KMS 2015'

import math as m
import numpy as np
import random


def mu():
    #mu = 2.9591220828559110225E-4  # Gravitational constant = [au^3/day^2]
    mu = 3.9864418E05
    return mu

#MU = 2.9591220828559110225E-4  # Gravitational constant = [au^3/day^2]
MU = mu()


def date(JD):
    a = int(round(JD)) + 32044
    b = (4 * a + 3) / 146097
    c = a - 146097 * b / 4
    d = (4 * c + 3) / 1461
    e = c - 1461 * d / 4
    m = (5 * e + 2) / 153
    day = e - (153 * m + 2) / e +1
    month = m + 3 - 12 * (m / 10)
    year = 100 * b + d - 4800 + m / 10
    date3 = np.array([day, month, year])
    return date3


def julian(day, month, year, hour, minute, second):
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    JDN = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    JD = float(JDN) + (float(hour) - 12.0) / 24.0 + float(minute) / 1440.0 + float(second) / 86400.0
    return JD


def py_angle(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cos = np.dot(v1, v2)
    sin = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sin, cos)


def mean(e, v):
    E = 2 * m.atan(m.sqrt((1 - e) / (1 + e)) * m.tan(v / 2))
    M = E - e * m.sin(E)
    if M < 0:
        M += 2 * m.pi
    elif M >= 2 * m.pi:
        M -= 2 * m.pi
    return M


def cross(orbit10):
    file_nodes_down = open('nodes_down.txt', 'a')
    file_nodes_up = open('nodes_up.txt', 'a')

    orbit10[0] = orbit10[0] / (1 - orbit10[1])
    count = 0
    while count != 4:
        orbit10[count + 2] = m.radians(orbit10[count + 2])
        count += 1

    rdown = orbit10[0] * (1 - m.pow(orbit10[1], 2)) / (1 + orbit10[1] * m.cos(3 * m.pi - orbit10[4]))
    nodedown2 = np.array([-rdown * m.cos(orbit10[3]), -rdown * m.sin(orbit10[3])])
    rup = orbit10[0] * (1 - m.pow(orbit10[1], 2)) / (1 + orbit10[1] * m.cos(2 * m.pi - orbit10[4]))
    nodeup2 = np.array([rup * m.cos(orbit10[3]), rup * m.sin(orbit10[3])])
    file_nodes_down.write(str(nodedown2[0]) + '\t' + str(nodedown2[1]) + '\t' + str(orbit10[9]) + '\n')
    file_nodes_up.write(str(nodeup2[0]) + '\t' + str(nodeup2[1]) + '\t' + str(orbit10[9]) + '\n')
    file_nodes_down.close()
    file_nodes_up.close()
    return 1


def kepler2descart(kep6):
    u = kep6[5] + kep6[4]
    par = kep6[0] * (1 - m.pow(kep6[1], 2))
    par_div_r = 1 + kep6[1] * m.cos(kep6[5])
    r = par / par_div_r
    n = m.sqrt(MU / par)
    vr = n * kep6[1] * m.sin(kep6[5])
    vn = n * par_div_r
    #v_sqr = MU * (2 * kep6[0] - r) / (r * kep6[0])
    p = []
    q = []
    xx6 = []

    cos_u = m.cos(u)
    sin_u = m.sin(u)
    cos_om = m.cos(kep6[3])
    sin_om = m.sin(kep6[3])
    cos_i = m.cos(kep6[2])
    sin_i = m.sin(kep6[2])

    p.append(cos_u * cos_om - sin_u * sin_om * cos_i)
    p.append(cos_u * sin_om + sin_u * cos_om * cos_i)
    p.append(sin_u * sin_i)
    q.append(-sin_u * cos_om - cos_u * sin_om * cos_i)
    q.append(-sin_u * sin_om + cos_u * cos_om * cos_i)
    q.append(cos_u * sin_i)
    count = 0
    while count != 3:
        xx6.append(r * p[count])
        count += 1

    count = 0
    while count != 3:
        xx6.append(p[count] * vr + q[count] * vn)
        count += 1
    return xx6


def descart2kepler(x6):
    x3 = np.array([x6[0], x6[1], x6[2]])
    vel3 = np.array([x6[3], x6[4], x6[5]])
    r = m.sqrt(np.dot(x3, x3))
    sigma3 = np.cross(x3, vel3)
    sigma = m.sqrt(np.dot(sigma3, sigma3))
    cos_i = sigma3[2] / sigma
    sin_i = m.sqrt(1 - m.pow(cos_i, 2))
    i = np.arctan2(sin_i, cos_i)
    sin_om = sigma3[0] / (sigma * sin_i)
    cos_om = - sigma3[1] / (sigma * sin_i)
    om = np.arctan2(sin_om, cos_om)
    par = m.pow(sigma, 2) / MU
    h = np.dot(vel3, vel3) - 2 * MU / r
    e = 1 + h * m.pow(sigma, 2) / m.pow(MU, 2)		#kvadrat ex
    a = par / (1 - e)
    e = m.sqrt(e)
    tg_v = np.dot(x3, vel3) * m.sqrt(par / MU) / (par - r)
    cos_v = (par-r)/(r*e)
    sin_v = tg_v * cos_v
    v = np.arctan2(sin_v, cos_v)
    tg_u = x3[2] / (sin_i * (x3[0] * cos_om + x3[1] * sin_om))
    sin_u = x3[2] / (r * sin_i)
    cos_u = sin_u / tg_u
    u = np.arctan2(sin_u, cos_u)
    if v != 0.0:
        w = u - v + 2 * m.pi
    else:
        w = u

    if w >= 2 * m.pi:
        w -= 2 * m.pi

    if om < 0:
        om += 2 * m.pi
    elif om >= 2 * m.pi:
        om -= 2 * m.pi

    if v < 0:
        v += 2 * m.pi
    elif v >= 2 * m.pi:
        v -= 2 * m.pi

    kep6 = np.array([a, e, i, om, w, v])
    return kep6


# date = 11.12.2015 7:00
def ejectoid(orbit10, nV=10, nOrb=1):

    file_particles = open('particles.txt', 'w')

    q0 = orbit10[0]
    e0 = orbit10[1]
    i0 = m.radians(orbit10[2])
    om0 = m.radians(orbit10[3])
    w0 = m.radians(orbit10[4])
    v0 = m.radians(orbit10[5])
    year = orbit10[6]
    JD0 = orbit10[7]
    number = orbit10[8]
    veject = orbit10[9]

    rc = 1
    dens=1.9
    mass=7.6E-06

    rpart = m.pow(0.75 * mass / (m.pi * dens), 1.0 / 3.0)  # radius = [cm]

    a0 = q0 / (1 - e0)

    p = a0 * (1 - m.pow(e0, 2))
    Q = a0 * (1 + e0)                       # aphelion

    kep60 = np.array([a0, e0, i0, om0, w0, v0])

    countV = 1
    while countV != (nV + 1):
        r = m.pow((1 / (random.random() * (m.pow(Q, -3.0) - m.pow(q0, -3.0)) + m.pow(q0, -3.0))), 1.0 / 3.0)
        v0 = (p - r) / (r * e0)
        v0 = m.acos(v0)
        if random.random() < 0.5:
            v0 = 2 * m.pi - v0

        kep60[5] = v0
        x60 = kepler2descart(kep60)
        x30 = ([x60[0], x60[1], x60[2]])
        countOrb = 1
        while countOrb != (nOrb + 1):
            cos_T = 1 - 2 * random.random()
            sin_T = m.pow((1 - m.pow(cos_T, 2)), 2)
            fi = 2 * m.pi * random.random()
            c = m.sqrt(1 / (dens * rpart * m.pow(np.dot(x30, x30), 9.0/8.0)) - 0.013 * rc) * m.sqrt(rc) * 656	# cm/sec
            c = c * (6.68458712E-9 / 1.15740741)	# au/day
            c3 = np.array([c * cos_T, c * sin_T * m.cos(fi), c * sin_T* m.sin(fi)])
            x6 = x60
            count = 0
            while count != 3:
                x6[count + 3] += c3[count]
                count += 1

            kep6 = descart2kepler(x6)

            M = 0.0
            JD = 0.0
            M = mean(kep6[1], kep6[5])
            M0 = 0.0
            n = m.sqrt(MU / m.pow(kep6[0], 3))
            if M >= m.pi:
                deltaM = M0 - M
            else:
                deltaM = M0 + M

            JD = (deltaM + JD0 * n) / n
            date3 = date(JD)

            kep6[0] = kep6[0] * (1 - kep6[1])
            kep6[2] = m.degrees(kep6[2])
            kep6[3] = m.degrees(kep6[3])
            kep6[4] = m.degrees(kep6[4])
            kep6[5] = m.degrees(kep6[5])

            countKep = 0
            while countKep != 6:
                file_particles.write(str(kep6[countKep]) + '\t')
                countKep += 1

            file_particles.write(str(date3[2]) + '\t' + str(JD) + '\t' + str(countV) + '\t' + str(m.degrees(v0)) + '\n')

            countOrb += 1

        countV += 1

    file_particles.close()


# test function date = 11.12.2015 7:00
def eject(q0=1.189, e0=0.619, i0=70.876, w0=171.354,
          om0=282.962, v0=0, JD=2457367.7916666665, rc=1, dens=1.9, mass=7.958E-03, nV=10):

    file_eject = open('eject.txt', 'w')

    i0 = m.radians(i0)
    om0 = m.radians(om0)
    w0 = m.radians(w0)

    rpart = m.pow(0.75 * mass / (m.pi * dens), 1.0 / 3.0)  # radius cm

    a0 = q0 / (1 - e0)

    p = a0 * (1 - m.pow(e0, 2))
    Q = a0 * (1 + e0)                       #apogelii

    kep60 = np.array([a0, e0, i0, om0, w0, v0])

    countV = 1
    while countV != (nV + 1):
        r = m.pow((1 / (random.random() * (m.pow(Q, -3.0) - m.pow(q0, -3.0)) + m.pow(q0, -3.0))), 1.0 / 3.0)
        v0 = (p - r) / (r * e0)
        v0 = m.acos(v0)
        if random.random() < 0.5:
            v0 = 2 * m.pi - v0

        kep60[5] = v0
        x60 = kepler2descart(kep60)
        x30 = ([x60[0], x60[1], x60[2]])
        cos_T = 1 - 2 * random.random()
        sin_T = m.pow((1 - m.pow(cos_T, 2)), 2)
        fi = 2 * m.pi * random.random()
        c = m.sqrt(1 / (dens * rpart * m.pow(np.dot(x30, x30), 9.0/8.0)) - 0.013 * rc) * m.sqrt(rc) * 656	#cm/sec
        c = c * (6.68458712E-9 / 1.15740741)	#au/day
        c3 = np.array([c * cos_T, c * sin_T * m.cos(fi), c * sin_T* m.sin(fi)])
        x6 = x60
        count = 0
        while count != 3:
            x6[count + 3] += c3[count]
            count += 1

        kep6 = descart2kepler(x6)
        M = mean(kep6[1], kep6[5])
        M0 = 0
        n = m.sqrt(MU / m.pow(kep6[0], 3))
        if M >= m.pi:
            deltaM = M0 - M
        else:
            deltaM = M0 + M

        JD = (deltaM + JD * n) / n

        kep6[2] = m.degrees(kep6[2])
        kep6[3] = m.degrees(kep6[3])
        kep6[4] = m.degrees(kep6[4])
        kep6[5] = m.degrees(kep6[5])
        kep6[5] = m.degrees(M)

        count = 0
        while count != 6:
            file_eject.write(str(kep6[count]).zfill(12) + '\t')
            count += 1
        file_eject.write(str(JD).zfill(12) + '\t' + str(countV).zfill(4) + '\t' + str(m.degrees(v0)).zfill(12) + '\n')

        countV += 1

    file_eject.close()