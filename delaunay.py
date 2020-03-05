#!/usr/bin/env python
# -*- coding=utf-8 -*-
'''
This part is for finding delaunay reduced cell
'''

import numpy as np
import os.path


def DelauExchBasis(b, uni_ax):
    if uni_ax == -1:
        flag = 1
        for k in range(0, 3):
            for i in range(k+1, 4):
                scal_prod = 0.0
                for j in range(0, 3):
                    scal_prod += b[k, j] * b[i, j]
                #print(scal_prod)
                if scal_prod > 0.00001:
                    for l in range(4):
                        if not (l == i or l == k):
                            for j in range(3):
                                b[l, j] += b[k, j]
                    for j in range(3):
                        b[k, j] = -b[k, j]
                    flag = 0
                    return flag, b
        return flag, b

    elif uni_ax != -1:
        flag = 1
        for k in range(0, 2):
            for i in range(k+1, 3):
                scal_prod = 0.0
                for j in range(0, 3):
                    scal_prod += b[k, j] * b[i, j]
                #print(scal_prod)
                if scal_prod > 0.00001:
                    for l in range(3):
                        if not (l == i or l == k):
                            for j in range(3):
                                b[l, j] += b[k, j]
                    for j in range(3):
                        b[k, j] = -b[k, j]
                    flag = 0
                    return flag, b
        return flag, b


def Delaunay(latt, uni_ax):
    if uni_ax == -1:
        flag = None
        # delaunay minimum sum
        b = np.zeros((4, 3))
        for i in range(3):
            for j in range(3):
                b[i, j] = latt[i, j]
        for i in range(3):
            b[3, i] = -latt[0, i] - latt[1, i] - latt[2, i]

        for count in range(100):
            #print('count: ', count)
            flag, b = DelauExchBasis(b, uni_ax)
            if flag:
                #print('Succeed!')
                break

        # search the 3 shortest lattice vec as basis vec
        v = np.zeros((7, 3))
        for i in range(0, 4):
            v[i] = b[i]
        v[4] = b[0] + b[1]
        v[5] = b[1] + b[2]
        v[6] = b[2] + b[0]

        #reduc_set = [(v[i], "%.10f" % np.linalg.norm(v[i])) for i in range(0, 7)]
        reduc_set = [(v[i], np.around(np.linalg.norm(v[i]), decimals=5)) for i in range(0, 7)]
        sortreduc = sorted(reduc_set, key=lambda setx: setx[1])
        #print('reduc_set:\n', reduc_set)
        #print('sortreduc:\n', sortreduc)
        reduc_b = np.zeros((3, 3))

        countp = 0
        vol_tolerance = 0.00001
        for i in range(2, 7):
            reduc_b[0] = sortreduc[0][0]
            reduc_b[1] = sortreduc[1][0]
            reduc_b[2] = sortreduc[i][0]
            if abs(np.linalg.det(reduc_b)) > vol_tolerance:
                break

        # determine whether the volume based on the shortest vec is appropriate or not
        volume = np.linalg.det(reduc_b)
        #print('check reb1', reduc_b)
        #print('old vol:\n', volume)
        if abs(volume) < vol_tolerance:
            print("Error: No Volume!")
            return 0, None, None
        if volume < 0.0:
            reduc_b = -reduc_b.copy()
        #print('check reb2', reduc_b)
        #print('new volume:\n', np.linalg.det(reduc_b))

        # determinate of delaunay transformation
        inv_latt = np.linalg.inv(latt)
        delauP = np.transpose(np.dot(reduc_b, inv_latt))   # // P is vertical vec
        #print('delaunayP1:\n', delauP)

        for i in range(3):
            for j in range(3):
                if delauP[i, j] < 0.0:
                    delauP[i, j] = int(delauP[i, j] - 0.5)
                else:
                    delauP[i, j] = int(delauP[i, j] + 0.5)
        #print('delaunayP2:\n', delauP)
        #print('det P:\n', np.linalg.det(delauP))

        if abs(np.linalg.det(delauP)) != 1:
            print("Error: The abs of determinant of the delaunay matrix should be 1")
            return 0, None, None
        #print(reduc_b)

        return flag, reduc_b, delauP

    elif uni_ax != -1:
        flag = None
        # delaunay minimum sum
        uni_vec = latt[uni_ax]

        b = np.zeros((3, 3))
        ct = 0
        for i in range(3):
            if i != uni_ax:
                b[ct] = latt[i]
                ct += 1

        for i in range(3):
            b[2, i] = -latt[0, i] - latt[1, i]

        for count in range(100):
            flag, b = DelauExchBasis(b, uni_ax)
            if flag:
                break
            #else:
                #print(count, b)

        # search the 3 shortest lattice vec as basis vec
        v = np.zeros((4, 3))
        for i in range(0, 3):
            v[i] = b[i]
        v[3] = b[0] + b[1]

        #reduc_set = [(v[i], "%.10f" % np.linalg.norm(v[i])) for i in range(0, 4)]
        reduc_set = [(v[i], np.around(np.linalg.norm(v[i]), decimals=5)) for i in range(0, 4)]
        sortreduc = sorted(reduc_set, key=lambda setx: setx[1])
        #print('reduc_set:\n', reduc_set)
        #print('sortreduc:\n', sortreduc)

        tp_reduc_b = np.zeros((3, 3))
        tp_reduc_b[0] = sortreduc[0][0]
        tp_reduc_b[1] = uni_vec
        #countp = 0
        vol_tolerance = 0.00001
        for i in range(1, 4):
            tp_reduc_b[2] = sortreduc[i][0]
            if abs(np.linalg.det(tp_reduc_b)) > vol_tolerance:
                break
            #else:
                #countp += 1
        # print('countp:\n', countp)
        # print('tp_reduc_b:\n', tp_reduc_b)
        reduc_b = np.zeros((3, 3))
        k = 0
        for i in range(3):
            if i != uni_ax:
                reduc_b[i] = tp_reduc_b[k]
                k += 2
            else:
                reduc_b[i] = uni_vec

        # determine whether the volume based on the shortest vec is appropriate or not
        volume = np.linalg.det(reduc_b)
        #print('old vol:\n', volume)
        if abs(volume) < vol_tolerance:
            print("Error: No Volume!")
            return 0, None, None
        if volume < 0.0:
            reduc_b[uni_ax] = -reduc_b[uni_ax]
        #volume = np.linalg.det (reduc_b)
        #print('new volume:\n', np.linalg.det(reduc_b))
        # determinate of delaunay transformation

        inv_latt = np.linalg.inv(latt)
        delauP = np.transpose(np.dot(reduc_b, inv_latt))  # // 3.13 this Q actually equals P:  (a'b'c') = (a b c) P   vertical vec

        #print('delaunayP1:\n', delauP)
        for i in range(3):
            for j in range(3):
                if delauP[i, j] < 0.0:
                    delauP[i, j] = int(delauP[i, j] - 0.5)
                else:
                    delauP[i, j] = int(delauP[i, j] + 0.5)
        #print('delaunayP2:\n', delauP)

        #print('det P:\n', np.linalg.det(delauP))

        if abs(np.linalg.det(delauP)) != 1:
            print("Error: The abs of determinant of the delaunay matrix should be 1")
            return 0, None, None
        #print(reduc_b)

        return flag, reduc_b, delauP


def ChangeOfBasis(coor, Q):     # // Q should be vertical vec
    #print('trans_matQ', Q)
    tpos = []
    for i in coor:
        tp = np.dot(np.array(i), np.transpose(Q))
        tp2 = tp.copy()     # Watch out! if you let "tp2 = tp", changing tp2 will result in tp changed, too
        for j in range(3):
            if tp2[j] < 0.0:
                tp2[j] = int(tp2[j] - 0.5)
            else:
                tp2[j] = int(tp2[j] + 0.5)
        tp = tp - tp2
        #print('3:\n', tp)
        for j in range(3):
            if tp[j] < 0.0:
                tp[j] += 1.0
        #print('3n:\n', tp)
        tpos.append(tp)
    return tpos

'''/// test
latt, pos, num, dictp = StructRead()
a, b, p = Delaunay(latt, -1)
print('delaunay:\n', b)
print('delaunay pos\n', ChangeOfBasis(pos, np.linalg.inv(p)))
c = spglib.delaunay_reduce(latt)
print('spgresult:\n', c)
print('volume spg:\n', abs(np.linalg.det(c)))
d, x2, x3 = spglib.find_primitive((latt, pos, num))
print(d)

/// test '''