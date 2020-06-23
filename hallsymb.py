#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Created on 2019-03-18
Update  on 2020-06-22
Author: Cecilia9999
GitHub: https://github.com/Cecilia9999/
"""

'''
    This file is for decompsing hall symbol
'''

import numpy as np
import csv
import spglib


E = np.identity(3)

latt_sym = {
    'P': [[0, 0, 0]],
    'A': [[0, 0, 0], [0, 1. / 2, 1. / 2]],
    'B': [[0, 0, 0], [1. / 2, 0, 1. / 2]],
    'C': [[0, 0, 0], [1. / 2, 1. / 2, 0]],
    'I': [[0, 0, 0], [1. / 2, 1. / 2, 1. / 2]],
    'R': [[0, 0, 0], [2. / 3, 1. / 3, 1. / 3], [1. / 3, 2. / 3, 2. / 3]],
    'H': [[0, 0, 0], [2. / 3, 1. / 3, 0], [1. / 3, 2. / 3, 0]],
    'F': [[0, 0, 0], [0, 1. / 2, 1. / 2], [1. / 2, 0, 1. / 2], [1. / 2, 1. / 2, 0]]
}

trans_sym = {
    'a': [ 1./2, 0, 0 ],
    'b': [ 0, 1./2, 0 ],
    'c': [ 0, 0, 1./2 ],
    'n': [ 1./2, 1./2, 1./2 ],
    'u': [ 1./4, 0, 0 ],
    'v': [ 0, 1./4, 0 ],
    'w': [ 0, 0, 1./4 ],
    'd': [ 1./4, 1./4, 1./4 ]
}

rot_mat = {
    '1x': [ [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ] ],
    '1y': [ [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ] ],
    '1z': [ [ 1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0, 1 ] ],
    '2x': [ [ 1, 0, 0 ],
            [ 0,-1, 0 ],
            [ 0, 0,-1 ] ],
    '2y': [ [-1, 0, 0 ],
            [ 0, 1, 0 ],
            [ 0, 0,-1 ] ],
    '2z': [ [-1, 0, 0 ],
            [ 0,-1, 0 ],
            [ 0, 0, 1 ] ],
    '3x': [ [ 1, 0, 0 ],
            [ 0, 0,-1 ],
            [ 0, 1,-1 ] ],
    '3y': [ [-1, 0, 1 ],
            [ 0, 1, 0 ],
            [-1, 0, 0 ] ],
    '3z': [ [ 0,-1, 0 ],
            [ 1,-1, 0 ],
            [ 0, 0, 1 ] ],
    '4x': [ [ 1, 0, 0 ],
            [ 0, 0,-1 ],
            [ 0, 1, 0 ] ],
    '4y': [ [ 0, 0, 1 ],
            [ 0, 1, 0 ],
            [-1, 0, 0 ] ],
    '4z': [ [ 0,-1, 0 ],
            [ 1, 0, 0 ],
            [ 0, 0, 1 ] ],
    '6x': [ [ 1, 0, 0 ],
            [ 0, 1,-1 ],
            [ 0, 1, 0 ] ],
    '6y': [ [ 0, 0, 1 ],
            [ 0, 1, 0 ],
            [-1, 0, 1 ] ],
    '6z': [ [ 1,-1, 0 ],
            [ 1, 0, 0 ],
            [ 0, 0, 1 ] ],
    '2px': [ [-1, 0, 0 ],     # b-c
             [ 0, 0,-1 ],
             [ 0,-1, 0 ] ],
    '2ppx': [ [-1, 0, 0 ],    # b+c
              [ 0, 0, 1 ],
              [ 0, 1, 0 ] ],
    '2py': [ [ 0, 0,-1 ],     # a-c
             [ 0,-1, 0 ],
             [-1, 0, 0 ] ],
    '2ppy': [ [ 0, 0, 1 ],    # a+c
              [ 0,-1, 0 ],
              [ 1, 0, 0 ] ],
    '2pz': [ [ 0,-1, 0 ],     # a-b
             [-1, 0, 0 ],
             [ 0, 0,-1 ] ],
    '2ppz': [ [ 0, 1, 0 ],    # a+b
              [ 1, 0, 0 ],
              [ 0, 0,-1 ] ],
    '3*': [ [ 0, 0, 1 ],     # a+b+c
            [ 1, 0, 0 ],
            [ 0, 1, 0 ] ]
    }


# hall_symbols = loadcsvfile()
def loadcsvfile():
    with open("../supports/spg.csv", "rt") as CF:
        cr = csv.reader(CF)
        halls = [[col[4], col[6]] for col in cr]
        #print(halls)
        CF.close()
    return halls
    
    
def LookHallSymbol(hallnum):
    hall_symbols = loadcsvfile()
    sgnum = hall_symbols[hallnum-1][0]
    LNATV = hall_symbols[hallnum-1][1]
    oper = []
    for j in LNATV.split():
        oper.append(j)
        # fulloper.append(oper)
        # print(fulloper)
    # print(oper)

    V = None
    if len(oper) > 3:
        if oper[-3].startswith('('):
            V1 = float(oper[-3].strip('('))
            V2 = float(oper[-2])
            V3 = float(oper[-1].strip(')'))
            V = np.array([V1, V2, V3])
            # print(V)
            oper = oper[:-3]
        # else:
            # print("no V")
            # pass
    # else:
        # print('pass')
        # pass

    L = oper[0]
    improR = 1
    N = oper[1:]
    if N[0].startswith('-'):
        improR = -1
        N[0] = N[0].strip('-')
    strN = N[0]
    preR = strN[0]

    T = []
    R1 = []
    T1 = []
    
    if len(strN) == 1:
        Axis = 'z'
        R1.append(improR * np.array(rot_mat[strN[0] + Axis]))   # default z
        R = improR * np.array(rot_mat[strN[0] + Axis])
    
    elif strN[1] in ['x', 'y', 'z']:
        R1.append(improR * np.array(rot_mat[strN[0] + strN[1]]))
        R = improR * np.array(rot_mat[strN[0] + strN[1]])
        Axis = strN[1]
        if len(strN) > 2:
            for j in strN[2:]:
                T.append(trans_sym[j])
    
    elif strN[0] == '3' and strN[1] == '*':
        Axis = '*'
        R1.append(improR * np.array(rot_mat[strN[0] + Axis]))  # 3*
        R = improR * np.array(rot_mat[strN[0] + Axis])
    
    else:
        Axis = 'z'
        R1.append(improR * np.array(rot_mat[strN[0] + Axis]))
        R = improR * np.array(rot_mat[strN[0] + Axis])
        
        for j in strN[1:]:
            if j in ['1', '2', '3', '4', '5']:
                T.append([0, 0, float(j)/float(strN[0])])
            else:
                T.append(trans_sym[j])

    T = np.array(T, dtype=float)
    count = 0
    eOper = []
    allT = np.zeros(3, dtype=float)
    for t in T:
        allT = allT + t
    eOper.append([count, R, allT, Axis])
    T1.append(allT)

    if len(N) > 1:
        # print(preR)
        N1 = oper[2:]

        for k, N in enumerate(N1):
            T = []
            if k == 0:
                if N.startswith('-'):
                    improR = -1
                    N = N.strip('-')
                
                strN = N
                if strN[0] == '2':
                    if len(strN) > 1 and strN[1] == '=':
                        Axis = 'z'
                        R1.append(improR * np.array(rot_mat[strN[0] + 'pp' + Axis]))  # 2ppz
                        R = improR * np.array(rot_mat[strN[0] + 'pp' + Axis])
                        if len(strN) > 2:
                            for j in strN[2:]:
                                T.append(trans_sym[j])
                    
                    elif preR == '2' or preR == '4':
                        Axis = 'x'
                        R1.append(improR * np.array(rot_mat[strN[0] + Axis]))  # 2x
                        R = improR * np.array(rot_mat[strN[0] + Axis])
                        if len(strN) > 1:
                            for j in strN[1:]:
                                T.append(trans_sym[j])
                    
                    elif preR == '3' or preR == '6':
                        Axis = 'z'
                        R1.append(improR * np.array(rot_mat[strN[0] + 'p' + Axis]))  # 2pz
                        R = improR * np.array(rot_mat[strN[0] + 'p' + Axis])
                        if len(strN) > 1:
                            for j in strN[1:]:
                                T.append(trans_sym[j])
                
                else:
                    Axis = 'z'
                    R1.append(improR * np.array(rot_mat[strN[0] + Axis]))  #
                    R = improR * np.array(rot_mat[strN[0] + Axis])
                    if len(strN) > 1:
                        for j in strN[1:]:
                            T.append(trans_sym[j])

            elif k == 1:
                if N.startswith('-'):
                    improR = -1
                    N = N.strip('-')
                strN = N
                
                if strN[0] == '3':
                    Axis = '*'
                    R1.append(improR * np.array(rot_mat[strN[0] + Axis]))  # 3*
                    R = improR * np.array(rot_mat[strN[0] + Axis])

                else:
                    Axis = 'z'
                    R1.append(improR * np.array(rot_mat[strN[0] + Axis]))
                    R = improR * np.array(rot_mat[strN[0] + Axis])
                
                if len(N) > 1:
                    for j in N[1:]:
                        T.append(trans_sym[j])
            
            else:
                if N.startswith('-'):
                    improR = -1
                    N = N.strip('-')
                
                strN = N
                Axis = 'z'
                R1.append(improR * np.array(rot_mat[strN[0] + Axis]))
                R = improR * np.array(rot_mat[strN[0] + Axis])
                if len(N) > 1:
                    for j in N[1:]:
                        T.append(trans_sym[j])

            T = np.array(T, dtype=float)
            allT = np.zeros(3, dtype=float)
            for t in T:
                allT = allT + t
            
            T1.append(allT)
            count = k + 1
            eOper.append([count, R, allT, Axis])

    return L, V, eOper, R1, T1


# recusive function
def mulGrouOper(GR, GT):      
    if not (GR[-1] == E).all():
        R = np.dot(GR[0], GR[-1])
        T = np.dot(GR[0], GT[-1]) + GT[0]
        GR.append(R)
        GT.append(T)
        GR, GT = mulGrouOper(GR, GT)
    return GR, GT


def GroupOperation(hallnum):
    L, V, eOper, R1, T1 = LookHallSymbol(hallnum)
    T0 = np.zeros(3, dtype=float)

    if L.startswith('-'):
        GR_2 = [E, -E]
        GT_2 = [T0, T0]
        L = L.strip('-')
    else:
        GR_2 = [E]
        GT_2 = [T0]

    for r, t in zip(R1, T1):
        GR_1, GT_1 = mulGrouOper([r], [t])
        GR_3 = []
        GT_3 = []
        for r1, t1 in zip(GR_1, GT_1):
            for r2, t2 in zip(GR_2, GT_2):
                r3 = np.dot(r1, r2)
                t3 = np.dot(r1, t2) + t1
                GR_3.append(r3)
                GT_3.append(t3)
        GR_2 = GR_3
        GT_2 = GT_3

    GR_new = [np.array(i, dtype=int) for i in GR_2]
    GT_new = []

    if V is not None:
        V1 = np.array(V/12.0)
        #print(V1)
        for r, t in zip(GR_new, GT_2):
            GT_new.append(t + np.dot((E - r), V1))
    else:
        GT_new = GT_2

    GR_conv = []
    GT_conv = []
    for l in latt_sym[L]:
        for r, t in zip(GR_new, GT_new):
            GR_conv.append(r)
            t = t + l
            for i, j in enumerate(t):
                if j >= 1.0:
                    t[i] = j - 1.0
                elif j <= -1.0:
                    t[i] = j + 1.0
            GT_conv.append(t)

    return GR_conv, GT_conv


# check LookHallSymbol()
def test_getInfoFromHallSym(num):
    L, V, eOper, R1, T1 = LookHallSymbol(408)
    #print(R1, T1)
    GR_111 = []
    GT_111 = []
    for r, t in zip(R1, T1):
        GR, GT = mulGrouOper([r], [t])
        #print(GR, GT)
        for i in GR:
            GR_111.append(i)
        for j in GT:
            GT_111.append(j)

    print(GR_111, GT_111)


# get symmetry operations from hall symbol
def test_symOperFromHallSym(num):
    GR, GT = GroupOperation(num)  
    print('rotation: ')
    for i in GR:
        print(i)
    print('translations: ')
    for i in GT:
        print(i)

    
# check operations of space group refering to spglib code
def teat_countAllOperations():
    count = 0
    for i in range(1, 531):
        GR, GT = GroupOperation(i)
        reference = spglib.get_symmetry_from_database(i)
        
        a1 = []     # symmetry operations from spglib code
        a2 = []     # symmetry operations from this code
        
        for j1, j2 in zip(reference['rotations'], reference['translations']):
            m = []
            for k1 in j1:
                m.append(list(k1))
            a1.append([m, list(j2)])

        for j1, j2 in zip(GR, GT):
            m = []
            for k1 in j1:
                m.append(list(k1))
            a2.append([m, list(j2)])
        # print(a2)

        a = [m for m in a1 if m not in a2]
        b = [n for n in a2 if n not in a1]

        if not ((a is not None) and (b is not None)):
            print('error:', i)
        else:
            count = count + 1
            # print('true')
        print(count)
