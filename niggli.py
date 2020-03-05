#!/usr/bin/env python
# -*- coding=utf-8 -*-
'''
This part is for finding delaunay reduced cell
Note: only for 3D lattice
'''

import numpy as np
import delaunay


count = 100
tolerance = 0

def InitialPara(latt):
    metric = np.dot(latt, np.transpose(latt))
    G = np.zeros(9)

    G[0] = metric[0, 0]
    G[4] = metric[1, 1]
    G[8] = metric[2, 2]
    G[1] = metric[0, 1]
    G[2] = metric[0, 2]
    G[3] = metric[1, 0]
    G[5] = metric[1, 2]
    G[6] = metric[2, 0]
    G[7] = metric[2, 1]

    A = G[0]
    B = G[4]
    C = G[8]
    alpha = 2.0 * G[5]
    beta = 2.0 * G[2]
    gamma = 2.0 * G[1]

    l = m = n = 0
    if alpha < -tolerance:
        l = -1
    if alpha > tolerance:
        l = 1
    if beta < -tolerance:
        m = -1
    if beta > tolerance:
        m = 1
    if gamma < -tolerance:
        n = -1
    if gamma > tolerance:
        n = 1

    return A, B, C, alpha, beta, gamma, l, m, n, G, latt


def SetPara(latt, trans_mat):

    new_latt = np.dot(np.transpose(trans_mat), latt)
    if not (new_latt == np.zeros((3, 3))).all():
        return InitialPara(new_latt)
    else:
        return None, None, None, None, None, None, None, None, None, None, None,


class Update:
    def __init__(self, A, B, C, alpha, beta, gamma, l, m, n, G, latt):
        self.A = A
        self.B = B
        self.C = C
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.l = l
        self.m = m
        self.n = n
        self.G = G
        self.latt = latt

    def update(self):
        return self.A, self.B, self.C, self.alpha, self.beta, self.gamma, self.l, self.m, self.n, self.G, self.latt


def Niggli(latt):

    A, B, C, alpha, beta, gamma, l, m, n, G, latt = InitialPara(latt)

    for c in range(count):

        if (A > B + tolerance) or ((not abs(A - B) > tolerance) and (abs(alpha) > abs(beta) + tolerance)):
            trans_mat = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]])
            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('1')

        if (B > C + tolerance) or ((not abs(B - C) > tolerance) and (abs(beta) > abs(gamma) + tolerance)):
            trans_mat = np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
            #A, B, C, alpha, beta, gamma, l, m, n, G, latt = SetPara(latt, trans_mat)
            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('2')
            continue

        if l * m * n == 1:
            i = j = k = 0
            if l == -1:
                i = -1
            else:
                i = 1
            if m == -1:
                j = -1
            else:
                j = 1
            if n == -1:
                k = -1
            else:
                k = 1
            trans_mat = np.array([[i, 0, 0], [0, j, 0], [0, 0, k]])

            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()

            print('3')


        if l == m == n == -1:
            print('4_1')
            pass
        elif l * m * n == 0 or l * m * n == -1:
            i = j = k = 1
            r = -1
            if l == 1:
                i = -1
            if l == 0:
                r = 0
            if m == 1:
                j = -1
            if m == 0:
                r = 1
            if n == 1:
                k = -1
            if n == 0:
                r = 2
            if i * j * k == -1:
                if r == 0:
                    i = -1
                elif r == 1:
                    j = -1
                elif r == 2:
                    k = -1
            trans_mat = np.array([[i, 0, 0], [0, j, 0], [0, 0, k]])

            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('4_2')

        if (abs(alpha) > B + tolerance) or ((not abs(B - alpha) > tolerance) and (2 * beta < gamma - tolerance)) or ((not abs(B + alpha) > tolerance) and (gamma < -tolerance)):
            if alpha > 0:
                trans_mat = np.array([[1, 0, 0], [0, 1, -1], [0, 0, 1]])
            elif alpha < 0:
                trans_mat = np.array([[1, 0, 0], [0, 1, 1], [0, 0, 1]])
            else:
                trans_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('5')
            continue

        if (abs(beta) > A + tolerance) or ((not abs(A - beta) > tolerance) and (2 * alpha < gamma - tolerance)) or (
                    (not abs(A + alpha) > tolerance) and (gamma < -tolerance)):
            if beta > 0:
                trans_mat = np.array([[1, 0, -1], [0, 1, 0], [0, 0, 1]])
            elif beta < 0:
                trans_mat = np.array([[1, 0, 1], [0, 1, 0], [0, 0, 1]])
            else:
                trans_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('6')
            continue

        if (abs(gamma) > A + tolerance) or ((not abs(A - gamma) > tolerance) and (2 * alpha < beta - tolerance)) or (
                    (not abs(A + gamma) > tolerance) and (beta < -tolerance)):
            if gamma > 0:
                trans_mat = np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]])
            elif gamma < 0:
                trans_mat = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
            else:
                trans_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('7')
            continue

        if (alpha + beta + gamma + A + B < -tolerance) or ((not abs(alpha + beta + gamma + A + B) > tolerance) and (2 * (A + beta) + gamma > tolerance)):
            trans_mat = np.array([[1, 0, 1], [0, 1, 1], [0, 0, 1]])

            A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1 = SetPara(latt, trans_mat)
            init = Update(A1, B1, C1, alpha1, beta1, gamma1, l1, m1, n1, G1, latt1)
            A, B, C, alpha, beta, gamma, l, m, n, G, latt = init.update()
            print('8')
            return latt

        print('repeat', c)
        return latt

    return None


'''
latt, pos, numb, dictp = delaunay.StructRead()
flag, reduc_b, delauP = delaunay.Delaunay(latt, -1)
print('red', reduc_b)
A = Niggli(reduc_b)
print(A)
'''