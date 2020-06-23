#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Created on 2019-04-06
Update  on 2020-06-23
Author: Cecilia9999
GitHub: https://github.com/Cecilia9999/
"""

'''
    This file is for handling input cell
'''

import os
import numpy as np
import primitive
import delaunay
import spgroup
import pointgp

A_mat = np.array ([[1, 0, 0],
                   [0, 1. / 2, -1. / 2],
                   [0, 1. / 2, 1. / 2]])

B_mat = np.array ([[1. / 2, 0, -1. / 2],
                   [0, 1, 0],
                   [1. / 2, 0, 1. / 2]])

C_mat = np.array ([[1. / 2, 1. / 2, 0],
                   [-1. / 2, 1. / 2, 0],
                   [0, 0, 1]])

R_mat = np.array ([[2. / 3, -1. / 3, -1. / 3],
                   [1. / 3, 1. / 3, -2. / 3],
                   [1. / 3, 1. / 3, 1. / 3]])

I_mat = np.array ([[-1. / 2, 1. / 2, 1. / 2],
                   [1. / 2, -1. / 2, 1. / 2],
                   [1. / 2, 1. / 2, -1. / 2]])

F_mat = np.array ([[0, 1. / 2, 1. / 2],
                   [1. / 2, 0, 1. / 2],
                   [1. / 2, 1. / 2, 0]])


def StructRead(filename):
    if os.path.isfile(filename):
        print("Using POSCAR to calculate...")
    else:
        print("Error!\nPlease prepare a POSCAR!")

    with open(filename, 'r') as POS:
        tmp = POS.readlines()
        tmp1 = tmp[1].split()
        latt_const = float(tmp1[0])
        # print(latt_const)
        
        tmp3 = []
        for i in tmp[2:5]:
            tmp2 = []
            for j in i.split():
                tmp2.append(float(j) * latt_const)
            tmp3.append(tmp2)
        
        lattice = np.array(tmp3)
        np.set_printoptions(formatter={'float': '{:0.10f}'.format})
        print('lattice:\n', lattice)

        species = tmp[5]
        tmp7 = species.split()  # [Si, O]
        tmp4 = [int(i) for i in tmp[6].split()]  # [4, 8]
        dictp = dict(zip(tmp7, tmp4))
        
        # print(dictp)
        numbers = []
        for i in dictp.keys():
            for k in range(dictp[i]):
                # numbers.append(GetAtom(i))  # [40, 40, 40, 40, 8, 8, 8, 8, 8, 8, 8, 8]
                numbers.append(i)  # [40, 40, 40, 40, 8, 8, 8, 8, 8, 8, 8, 8]
        print('numbers:\n', numbers)

        tmp5 = tmp[7].split()
        m = 0
        if tmp5[0].startswith('D') or tmp5[0].startswith('d') or tmp5[0].startswith('C') or tmp5[0].startswith('c'):
            m = 8
        elif tmp5[0].startswith('S') or tmp5[0].startswith('s'):
            m = 9

        pos = []
        for j in tmp[m:]:
            tmp6 = []
            for k in j.split():
                if k != "\n":
                    tmp6.append(float(k))
                if len(tmp6) == 3:
                    break
            pos.append(tmp6)
            if len(pos) == len(numbers):
                break

        positions = [i for i in pos if i != []]
        print('positions:')
        print(positions)

        return lattice, positions, numbers, dictp


# determine if the input cell is reasonable or not (check distortion)
def AtomsTooClosed(latt, pos, num):
    tolerance = 0.00001
    cel_size = len(pos)
    check_pos = []
    for i in range(cel_size):
        tp2 = np.array(pos[i]).copy()
        tp_pos = np.array(pos[i]).copy()
        for j in range(3):
            if tp_pos[j] < 0.0:
                tp_pos[j] = int(tp_pos[j] - 0.5)
            else:
                tp_pos[j] = int(tp_pos[j] + 0.5)
        tp2 = tp2 - tp_pos
        # print(tp2)
        for j in range(3):
            if tp2[j] < 0.0:
                tp2[j] += 1.0
        check_pos.append(tp2)
    
    # print(check_pos)
    for i in range(cel_size):
        for j in range(cel_size):
            if i == j:
                continue
            if num[i] == num[j]:
                if primitive.VecIsOverlap(np.array(check_pos[i]), np.array(check_pos[j]), latt, tolerance):
                    return 1, None

    return 0, check_pos


def StandardPrimCell(latt, pos, num, center):
    flag, ch_pos = AtomsTooClosed(latt, pos, num)
    if flag:
        return None, None, None

    if center is None:
        return None, None, None

    std_prim_latt = np.zeros((3, 3))
    std_prim_pos = []
    std_prim_num = []

    o_trans_mat = np.identity(3)
    inv_trans_mat = np.linalg.inv(o_trans_mat)
    trans_mat = np.zeros((3, 3))

    if center == 'PRIMITIVE':
        trans_mat = inv_trans_mat.copy()
    elif center == 'BODY':
        trans_mat = np.dot(inv_trans_mat, I_mat)
    elif center == 'FACE':
        trans_mat = np.dot(inv_trans_mat, F_mat)
    elif center == 'A_FACE':
        trans_mat = np.dot(inv_trans_mat, A_mat)
    elif center == 'C_FACE':
        trans_mat = np.dot(inv_trans_mat, C_mat)
    elif center == 'R_CENTER':
        trans_mat = np.dot(inv_trans_mat, R_mat)

    prim_latt = np.dot(np.transpose(trans_mat), latt)
    # t = np.dot(trans_matrix, trans_mat)

    std_prim_latt, std_prim_num, std_prim_pos = primitive.AtomPosConv2Prim(prim_latt, latt, pos, num)  
    # pos --> check_pos?

    if not (std_prim_latt == np.zeros((3, 3))).all():
        # adjust positions to list format
        new_std_prim_pos = []
        for i in std_prim_pos:
            tp = []
            for j in i:
               tp.append(j)
            new_std_prim_pos.append(tp)
        print('standard lattice:\n', std_prim_latt)
        print('standard numbers:\n', std_prim_num)
        print('standard positions:\n', new_std_prim_pos)
        return std_prim_latt, std_prim_num, new_std_prim_pos
    else:
        print('Error! Cannot change to the standard primitive cell!')
        return None, None, None


def GetStandardPos(std_latt, pos, num, trans_matrix, center, is_prim, origin, prim_flag=False):
    # Get standard positions
    tp_positions = []
    positions = []
    std_num = []
    std_pos = []
    if is_prim:
        tp1_positions = delaunay.ChangeOfBasis(pos, np.linalg.inv(trans_matrix))
        tp = np.zeros(3)
        for i in tp1_positions:
            tp_positions.append(np.array(i) + origin)  # // Xc~ = Q * Xp + Pc

        shift = [[0, 0, 0]]
        if center == 'A_FACE':
            shift.append([0, 0.5, 0.5])

        elif center == 'B_FACE':
            shift.append([0.5, 0, 0.5])

        elif center == 'C_FACE':
            shift.append([0.5, 0.5, 0])

        elif center == 'BODY':
            shift.append([0.5, 0.5, 0.5])

        elif center == 'R_CENTER':
            shift.append([2. / 3, 1. / 3, 1. / 3])
            shift.append([1. / 3, 2. / 3, 2. / 3])

        elif center == 'FACE':
            shift.append([0, 0.5, 0.5])
            shift.append([0.5, 0, 0.5])
            shift.append([0.5, 0.5, 0])

        # new_pos = np.zeros(3)
        for j in range(len(tp_positions)):
            for i in range(len(shift)):
                new_pos = np.array(tp_positions[j]) + np.array(shift[i])
                tp2 = new_pos.copy()
                for k in range(3):
                    if tp2[k] < 0.0:
                        tp2[k] = int(tp2[k] - 0.5)
                    else:
                        tp2[k] = int(tp2[k] + 0.5)
                new_pos = new_pos - tp2
                for k in range(3):
                    if new_pos[k] < 0.0:
                        new_pos[k] += 1.0
                positions.append(new_pos)  # // Xc = Xc~ + shift vec
        # print('cek1', positions)

        for j in num:
            for i in range(len(shift)):
                std_num.append(j)

        # adjust positions to list format
        for i in positions:
            tp = []
            for j in i:
                tp.append(j)
            std_pos.append(tp)

    else:
        # input is conventional
        tp_positions = delaunay.ChangeOfBasis(pos, np.linalg.inv(trans_matrix))
        for i in tp_positions:
            positions.append(np.array(i) + origin)

        # adjust positions to list format
        for i in positions:
            tp = []
            for j in i:
                tp.append(j)
            std_pos.append(tp)

        std_num = num.copy()
        # print('cek1', positions)

    if is_prim or prim_flag:
        print('standard conventional lattice:\n', std_latt)
        print('standard conventional numbers:\n', std_num)
        print('standard conventional positions:\n', std_pos)
        return StandardPrimCell(std_latt, std_pos, std_num, center)
    else:
        return std_latt, std_num, std_pos


def GetStandardCell(latt, pos, num, trans_matrix, center, is_prim, pgt_type, origin, prim_flag=False):

    tolerance = 0.00001
    transf = np.zeros((3, 3))
    std_latt = np.zeros((3, 3))
    la = np.linalg.norm(latt[0])
    lb = np.linalg.norm(latt[1])
    lc = np.linalg.norm(latt[2])
    print(la, lb, lc)
    alpha = np.arccos(np.dot(latt[1], latt[2]) / (lb * lc))
    beta = np.arccos(np.dot(latt[2], latt[0]) / (lc * la))
    gamma = np.arccos(np.dot(latt[0], latt[1]) / (la * lb))
    tp_list = [[latt[0], la, 0], [latt[1], lb, 1], [latt[2], lc, 2]]
    srt_l = sorted(tp_list, key=lambda k: k[1])
    print('srt_l', srt_l)

    if pgt_type[4] == 'CUBIC' or pgt_type[4] == 'ORTHO':
        if center == 'C_FACE':
            transf[2] = [0, 0, 1]
            srt_l = sorted(tp_list[:2], key=lambda k: k[1])
            for i in range(2):
                transf[i][srt_l[i][2]] = 1
            std_latt[0] = [srt_l[0][1], 0, 0]
            std_latt[1] = [0, srt_l[1][1], 0]
            std_latt[2] = [0, 0, lc]
        elif center == 'A_FACE':    # // change to C-centering to match Setyawan/Curtarolo convention
            transf[2] = [1, 0, 0]
            srt_l = sorted(tp_list[1:], key=lambda k: k[1])
            for i in range(2):
                transf[i][srt_l[i][2]] = 1
            std_latt[0] = [srt_l[0][1], 0, 0]
            std_latt[1] = [0, 0, srt_l[1][1]]
            std_latt[2] = [0, 0, la]
        else:
            for i in range(3):
                transf[i][srt_l[i][2]] = 1
            for i in range(3):
                std_latt[i][i] = srt_l[i][1]

    elif pgt_type[4] == 'TETRA':
        nla = srt_l[0][1]
        nlb = srt_l[1][1]
        nlc = srt_l[2][1]
        for i in range(3):
            transf[i][srt_l[i][2]] = 1
        if abs(nlb - nlc) < tolerance < abs(nla - nlc):
            std_latt[0] = [nlc, 0, 0]
            std_latt[1] = [0, nlb, 0]
            std_latt[2] = [0, 0, nla]
            transf = np.dot([[0, 0, 1], [0, 1, 0], [1, 0, 0]], transf)
        else:
            std_latt[0] = [nla, 0, 0]
            std_latt[1] = [0, nlb, 0]
            std_latt[2] = [0, 0, nlc]

    elif pgt_type[4] == 'HEXA' or pgt_type[4] == 'TRIGO':
        # for the conventional cell representation, show the rhombohedral lattices as hexagonal
        nla = nlb = nlc = 0
        
        # a = b = c rhombohedral
        if np.all(np.abs([la - lb, lc - lb, la - lc]) < 0.001):     
            tp_latt = np.zeros((3, 3))
            tp_latt[0] = latt[0] - latt[1]
            tp_latt[1] = latt[1] - latt[2]
            tp_latt[2] = latt[0] + latt[1] + latt[2]
            nla = np.linalg.norm(tp_latt[0])
            nlb = np.linalg.norm(tp_latt[1])
            nlc = np.linalg.norm(tp_latt[2])
            tp_list = [[tp_latt[0], nla, 0], [tp_latt[1], nlb, 1], [tp_latt[2], nlc, 2]]
            srt_l = sorted(tp_list, key=lambda k: k[1])
            la = srt_l[0][1]
            lb = srt_l[1][1]
            lc = srt_l[2][1]

        if abs(lb - lc) < 0.001:
            std_latt[0] = [lc / 2, -lc * np.sqrt(3) / 2, 0]
            std_latt[1] = [lc / 2, lc * np.sqrt(3) / 2, 0]
            std_latt[2] = [0, 0, la]
        
        else:
            std_latt[0] = [la / 2, -la * np.sqrt(3) / 2, 0]
            std_latt[1] = [la / 2, la * np.sqrt(3) / 2, 0]
            std_latt[2] = [0, 0, lc]
        
        transf = np.eye(3, 3)

    elif pgt_type[4] == 'MONOCLI':
        # C-setting
        # Please check A/B-FACE, are they needed to be set?
        if center == 'C_FACE':
            transf[2] = [0, 0, 1]
            srt_l = sorted(tp_list[0:2], key=lambda k: k[1])
            nla = srt_l[0][1]
            nlb = srt_l[1][1]
            nlc = lc
            tpp_latt = latt.copy()
            # for t in itertools.permutations (list (range (2)), 2):   
            # 2, 0, 1;  0, 1, 2
            for t in [[0, 1], [1, 0]]:
                tp_latt = np.zeros((3, 3))
                std_latt = np.zeros((3, 3))
                nnla = np.linalg.norm(tpp_latt[t[0]])
                nnlb = np.linalg.norm(tpp_latt[t[1]])
                nnlc = np.linalg.norm(tpp_latt[2])
                nalpha = np.arccos(np.dot(tpp_latt[t[1]], tpp_latt[2]) / (nnlb * nnlc))
                nbeta = np.arccos(np.dot(tpp_latt[2], tpp_latt[t[0]]) / (nnlc * nnla))
                ngamma = np.arccos(np.dot(tpp_latt[t[0]], tpp_latt[t[1]]) / (nnla * nnlb))
                
                if (ngamma / np.pi * 180) == 90:
                    continue
                
                elif (ngamma / np.pi * 180) > 90:
                    # if the alpha angle is > 90 we invert a and b to get
                    # an angle < 90
                    tp_latt[0] = -tpp_latt[t[0]]
                    tp_latt[1] = -tpp_latt[t[1]]
                    tp_latt[2] = tpp_latt[2]
                    transf[0][t[0]] = -1
                    transf[1][t[1]] = -1
                    transf[2][2] = 1
                    nnla = np.linalg.norm(tp_latt[0])
                    nnlb = np.linalg.norm(tp_latt[1])
                    nnlc = np.linalg.norm(tp_latt[2])
                    ngamma = np.arccos(np.dot(tp_latt[0], tp_latt[1]) / (nnla * nnlb))
                    std_latt = np.array([[nnla, 0, 0],
                                        [nnlb * np.cos(ngamma), nnlb * np.sin(ngamma), 0],
                                        [0, 0, nnlc]])
                    tpp_latt = tp_latt.copy()
                    break

                elif (ngamma / np.pi * 180) < 90:
                    transf[0][t[0]] = 1
                    transf[1][t[1]] = 1
                    transf[2][2] = 1
                    std_latt = np.array([[nnla, 0, 0],
                                        [nnlb * np.cos(ngamma), nnlb * np.sin(ngamma), 0],
                                        [0, 0, nnlc]])
                    break

            if (std_latt == np.zeros((3, 3))).all():
                # this if is to treat the case
                # where alpha = 90 (but we still have a monoclinic
                std_latt = np.array([[nla, 0, 0],
                                    [0, nlb, 0],
                                    [0, 0, nlc]])
                for i in range(len(srt_l)):
                    transf[i][srt_l[i][2]] = 1

        else:
            # try all permutations of the axis
            # keep the ones with the non-90 BETA=beta
            # and a<c
            tpp_latt = latt.copy()
            for t in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]:
                # for t in itertools.permutations(list(range(3)), 3):
                tp_latt = np.zeros((3, 3))
                std_latt = np.zeros((3, 3))
                nnla = np.linalg.norm(tpp_latt[t[0]])
                nnlb = np.linalg.norm(tpp_latt[t[1]])
                nnlc = np.linalg.norm(tpp_latt[t[2]])
                nalpha = np.arccos(np.dot(tpp_latt[t[1]], tpp_latt[t[2]]) / (nnlb * nnlc))
                nbeta = np.arccos(np.dot(tpp_latt[t[2]], tpp_latt[t[0]]) / (nnlc * nnla))
                ngamma = np.arccos(np.dot(tpp_latt[t[0]], tpp_latt[t[1]]) / (nnla * nnlb))
                
                print('SDS', nalpha, nbeta, ngamma)
                print('long:', nnla, nnlb, nnlc, t)
                print('first angle: \n', (nbeta / np.pi * 180))
                
                if (nbeta / np.pi * 180) == 90 or nnla > nnlc:
                    continue
                
                elif (nbeta / np.pi * 180) > 90 and nnla < nnlc:
                    tp_latt[0] = -tpp_latt[t[0]]
                    tp_latt[1] = tpp_latt[t[1]]
                    tp_latt[2] = -tpp_latt[t[2]]
                    transf[0][t[0]] = -1
                    transf[1][t[1]] = 1
                    transf[2][t[2]] = -1
                    nnla = np.linalg.norm(tp_latt[0])
                    nnlb = np.linalg.norm(tp_latt[1])
                    nnlc = np.linalg.norm(tp_latt[2])
                    nbeta = np.arccos(np.dot(tp_latt[2], tp_latt[0]) / (nnlc * nnla))

                    # alpha = np.pi * nalpha / 180
                    std_latt = np.array([[nnla, 0, 0],
                                        [0, nnlb, 0],
                                        [nnlc * np.cos(nbeta), 0, nnlc * np.sin(nbeta)]])
                    tpp_latt = tp_latt.copy()
                    break

                elif (nbeta / np.pi * 180) < 90 and nnla < nnlc:
                    transf[0][t[0]] = 1
                    transf[1][t[1]] = 1
                    transf[2][t[2]] = 1
                    std_latt = np.array([[nnla, 0, 0],
                                         [0, nnlb, 0],
                                         [nnlc * np.cos(nbeta), 0, nnlc * np.sin(nbeta)]])
                    break

            if (std_latt == np.zeros((3, 3))).all():
                # this if is to treat the case
                # where alpha==90 (but we still have a monoclinic sg
                std_latt = np.array([[srt_l[0][1], 0, 0],
                                    [0, srt_l[1][1], 0],
                                    [srt_l[2][1] * np.cos(beta), 0, srt_l[2][1] * np.sin(beta)]])

                for i in range(3):
                    transf[i][srt_l[i][2]] = 1
            '''
            # if international_monoclinic:
                # The above code makes alpha the non-right angle.
                # The following will convert to proper international convention
                # that beta is the non-right angle.
            
            op = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
            transf = np.dot(op, transf)
            #std_latt = np.dot(op, std_latt)
            std_latt = np.dot(np.transpose(std_latt), op)
            #std_latt = np.dot(op, latt)
            beta = np.arccos(np.dot(std_latt[2], std_latt[0]) / (np.linalg.norm(std_latt[2]) * np.linalg.norm(std_latt[0])))
            print('beta: ', beta)
            print(std_latt)
            if beta < 90:
                op = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                transf = np.dot(op, transf)
                #std_latt = np.dot(op, std_latt)
                std_latt = np.dot(np.transpose(std_latt), op)
            print(std_latt)
            '''

    elif pgt_type[4] == 'TRICLI':
        std_latt = latt
        transf = np.identity(3)
    '''
    elif pgt_type[4] == 'TRIGO':
        std_latt = latt
        transf = np.identity(3)
    '''
    trans_matrix = np.dot(trans_matrix, np.transpose(transf))
    # print('check ', trans_matrix)

    if not (std_latt == np.zeros((3, 3))).all():
        std_lattice, std_numbers, std_positions = GetStandardPos(std_latt, pos, num, trans_matrix, center, is_prim, origin)
        return std_lattice, std_numbers, std_positions
    else:
        print('Error! Cannot change to the standard conventional cell!')
        return None, None, None


if __name__ = "__main__":
    filename = input("Please input the POSCAR name: ")
    filename = './testfile/' + filename
    o_latt, o_pos, o_num, dictp = StructRead(filename)
    prim_latt, p_num, p_pos, is_prim = primitive.GetPrimLattice(o_latt, o_pos, o_num, dictp)
    symmetry, pgt_type, axis = pointgp.GetPointGroup(prim_latt, p_num, p_pos, dictp)  
    # // pgt_type = [pgnum, pgtable, symbol, schoenf, holo, Laue]
    spgnum, hallnum, origin, new_latt, trans_matrix, center = spgroup.GetSpaceGroup(prim_latt, symmetry, pgt_type, axis)

    print('spg number: {}\n'.format(spgnum))
    print('hall num: {}\n'.format(hallnum))
    print('origin: {}\n'.format(origin))
    print('bravais lattice: {}\n'.format(new_latt))
    # print('trans matrix:\n', trans_matrix)
    print('point group: {}\n'.format(pgt_type[0]))
    print('schoenf symbol: {}\n'.format(pgt_type[2]))
    print('holohedry: {}\n'.format(pgt_type[3]))
    print('laue: {}\n'.format(pgt_type[4]))

    if is_prim:
        std_latt, std_num, std_pos = GetStandardCell(new_latt, p_pos, o_num, trans_matrix, center, is_prim, pgt_type, origin)
        # std_latt, std_num, std_pos = GetStandardPos(new_latt, p_pos, o_num, trans_matrix, center, is_prim, origin)
    else:
        trans_matrix2 = np.dot(np.transpose(np.dot(prim_latt, np.linalg.inv(o_latt))), trans_matrix)
        std_latt, std_num, std_pos = GetStandardCell(new_latt, o_pos, o_num, trans_matrix2, center, is_prim, pgt_type, origin)
        # std_latt, std_num, std_pos = GetStandardPos(new_latt, o_pos, o_num, trans_matrix2, center, is_prim, origin)

    print('lattice: {}\n'.format(std_latt))
    print('numbers: {}\n'.format(std_num))
    print('positions: {}\n'.format(std_pos))

    output = filename + "-stand-test"
    with open(output, 'w') as ST:
        ST.write("POSCAR-STAN\n" + "1.0\n")
        for i in range(3):
            for j in range(3):
                ST.write(str('{:.10f}'.format(std_latt[i][j])) + ' ')
            ST.write('\n')

        a = []
        for i in std_num:
            if i not in a:
                ST.write("  " + i + " ")
                a.append(i)
        ST.write('\n')

        b = []
        count = 0
        for j in a:
            for i in std_num:
                if j == i:
                    count += 1
            b.append(count)
            count = 0

        for i in b:
            ST.write("  " + str(i) + " ")
        ST.write('\n')

        ST.write("Direct: \n")

        for i in std_pos:
            for j in range(3):
                ST.write("  " + str('{:.12f}'.format(i[j])) + " ")
            ST.write('\n')
        ST.close()
