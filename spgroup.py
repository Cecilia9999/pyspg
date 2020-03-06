#!/usr/bin/env python
# -*- coding=utf-8 -*-

'''
This part is for finding space group
'''

import numpy as np
import delaunay
import primitive
import pointgp
import csv
import matchHall
import niggli

spg_to_hall_num = [
        1, 2, 3, 6, 9, 18, 21, 30, 39, 57,
        60, 63, 72, 81, 90, 108, 109, 112, 115, 116,
        119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
        155, 161, 164, 170, 173, 176, 182, 185, 191, 197,
        203, 209, 212, 215, 218, 221, 227, 228, 230, 233,
        239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
        290, 292, 298, 304, 310, 313, 316, 322, 334, 335,
        337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
        355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
        368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
        378, 379, 380, 381, 382, 383, 384, 385, 386, 387,
        388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
        398, 399, 400, 401, 402, 404, 406, 407, 408, 410,
        412, 413, 414, 416, 418, 419, 420, 422, 424, 425,
        426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
        440, 441, 442, 443, 444, 446, 447, 448, 449, 450,
        452, 454, 455, 456, 457, 458, 460, 462, 463, 464,
        465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
        475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
        485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
        495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
        507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
        517, 518, 520, 521, 523, 524, 525, 527, 529, 530
    ]


def GetCenter(trans_matP, laue):   # //trans_matP is vertical
    identity = np.array([[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]])
    mono_i2c = np.array([[1, 0, -1],
                         [0, 1, 0],
                         [1, 0, 0]])
    mono_a2c = np.array([[0, 0, 1],
                         [0, -1, 0],
                         [1, 0, 0]])
    rho_obv = np.array([[2./3, -1./3, -1./3],
                       [1./3, 1./3, -2./3],
                       [1./3, 1./3, 1./3]])
    rho_rev = np.array([[1./3, -2./3, 1./3],
                       [2./3, -1./3, -1./3],
                       [1./3, 1./3, 1./3]])
    a2c = np.array([[0, 0, 1],
                    [1, 0, 0],
                    [0, 1, 0]])
    b2c = np.array([[0, 1, 0],
                    [0, 0, 1],
                    [1, 0, 0]])

    correct_mat = np.identity(3)
    det = abs(np.linalg.det(trans_matP))
    #print(det)
    center = None
    # // see TABLE III
    if det == 1:
        center = 'PRIMITIVE'

    elif det == 2:
        center = 'PRIMITIVE'  # // complete get center
        flag = []
        for i in range(3):
            # // A center
            # // exist one row, |a0| + |a1| + |a2| = 1 and |a0| = 1
            if trans_matP[i][1] == 0 and trans_matP[i][2] == 0 and abs(trans_matP[i][0]) == 1:
                center = 'A_FACE'
                break
            # // B center
            elif trans_matP[i][0] == 0 and trans_matP[i][2] == 0 and abs(trans_matP[i][1]) == 1:
                center = 'B_FACE'
                break
            # // C center
            # // exist one row, |a0| + |a1| + |a2| = 1 and |a2| = 1
            elif trans_matP[i][0] == 0 and trans_matP[i][1] == 0 and abs(trans_matP[i][2]) == 1:
                center = 'C_FACE'
                break
            # // Body center
            # // for each row, |a0| + |a1| + |a2| = 2
            elif (abs(trans_matP[i][0]) + abs(trans_matP[i][1]) + abs(trans_matP[i][2])) == 2:
                flag.append(True)
                if len(flag) == 3:
                    center = 'BODY'
                    break

        if center == 'A_FACE':
            center = 'C_FACE'
            print('A to C...')
            if laue == 'LAUE2M':
                correct_mat = mono_a2c.copy()
            else:
                correct_mat = a2c.copy()
        elif center == 'B_FACE':
            center = 'C_FACE'
            print('B to C...')
            correct_mat = b2c.copy()
        elif center == 'BODY' and laue == 'LAUE2M':
            center = 'C_FACE'
            print('Monoclinic I to C...')
            correct_mat = mono_i2c.copy()

    elif det == 3:
        center = 'R_CENTER'
        tp_mat1 = np.dot(trans_matP, rho_obv)   # // M' M
        tp_mat = tp_mat1.copy()
        # // determine whether in obverse setting
        # // see the tp_mat is int mat or not
        for j in range(3):
            for k in range(3):
                if tp_mat[j, k] < 0.0:
                    tp_mat[j, k] = int(tp_mat[j, k] - 0.5)
                else:
                    tp_mat[j, k] = int(tp_mat[j, k] + 0.5)

        tp_flag = True
        for j in range(3):
            for k in range(3):
                if abs(tp_mat[j, k] - tp_mat1[j, k]) > 0.00001:
                    tp_flag = False
                    break
                else:
                    tp_flag = True
        if tp_flag:
            correct_mat = rho_obv.copy()
        else:
            # // determine whether in reverse setting
            # // see the tp_mat is int mat or not
            tp_mat2 = np.dot(trans_matP, rho_rev)
            tp_mat = tp_mat2.copy()
            for j in range(3):
                for k in range(3):
                    if tp_mat[j, k] < 0.0:
                        tp_mat[j, k] = int(tp_mat[j, k] - 0.5)
                    else:
                        tp_mat[j, k] = int(tp_mat[j, k] + 0.5)

            tp_flag = True
            for j in range(3):
                for k in range(3):
                    if abs(tp_mat[j, k] - tp_mat2[j, k]) > 0.00001:
                        tp_flag = False
                        break
                    else:
                        tp_flag = True
                    if tp_flag:
                        correct_mat = rho_rev.copy()
                    else:
                        print('Cannot find correction matrix!')

    elif det == 4:
        center = 'FACE'

    return center, correct_mat


def GetConvSymm(center, trans_mat, symmetry):
    conv_symm = []
    for sy in symmetry:
        #tp_rot1 = np.dot(trans_mat, sy[0])
        #conv_rot = np.dot(tp_rot1, np.linalg.inv(trans_mat))  # // new_rot W'= P^-1 W P   P: vertical
        tp_rot1 = np.dot(np.linalg.inv(trans_mat), sy[0])
        conv_rot = np.dot(tp_rot1, trans_mat)  # // new_rot W'= P^-1 W P
        # // let W' with elements are all integers
        for j in range(3):
            for k in range(3):
                if conv_rot[j, k] < 0.0:
                    conv_rot[j, k] = int(conv_rot[j, k] - 0.5)
                else:
                    conv_rot[j, k] = int(conv_rot[j, k] + 0.5)

        conv_trans = np.dot(np.linalg.inv(trans_mat), sy[1])
        #print(conv_symm)
        conv_symm.append([np.array(conv_rot, dtype=int), conv_trans])

    #print('check conv_symm:', len(conv_symm))

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
        shift.append([2./3, 1./3, 1./3])
        shift.append([1./3, 2./3, 2./3])

    elif center == 'FACE':
        shift.append([0, 0.5, 0.5])
        shift.append([0.5, 0, 0.5])
        shift.append([0.5, 0.5, 0])

    new_conv_rot = np.zeros((3, 3))
    new_conv_trans = np.zeros(3)
    new_conv_sym = []

    for i in range(len(shift)):
        for j in range(len(conv_symm)):
            new_conv_rot = conv_symm[j][0].copy()
            new_conv_trans = conv_symm[j][1] + shift[i]
            tp2 = new_conv_trans.copy()
            for k in range(3):
                if tp2[k] < 0.0:
                    tp2[k] = int(tp2[k] - 0.5)
                else:
                    tp2[k] = int(tp2[k] + 0.5)
            new_conv_trans = new_conv_trans - tp2
            for k in range(3):
                if new_conv_trans[k] < 0.0:
                    new_conv_trans[k] += 1.0
            new_conv_sym.append([new_conv_rot, new_conv_trans])
            #print('check int2', new_conv_sym[0])
    #print('check new_conv_symm: ', len(new_conv_sym))
    return new_conv_sym


def GetSpaceData(hall_num):
    spgtype = None
    with open("../supports/newspg.csv", "rt") as CF:
        cr = csv.reader(CF)
        for i, j in enumerate(cr):
            if i == hall_num:
                choice = j[2]
                schoenf = j[4]
                hall = j[5]
                international = j[6]
                intern_full = j[7]
                pg_num = j[-1]
                spgtype = [choice, schoenf, hall, international, intern_full, int(pg_num)]
                break
        CF.close()

    return spgtype


def GetHallNumber(conv_latt, spgn, pgt_type, conv_symm, center):
    change_of_b_501 = np.array([[0, 0, 1],
                                [0, -1, 0],
                                [1, 0, 0]])

    hR_to_hP = np.array([[1, 0, 1],
                         [-1, 1, 1],
                         [0, -1, 1]])

    hall_num = spg_to_hall_num[spgn]
    spgtype = GetSpaceData(hall_num)
    num_hall_type = spg_to_hall_num[spgn] - spg_to_hall_num[spgn - 1]
    #print('check', pgt_type[0], spgtype[5])
    if spgtype[5] != pgt_type[0]:
        return 0, None, None, None

    if pgt_type[4] == 'MONOCLI':
        flag, new_latt, origin = matchHall.MatchHallMono(conv_latt, hall_num, num_hall_type, conv_symm, center)
        if flag:
            return 1, new_latt, origin, spgn

    if pgt_type[4] == 'ORTHO':
        if spgn in [48, 50, 59, 68, 70]:
            num_hall_type = num_hall_type / 2

        if num_hall_type == 1:
            num_free_ax = 6
            #print('num_hall_type 1')
            flag, new_latt, origin = matchHall.MatchHallOrtho(conv_latt, hall_num, conv_symm, center, num_free_ax)
            if flag:
                return 1, new_latt, origin, spgn

        elif num_hall_type == 2:
            num_free_ax = 3
            #print('num_hall_type 2')
            flag, new_latt, origin = matchHall.MatchHallOrtho(conv_latt, hall_num, conv_symm, center, num_free_ax)
            if flag:
                return 1, new_latt, origin, spgn

        elif num_hall_type == 3:
            #print('num_hall_type 3_1')
            tp_conv_latt = conv_latt.copy()
            num_free_ax = 0
            flag, new_latt, origin = matchHall.MatchHallOrtho(tp_conv_latt, spg_to_hall_num[spgn - 1], conv_symm,
                                                                center, num_free_ax)
            # // strange!!!!!!!!!!!need to be changed!!!!!!!!!
            if flag:
                return 1, new_latt, origin, spgn-1

            #print('num_hall_type 3_2')
            trans_m = np.transpose(np.dot(tp_conv_latt, np.linalg.inv(conv_latt)))
            tp_conv_symm = GetConvSymm('PRIMITIVE', trans_m, conv_symm)
            num_free_ax = 2
            flag, new_latt, origin = matchHall.MatchHallOrtho(tp_conv_latt, hall_num, tp_conv_symm,
                                                               center, num_free_ax)
            if flag:
                return 1, new_latt, origin, spgn

        elif num_hall_type == 6:
            #print('num_hall_type 6')
            num_free_ax = 1
            flag, new_latt, origin = matchHall.MatchHallOrtho(conv_latt, hall_num, conv_symm, center, num_free_ax)
            if flag:
                return 1, new_latt, origin, spgn

    if pgt_type[4] == 'CUBIC':
        flag, origin = matchHall.MatchHall(conv_latt, hall_num, conv_symm, center)
        if flag:
            #print('check flag', flag)
            return 1, conv_latt, origin, spgn
        if hall_num == 501:
            tp_conv_latt = np.dot(np.transpose(change_of_b_501), conv_latt)
            tp_conv_symm = GetConvSymm('PRIMITIVE', change_of_b_501, conv_symm)
            flag1, origin = matchHall.MatchHall(tp_conv_latt, hall_num, tp_conv_symm, 'PRIMITIVE')
            if flag1:
                return 1, tp_conv_latt, origin, spgn

    if pgt_type[4] == 'TRIGO':
        if center == 'R_CENTER':
            if hall_num in [433, 436, 444, 450, 452, 458, 460]:
                tp_conv_latt = np.dot(np.transpose(hR_to_hP), conv_latt)
                tp_conv_symm = GetConvSymm('R_CENTER', hR_to_hP, conv_symm)
                flag, origin = matchHall.MatchHall(tp_conv_latt, hall_num, tp_conv_symm, 'PRIMITIVE')
                if flag:
                    return 1, tp_conv_latt, origin, spgn
            else:
                flag, origin = matchHall.MatchHall(conv_latt, hall_num, conv_symm, 'PRIMITIVE')
                if flag:
                    return 1, conv_latt, origin, spgn

    flag, origin = matchHall.MatchHall(conv_latt, hall_num, conv_symm, center)
    if flag:
        return 1, conv_latt, origin, spgn

    return 0, None, None, None


def GetSpaceGroup(prim_latt, symmetry, pgt_type, axis):
    # // search hall number
    # // standardize
    #o_latt, o_pos, o_num, dictp = delaunay.StructRead()
    #prim_latt, num, pos = primitive.GetPrimLattice(o_latt, o_pos, o_num, dictp)
    #symmetry, pgt_type, axis = pointgp.GetPointGroup(prim_latt, num, pos, dictp)  # // pgt_type = [pgnum, pgtable, symbol, schoenf, holo, Laue]

    laue = pgt_type[5]
    trans_matP = (np.transpose(axis)).copy()   # // P : vertical vectors
    # print(axis, trans_matP)
    if pgt_type[0] == 0:
        print('Error! No symmetry!')
        return None, None, None, None, None, None    # // return 0

    #conv_latt = np.zeros((3, 3))
    conv_latt = np.dot(axis, prim_latt)   # // axis = np.transpose(trans_matP)

    if laue == 'LAUE1':      # // TRICLINIC

        #   // get new transformation matrix
        #   // change conv to min using niggli, then find trans_nat of prim to min
        nigg = niggli.Niggli(conv_latt)
        min_latt = np.zeros((3, 3))
        if not (nigg == np.zeros((3, 3))).all():
            for i in range(3):
                for j in range(3):
                    min_latt[i, j] = nigg[i, j]

        if np.linalg.det(min_latt) < 0:
            min_latt = -min_latt

        inv_latt = np.linalg.inv(np.transpose(prim_latt))
        trans_matP = np.dot(inv_latt, np.transpose(min_latt))   # // here, P is vertical vec
        for i in range(3):
            for j in range(3):
                if trans_matP[i, j] < 0.0:
                    trans_matP[i, j] = int(trans_matP[i, j] - 0.5)
                else:
                    trans_matP[i, j] = int(trans_matP[i, j] + 0.5)
        # // return trans_matP

    if laue == 'LAUE2M':        # // MONOCLINIC

        #   // get new transformation matrix
        #   // change conv to min using delaunay, then find trans_nat of prim to min
        flag, min_latt, delauP = delaunay.Delaunay(conv_latt, 1)   # // delauP is vertical
        if flag:
            inv_latt = np.linalg.inv(prim_latt)
            trans_matP = np.transpose(np.dot(min_latt, inv_latt))
            for i in range(3):
                for j in range(3):
                    if trans_matP[i, j] < 0.0:
                        trans_matP[i, j] = int(trans_matP[i, j] - 0.5)
                    else:
                        trans_matP[i, j] = int(trans_matP[i, j] + 0.5)
        # // return trans_matP

    center = None
    center, correct_mat = GetCenter(trans_matP, laue)

    #print(trans_matP)
    print('correct matrix: \n', correct_mat)
    print('center: \n', center)

    new_trans_mat = trans_matP.copy()
    if center is not None:
        new_trans_mat = np.dot(trans_matP, correct_mat)  # // trans_P, correct_mat and new_trans_mat are vertical vecs
        conv_latt = np.dot(np.transpose(new_trans_mat), prim_latt)  # // conv_latt : horizontal vectors

    #print('conven latt1:\n', conv_latt)
    print('new trans mat: \n', new_trans_mat)
    if center == 'R_CENTER':
        conv_symm = GetConvSymm('PRIMITIVE', new_trans_mat, symmetry)
    else:
        conv_symm = GetConvSymm(center, new_trans_mat, symmetry)

    #print('check numsym:', len(conv_symm))
    spgnum = None
    hallnum = None
    origin = np.zeros(3)
    new_latt = np.zeros((3, 3))

    if conv_symm != []:
        for i in range(231):
            flag, new_latt, origin, spgnum = GetHallNumber(conv_latt, i, pgt_type, conv_symm, center)
            if flag:
                print('Congratulations! Find space group!')
                hallnum = spg_to_hall_num[spgnum]
                return spgnum + 1, hallnum, origin, new_latt, new_trans_mat, center
                #break

        print('Error! Cannot match any space group data!')
        return None, None, None, None, None, None

    '''
    print('find sym!')
    print('spg number: ', spgnum+1)
    print('hall num: ', hallnum)
    print(new_trans_mat)
    '''