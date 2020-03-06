#!/usr/bin/env python
# -*- coding=utf-8 -*-

'''
This part is for finding primitive cell
'''


import numpy as np
import cardinality as cd
import delaunay


# determine two equivalent vectors or not
def VecIsOverlap(vec1, vec2, latt, tolerance):
    if tolerance is None:
        tolerance = 0.00001
    vec_diff = vec1 - vec2
    tp_diff = vec_diff.copy()
    for i in range(3):
        if tp_diff[i] < 0.0:
            tp_diff[i] = int(tp_diff[i] - 0.5)
        else:
            tp_diff[i] = int(tp_diff[i] + 0.5)
    nvec_diff = vec_diff - tp_diff
    if np.linalg.norm(np.dot(np.transpose(latt), nvec_diff)) < tolerance:
        return 1                                  # overlap
    else:
        return 0


def GetDict(pos, num, dictp):
    # initial_vol = abs(np.linalg.det(latt))
    # get the atom type with smallest numbers
    l_dic = []
    for i, k in zip(dictp.keys(), dictp.values()):
        l_dic.append([i, k])
    #print('l_dic:', l_dic)
    srtl_dic = sorted(l_dic, key=lambda setx: setx[1])    # this is sorted list of all atoms number and site (important!)
    atypekey = srtl_dic[0][0]  # this is the atom type with smallest numbers
    #print('min_atom_type:\n', atypekey)
    #print('sorted dict:\n', srtl_dic)

    # get atomic positions
    #print('old pos: ', pos)
    pos = [list(i) for i in pos]
    #print('new pos: ', pos)
    dict1 = {}
    for i, k in enumerate(num):
        if k not in dict1.keys():
            dict1[k] = pos[i]
        else:
            dict1[k] = dict1[k] + pos[i]
    #print('old dict1:', dict1)
    for k in dict1.keys():
        count = 0
        tmp = []
        tp = []
        for j in dict1[k]:
            if count != 3:
                tp.append(j)
                count += 1
            elif count == 3:
                tmp.append(np.array(tp))
                count = 1
                tp = []
                tp.append(j)
        tmp.append(np.array(tp))
        dict1[k] = tmp
    #print('dict1', dict1)  # dict of {atom type, atom sites}
    '''
    dict11 = {}
    for k in srtl_dic:
        if k[0] not in dict11.keys():
            dict11[k[0]] = dict1[k[0]]
    print('new dict1', dict11)
    '''
    cel_size = len(num)     # // total numbers of atoms
    #print('atypekey: ', atypekey)
    #print('dict1: ', dict1)
    return atypekey, dict1, cel_size


# check the atom sites + trans
# set the first atom site with minimum number of atoms as origin
# get all difference between the atom sites with same type and the origin
def SymPureTrans2(is_found, trans, cel_size, dict1, atypekey, pos, latt, aty_i):
    tolerance = 0.00001
    num_tr = 0
    is_found_2 = is_found
    for i in range(cel_size):
        if not is_found_2[i]:       # // is_found = 0, continue next loop; otherwise(found), next step
            continue
        init_atom = i               # // find the equivalent atom of atom[i(init_atom)]

        for j in range(cel_size):
            vec = np.array(pos[init_atom]) + trans
            for k in range(cel_size):
                if k < aty_i:
                    continue
                if VecIsOverlap(vec, dict1[atypekey][k - aty_i], latt, tolerance):
                    if not is_found[k]:
                        is_found[k] = 1
                        num_tr += 1
                    init_atom = k       # // continue find the equivalent atom of atom[k]
                    break
            if init_atom == i:      # // do loop until finding all the equivalent atom[i]
                break
    return num_tr


# check total atom overlap
# choose the second number of atoms to check
def CheckTotalOverlap(rot, trans, dict1, latt):
    tolerance = 0.00001
    check_size = 0
    type_sorted = []
    if len(dict1.keys()) < 3:
        check_size = len(dict1.keys())
    else:
        check_size = 3

    for i in dict1.keys():
        type_sorted.append(i)   # // [Zr, O]  actually, it is num

    #print('type_sorted:\n', type_sorted)
    for i in range(check_size):
        type_name = type_sorted[i]  # // Zr, O
        #j1 = np.zeros(3)
        for j in dict1[type_name]:
            j1 = np.dot(rot, j) + trans     # // new pos
            tp_found = 0
            for k in dict1[type_name]:      #// for each j, find whether there is a k overlap or not
                if VecIsOverlap(j1, k, latt, tolerance):
                    tp_found = 1
                    break
            if tp_found == 0:
                return 0
    return 1


def SymPureTrans1(is_found, rot, ori, cel_size, dict1, atypekey, num, pos, latt):  # // it's better to use only dict1
    num_tr = 0

    aty_i = 0
    for k in dict1.keys():
        if k != atypekey:
            aty_i += len(dict1[k])
        elif k == atypekey:
            break

    #trans = np.zeros(3)
    for i in range(cel_size):
        if is_found[i]:
            continue
        if num[i] != atypekey:
            continue
        #print(dict1, atypekey, ori)
        #trans = dict1[atypekey][i] - ori
        trans = dict1[atypekey][i - aty_i] - ori
        # // add CheckTotalOverlap(trans, dict1)
        if CheckTotalOverlap(rot, trans, dict1, latt):
            is_found[i] = 1
            num_tr += 1
            if (rot == np.identity(3)).all():
                num_tr += SymPureTrans2(is_found, trans, cel_size, dict1, atypekey, pos, latt, aty_i)
            break
            '''
            vec = np.array(pos[i]) + trans
            for k in range(cel_size):
                if VecIsOverlap(vec, dict1[atypekey][k], latt, 0.00001):
                    is_found[i] = 1
                    num_tr += 1
                    num_tr += SymPureTrans2(is_found, trans, cel_size, dict1, atypekey, pos, latt)
                    break
            '''
    return num_tr


def SymTrans(rot, dictp, num, pos, latt):
    atypekey, dict1, cel_size = GetDict(pos, num, dictp)   # get the atom type with minimum numbers

    #num_tr = 0
    #ori = dict1[atypekey][0]
    ori = np.dot(rot, dict1[atypekey][0])
    vec_diff1 = []
    tr = []
    #print('origin site:\n', ori)
    is_found = []   # flag of whether find equivalent atom or not
    for i in range(cel_size):
        is_found.append(0)
    #print(is_found)

    num_tr = SymPureTrans1(is_found, rot, ori, cel_size, dict1, atypekey, num, pos, latt)
    #print(num_tr)
    if num_tr != 0:
        for i in range(cel_size):
            if is_found[i]:
                tp_tr = np.array(pos[i]) - ori
                tp_tr2 = tp_tr.copy()    # // Watch out! Don't use "tp_tr2 = tp_tr", otherwise they will change at the same time
                for k in range(3):
                    if tp_tr2[k] < 0.0:
                        tp_tr2[k] = int(tp_tr2[k] - 0.5)
                    else:
                        tp_tr2[k] = int(tp_tr2[k] + 0.5)
                tp_tr = tp_tr - tp_tr2
                for k in range(3):
                    if tp_tr[k] < 0.0:
                        tp_tr[k] += 1.0
                tr.append(tp_tr)

    return tr


def check_volume1(prim_vec, latt):
    initial_vol = abs(np.linalg.det(latt))
    vol_tolerance = 0.00001
    tmp_latt = np.zeros((3, 3))
    min_vec = np.zeros((3, 3))
    flag3 = 'Not Found'
    for i in range(0, len(prim_vec)-2):
        for j in range(i+1, len(prim_vec)-1):
            for k in range(j+1, len(prim_vec)):
                if flag3 == 'Not Found':
                    tmp_latt[0] = np.dot(prim_vec[i], latt)
                    tmp_latt[1] = np.dot(prim_vec[j], latt)
                    tmp_latt[2] = np.dot(prim_vec[k], latt)
                    tmp_vol = abs(np.linalg.det(tmp_latt))
                    #print('tmp_latt:\n', tmp_latt)
                    #print('tmpvol:\n', tmp_vol)
                    #print('initvol:\n', initial_vol)
                    #print('cardinality\n', cd.count(prim_vec))
                    if tmp_vol > vol_tolerance:
                        v = initial_vol / tmp_vol
                        if v < 0.0:
                            v = int(v - 0.5)
                        else:
                            v = int(v + 0.5)
                        if v == (cd.count(prim_vec) - 2):  #primitive satisfy, but convential one not satisify
                            min_vec[0] = prim_vec[i]       # coventional smallset is ~32. intial vol~127.8 if(127/32), the size~4
                            min_vec[1] = prim_vec[j]
                            min_vec[2] = prim_vec[k]
                            return check_volume2(min_vec, prim_vec, latt)   #flag3 = 'Found'
    return flag3, None  # flag3 = 'Not Found'



# 为了检查满足体积条件的prim_latt是不是本身就是input cell的lattice vec
#  如果我们找出来的prim_latt就是输入的，说明输入的本身就是个原胞
def check_volume2(min_vec, prim_vec, latt):
    tp_vec = np.transpose(min_vec)   # // tp_vec : vertical
    inv_tp_vec = np.linalg.inv(tp_vec)
    #print(inv_tp_vec)

    for i in range(3):
        for j in range(3):
            if inv_tp_vec[i, j] < 0.0:
                inv_tp_vec[i, j] = int(inv_tp_vec[i, j] - 0.5)
            else:
                inv_tp_vec[i, j] = int(inv_tp_vec[i, j] + 0.5)

    print('here:\n', np.linalg.det(inv_tp_vec))
    #print(inv_tp_vec)

    if abs(np.linalg.det(inv_tp_vec)) == (cd.count(prim_vec) - 2):
        print('Found')
    else:
        print('Warning! Not completely found primitive vectors!')

    # min_prim_latt = np.dot(latt, np.linalg.inv(inv_tp_vec))  # // tp_vec : vertical
    min_prim_latt = np.dot(np.transpose(np.linalg.inv(inv_tp_vec)), latt)  # // min_prim_latt : horizontal
    flag3 = 'Found'

    return flag3, min_prim_latt


# conventional changes to primitive
# find new atom position (for primitive)
def AtomPosConv2Prim(prim_latt, latt, pos, num):
    #print('check num: ', num)
    latt_vol = np.linalg.det(latt)
    prim_vol = np.linalg.det(prim_latt)
    ratio = latt_vol/prim_vol
    #print('ratio1:\n', ratio)
    cel_size1 = len(pos)
    #print('length of pos: ', len(pos))
    # if abs(ratio) >= 1:
    if ratio < 0.0:
        ratio = int(ratio - 0.5)
    else:
        ratio = int(ratio + 0.5)
    ratio = abs(ratio)
    #print('ratio2:\n', ratio)

    inv_tp_vec = np.linalg.inv(latt)
    P = np.transpose(np.dot(prim_latt, inv_tp_vec))  # P, Q is a vertical vec
    Q = np.linalg.inv(P)
    print('check Q1: ', np.linalg.det(Q))
    #deterQ = None
    #if abs(np.linalg.det(Q)) >= 1:
    for i in range(3):
        for j in range(3):
            if Q[i, j] < 0.0:
                Q[i, j] = int(Q[i, j] - 0.5)
            else:
                Q[i, j] = int(Q[i, j] + 0.5)
    #print('check Q2: ', np.linalg.det(Q))

    if abs(np.linalg.det(Q)) != ratio:  # P: big(conv) --> small(prim). so P is no larger than 1, Q is no smaller than 1
        print('Error1')
        return None, None, None

    if (cel_size1 / ratio) * ratio != cel_size1:
        print('Error2')
        return None, None, None

    tp_new_pos = delaunay.ChangeOfBasis(pos, Q)
    #print('temp_postion1:\n', tp_new_pos)

    new_site = []
    tolerance = 0.00001
    for i in range(cel_size1):
        new_site.append(i)     # // note: maybe index is out of range if it isn't initialized with len(cel_size)
    #print('new site: ', new_site)
    for cnt in range(100):
        for i in range(cel_size1):
            for j in range(cel_size1):
                if num[i] == num[j]:
                    if VecIsOverlap(np.array(tp_new_pos[i]), np.array(tp_new_pos[j]), prim_latt, tolerance):
                        if new_site[j] == j:
                            new_site[i] = j    # // change the equiv to the first showed one
                            break
    print('new site:\n', new_site)

    new_num = []
    for i in sorted(set(new_site)):
        new_num.append(num[i])  # // ['Zr', 'O', 'O']
    print('check!', new_num)

    tp_p = []
    for i in new_site:
        tp_p.append(tp_new_pos[i])      # // atom sites according to the new_site ['Zr', 'O', 'O']
    #print('temp_postion2:\n', tp_p)

    # // boundary trim
    new_posi = []
    for i in range(cel_size1):
        new_posi.append(np.zeros(3))

    for i in range(cel_size1):
        for j in range(3):
            if abs(tp_new_pos[i][j] - tp_p[i][j]) > 0.5:
                if tp_new_pos[i][j] < tp_p[i][j]:
                    new_posi[new_site[i]][j] += tp_new_pos[i][j] + 1
                else:
                    new_posi[new_site[i]][j] += tp_new_pos[i][j] - 1
            else:
                new_posi[new_site[i]][j] = tp_new_pos[i][j]

        # // average of the overlapped two atoms
        ave_pos = cel_size1 / len(new_posi)
        for j in range(len(new_posi)):
            #tp2 = np.zeros(3)
            for k in range(3):
                new_posi[j][k] = new_posi[j][k] / ave_pos
            tp2 = new_posi[j].copy()    # // Watch out! Don't use "tp2 = new_posi[j]", otherwise they will change at the same time
            for k in range(3):
                if tp2[k] < 0.0:
                    tp2[k] = int(tp2[k] - 0.5)
                else:
                    tp2[k] = int(tp2[k] + 0.5)
            new_posi[j] = new_posi[j] - tp2
            for k in range(3):
                if new_posi[j][k] < 0.0:
                    new_posi[j][k] += 1.0

    new_pos = []
    for i in sorted(set(new_site)):
        new_pos.append(new_posi[i])

    #print('new numbers:\n', new_num)
    #print('new positions:\n', [list(i) for i in new_pos])

    return prim_latt, new_num, new_pos



# primitive changes to conventional 
def AtomPosPrim2Conv(prim_latt, latt, pos, num, cel_size):
    latt_vol = np.linalg.det(latt)
    prim_vol = np.linalg.det(prim_latt)
    ratio = latt_vol/prim_vol
    print('ratio1:\n', ratio)
    #cel_size = len(pos)

    ratio = np.around(ratio, decimals=3)
    ratio = abs(ratio)
    print('ratio2:\n', ratio)

    inv_tp_vec = np.linalg.inv(latt)
    P = np.transpose(np.dot(prim_latt, inv_tp_vec))  # P, Q is a vertical vec
    Q = np.linalg.inv(P)
    print('check Q1: ', np.linalg.det(Q))

    deterQ = np.around(np.linalg.det(Q), decimals=3)
    print('check Q2: ', deterQ)

    if abs(deterQ) != ratio:  # P: big(conv) --> small(prim). so P is no larger than 1, Q is no smaller than 1
        print('Error1')
        return None, None, None

    if (cel_size / ratio) * ratio != cel_size:
        print('Error2')
        return None, None, None

    tp_new_pos = delaunay.ChangeOfBasis(pos, Q)
    #print('temp_postion1:\n', tp_new_pos)

    new_site = []
    tolerance = 0.00001
    for i in range(cel_size):
        new_site.append(i)     # // note: maybe index is out of range if it isn't intialized with len(cel_size)
    for cnt in range(100):
        for i in range(cel_size):
            for j in range(cel_size):
                if num[i] == num[j]:
                    if VecIsOverlap(np.array(tp_new_pos[i]), np.array(tp_new_pos[j]), prim_latt, tolerance):
                        if new_site[j] == j:
                            new_site[i] = j    # // change the equiv to the first showed one
                            break
    #print('new site:\n', new_site)

    new_num = []
    for i in set(new_site):
        new_num.append(num[i])  # // ['Zr', 'O', 'O']
    print('check!', new_num)

    tp_p = []
    for i in new_site:
        tp_p.append(tp_new_pos[i])      # // atom sites according to the new_site ['Zr', 'O', 'O']
    #print('temp_postion2:\n', tp_p)

    # // boundary trim
    new_posi = []
    for i in range(cel_size):
        new_posi.append(np.zeros(3))

    for i in range(cel_size):
        for j in range(3):
            if abs(tp_new_pos[i][j] - tp_p[i][j]) > 0.5:
                if tp_new_pos[i][j] < tp_p[i][j]:
                    new_posi[new_site[i]][j] += tp_new_pos[i][j] + 1
                else:
                    new_posi[new_site[i]][j] += tp_new_pos[i][j] - 1
            else:
                new_posi[new_site[i]][j] = tp_new_pos[i][j]

        # // average of the overlapped two atoms
        ave_pos = cel_size / len(new_posi)
        for j in range(len(new_posi)):
            #tp2 = np.zeros(3)
            for k in range(3):
                new_posi[j][k] = new_posi[j][k] / ave_pos
            tp2 = new_posi[j].copy()    # // Watch out! Don't use "tp2 = new_posi[j]", otherwise they will change at the same time
            for k in range(3):
                if tp2[k] < 0.0:
                    tp2[k] = int(tp2[k] - 0.5)
                else:
                    tp2[k] = int(tp2[k] + 0.5)
            new_posi[j] = new_posi[j] - tp2
            for k in range(3):
                if new_posi[j][k] < 0.0:
                    new_posi[j][k] += 1.0

    new_pos = []
    for i in set(new_site):
        new_pos.append(new_posi[i])

    print('new numbers:\n', new_num)
    print('new positions:\n', [list(i) for i in new_pos])

    return prim_latt, new_num, new_pos



def GetPrimLattice(latt, pos, num, dictp):
    is_prim = None
    rot = np.identity(3)
    tr = SymTrans(rot, dictp, num, pos, latt)
    #print('tr:\n', tr)
    if tr is not None:
        if len(tr) == 1:
            print('Input cell is primitive!')
            is_prim = 1  # // primitive cell flag
            #print('primitive lattice:\n', latt)  # // delaunay of primitive lattice vec
            #print('primitive numbers:\n', num)
            #print('primitive cell positions:\n', pos)  # atom postion
            #return latt, num, pos
            #// delaunay change part
            flag, reduc_b, delauP = delaunay.Delaunay(latt, -1)
            if flag:
                new_pos = delaunay.ChangeOfBasis(pos, np.linalg.inv(delauP))
                print('primitive lattice (delaunay):\n', reduc_b)      # // delaunay of primitive lattice vec
                print('primitive cell positions:\n', [list(i) for i in new_pos])  # atom postion
                return reduc_b, num, new_pos, is_prim
        else:
            print('Input cell is not primitive! Changing...')
            is_prim = 0  # // convetional cell flag
            ori_tr = [np.array([1.000000000, 0.0000000000, 0.0000000000]), np.array([0.000000000, 1.0000000000, 0.0000000000]), np.array([0.000000000, 0.0000000000, 1.0000000000])]
            #print('original:\n', ori_tr)
            prim_vec = tr[1:] + ori_tr      # delete the first vec in tr (ori-ori)
            print('all primitive translations:\n', prim_vec)
            #min_prim_latt = np.zeros((3, 3))
            flag3, min_prim_latt = check_volume1(prim_vec, latt)         # // delaunay of primitive lattice vec
            if flag3 == 'Found':
                #print('final primitive lattice\n', min_prim_latt)
                #prim_latt, new_num, new_pos = AtomPosConv2Prim(min_prim_latt, latt, pos, num)
                #return prim_latt, new_num, new_pos
                #print('new numbers:\n', new_num)
                #print('new positions:\n', new_pos)
                # // delaunay change part
                flag, reduc_b, delaup = delaunay.Delaunay(min_prim_latt, -1)
                if flag:
                    print('primitive lattice (delaunay):\n', reduc_b)
                    reduc_latt, new_num, new_pos = AtomPosConv2Prim(reduc_b, latt, pos, num)
                    print('primitive cell positions:\n', [list(i) for i in new_pos])  # atom postion
                    return reduc_latt, new_num, new_pos, is_prim



def test_delaunay():
    flag, reduc_b, delauP = Delaunay(latt, -1)
    print(reduc_b)
    print(delauP)
