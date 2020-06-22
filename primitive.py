#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Created on 2019-03-06
Update  on 2020-06-22
Author: Cecilia9999
GitHub: https://github.com/Cecilia9999/
"""

'''
    This file is for finding primitive cell
'''

import numpy as np
import cardinality as cd
import delaunay

 
def VecIsOverlap(vec1, vec2, latt, tolerance=0.00001):
    """
        Determine if two vectors are equivalent or not
    
        input  :    vec1, vec2  (np.array)
                    latt        (3x3 np.array)
                    tolerance   (default = 0.00001)
        return :    True / False
    """
    #if tolerance is None:   
    #    tolerance = 0.00001
        
    vec_diff = vec1 - vec2
    tp_diff = vec_diff.copy()
    
    for i in range(3):
        if tp_diff[i] < 0.0:
            tp_diff[i] = int(tp_diff[i] - 0.5)
        else:
            tp_diff[i] = int(tp_diff[i] + 0.5)
    
    nvec_diff = vec_diff - tp_diff
    
    if np.linalg.norm(np.dot(np.transpose(latt), nvec_diff)) < tolerance:
        return 1                                  
    
    return 0


def GetDict(pos, num, dictp):
    """ 
        Get atomic information from input
    
        input :     POSCAR
        return:     atypekey: atom type (int)
                    dict    : {atom type, atom sites}
                    
    """
    # initial_vol = abs(np.linalg.det(latt))
    # get the atom type with smallest numbers
    l_dic = []
    for i, k in zip(dictp.keys(), dictp.values()):
        l_dic.append([i, k])

    # sorted list of all atoms number and site 
    srtl_dic = sorted(l_dic, key=lambda setx: setx[1])    
    
    # atom type with smallest numbers
    atypekey = srtl_dic[0][0]  
    
    # get atomic positions
    pos = [list(i) for i in pos]

    dict1 = {}
    for i, k in enumerate(num):
        if k not in dict1.keys():
            dict1[k] = pos[i]
        else:
            dict1[k] = dict1[k] + pos[i]

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
        
        # dict : {atom type, atom sites}
        dict1[k] = tmp      
        
    # total numbers of atoms
    cel_size = len(num)     
    
    return atypekey, dict1, cel_size



def SymPureTrans2(is_found, trans, cel_size, dict1, atypekey, pos, latt, aty_i):
    """ 
        check the atom sites + trans
        set the first atom site with minimum number of atoms as origin
        get all difference between the atom sites with same type and the origin
    """
    tolerance = 0.00001
    num_tr = 0
    is_found_2 = is_found
    for i in range(cel_size):
        
        # is_found = 0, continue next loop; else (found), next step
        if not is_found_2[i]:       
            continue
         # find the equivalent atom of atom[i(init_atom)]
        init_atom = i              

        for j in range(cel_size):
            vec = np.array(pos[init_atom]) + trans
            for k in range(cel_size):
                if k < aty_i:
                    continue
                if VecIsOverlap(vec, dict1[atypekey][k - aty_i], latt, tolerance):
                    if not is_found[k]:
                        is_found[k] = 1
                        num_tr += 1
                    
                    # continue find the equivalent atom of atom[k]
                    init_atom = k       
                    break
                    
            # loop until finding all the equivalent atom[i]
            if init_atom == i:      
                break
                
    return num_tr



def CheckTotalOverlap(rot, trans, dict1, latt):
    """
        check total atom overlap
        choose the second number of atoms to check
    """
    tolerance = 0.00001
    check_size = 0
    type_sorted = []
    if len(dict1.keys()) < 3:
        check_size = len(dict1.keys())
    else:
        check_size = 3
    
    # test: [Zr, O] corresponding number 
    for i in dict1.keys():  
        type_sorted.append(i)                         

    #print('type_sorted:\n', type_sorted)
    for i in range(check_size):
        type_name = type_sorted[i]   # test: Zr, O
        
        for j in dict1[type_name]:
            
            # new position np.zeros(3)
            j1 = np.dot(rot, j) + trans     
            tp_found = 0
            
            # for each j, find if there is a k overlap or not
            for k in dict1[type_name]:      
                if VecIsOverlap(j1, k, latt, tolerance):
                    tp_found = 1
                    break
                    
            if tp_found == 0:
                return 0
                
    return 1



def SymPureTrans1(is_found, rot, ori, cel_size, dict1, atypekey, num, pos, latt):  
    """
        should update (parameters in func are too many)
    """
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
        
        trans = dict1[atypekey][i - aty_i] - ori
        
        # add CheckTotalOverlap(trans, dict1)
        if CheckTotalOverlap(rot, trans, dict1, latt):
            is_found[i] = 1
            num_tr += 1
            if (rot == np.identity(3)).all():
                num_tr += SymPureTrans2(is_found, trans, cel_size, dict1, atypekey, pos, latt, aty_i)
            
            break
            
            ''' OLD TYPE
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
    # get the atom type with minimum numbers
    atypekey, dict1, cel_size = GetDict(pos, num, dictp)   

    #num_tr = 0
    #ori = dict1[atypekey][0]
    ori = np.dot(rot, dict1[atypekey][0])
    vec_diff1 = []
    tr = []
    #print('origin site:\n', ori)
   
    # flag of whether find equivalent atom or not
    is_found = []   
    for i in range(cel_size):
        is_found.append(0)
    #print(is_found)
    
    # call SymPureTrans1 firstly
    num_tr = SymPureTrans1(is_found, rot, ori, cel_size, dict1, atypekey, num, pos, latt)

    if num_tr != 0:
        for i in range(cel_size):
            if is_found[i]:
                tp_tr = np.array(pos[i]) - ori
                
                # Note: Don't use "tp_tr2 = tp_tr"
                tp_tr2 = tp_tr.copy()   
                
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
                            
                        # only primitive cell will satisfy, but not conventional one
                        # coventional smallset is ~32. initial vol~127.8 if(127/32), the size~4
                        if v == (cd.count(prim_vec) - 2):  
                            min_vec[0] = prim_vec[i]       
                            min_vec[1] = prim_vec[j]
                            min_vec[2] = prim_vec[k]
                            return check_volume2(min_vec, prim_vec, latt)   # flag3 = 'Found'
    
    return flag3, None  # flag3 = 'Not Found'



# 为了检查满足体积条件的 primitive lattice 是不是本身就是 input cell 的 lattice vector
# 如果我们找出来的 primitive lattice 就是输入的，说明输入的本身就是个原胞
def check_volume2(min_vec, prim_vec, latt):
    """
        To check if input cell is primitive or not
        input :  min_vec  :   minimum vetors of an atom
                 prim_vec :   basic vectors of primivtive cell
                 latt     :   lattice vectors
        return:  flage3   :   True: input cell is primitive
                 min_prim_latt: 
    """
    tp_vec = np.transpose(min_vec)      # tp_vec : vertical
    inv_tp_vec = np.linalg.inv(tp_vec)
    #print(inv_tp_vec)

    for i in range(3):
        for j in range(3):
            if inv_tp_vec[i, j] < 0.0:
                inv_tp_vec[i, j] = int(inv_tp_vec[i, j] - 0.5)
            else:
                inv_tp_vec[i, j] = int(inv_tp_vec[i, j] + 0.5)

    print('inv_tp_vec:\n', np.linalg.det(inv_tp_vec))
    #print(inv_tp_vec)

    if abs(np.linalg.det(inv_tp_vec)) == (cd.count(prim_vec) - 2):
        print('Found')
    else:
        print('Warning! Not completely found primitive vectors!')

    # min_prim_latt = np.dot(latt, np.linalg.inv(inv_tp_vec))
    # inv_tp_vec : vertical
    # min_prim_latt : horizontal
    min_prim_latt = np.dot(np.transpose(np.linalg.inv(inv_tp_vec)), latt)  
    flag3 = 'Found'

    return flag3, min_prim_latt



def AtomPosConv2Prim(prim_latt, latt, pos, num):
    """
        Conventional lattice changes to primitive lattice
        Find correspoing new atomic positions (for primitive)
    """

    latt_vol = np.linalg.det(latt)
    prim_vol = np.linalg.det(prim_latt)
    ratio = latt_vol/prim_vol

    cel_size1 = len(pos)

    if ratio < 0.0:
        ratio = int(ratio - 0.5)
    else:
        ratio = int(ratio + 0.5)
    
    ratio = abs(ratio)
    #print('ratio2:\n', ratio)

    inv_tp_vec = np.linalg.inv(latt)
    P = np.transpose(np.dot(prim_latt, inv_tp_vec))  # P, Q : vertical vectors
    Q = np.linalg.inv(P)
    # print('check Q1: ', np.linalg.det(Q))

    for i in range(3):
        for j in range(3):
            if Q[i, j] < 0.0:
                Q[i, j] = int(Q[i, j] - 0.5)
            else:
                Q[i, j] = int(Q[i, j] + 0.5)
    # print('check Q2: ', np.linalg.det(Q))
    
    # P: big(conv) --> small(prim)
    # P is no larger than 1, but Q is no smaller than 1
    if abs(np.linalg.det(Q)) != ratio: 
        print('Error1')
        return None, None, None

    if (cel_size1 / ratio) * ratio != cel_size1:
        print('Error2')
        return None, None, None

    tp_new_pos = delaunay.ChangeOfBasis(pos, Q)
    #print('temp_postion1:\n', tp_new_pos)

    new_site = []
    tolerance = 0.00001
    
    # note: index maybe out of range if it isn't initialized with len(cel_size)
    for i in range(cel_size1):
        new_site.append(i)     
    
    for cnt in range(100):
        for i in range(cel_size1):
            for j in range(cel_size1):
                if num[i] == num[j]:
                    if VecIsOverlap(np.array(tp_new_pos[i]), np.array(tp_new_pos[j]), prim_latt, tolerance):
                        # change the equivalent to the first showed one
                        if new_site[j] == j:
                            new_site[i] = j    
                            break
    # print('new site:\n', new_site)

    new_num = []
    for i in sorted(set(new_site)):
        new_num.append(num[i])          # test: ['Zr', 'O', 'O']
    
    # print('check!', new_num)

    tp_p = []
    for i in new_site:
        tp_p.append(tp_new_pos[i])      # test: atomic sites according to the new_site ['Zr', 'O', 'O']
  
    # boundary trim
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

        # average of the overlapped two atoms
        ave_pos = cel_size1 / len(new_posi)
        for j in range(len(new_posi)):
            #tp2 = np.zeros(3)
            for k in range(3):
                new_posi[j][k] = new_posi[j][k] / ave_pos
                
            # Note: Don't use "tp2 = new_posi[j]" (change at the same time)
            tp2 = new_posi[j].copy()   
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

    print('new numbers:\n', new_num)
    print('new positions:\n', [list(i) for i in new_pos])

    return prim_latt, new_num, new_pos



def AtomPosPrim2Conv(prim_latt, latt, pos, num, cel_size):
    """
        Primitive lattice changes to conventional lattice 
    """
    latt_vol = np.linalg.det(latt)
    prim_vol = np.linalg.det(prim_latt)
    ratio = latt_vol/prim_vol
    # print('ratio1:\n', ratio)

    ratio = np.around(ratio, decimals=3)
    ratio = abs(ratio)
    # print('ratio2:\n', ratio)

    inv_tp_vec = np.linalg.inv(latt)
    P = np.transpose(np.dot(prim_latt, inv_tp_vec))  # P, Q is a vertical vec
    Q = np.linalg.inv(P)
    # print('check Q1: ', np.linalg.det(Q))

    deterQ = np.around(np.linalg.det(Q), decimals=3)
    # print('check Q2: ', deterQ)
    
    # P: big(conv) --> small(prim)
    # P is no larger than 1, Q is no smaller than 1
    if abs(deterQ) != ratio:  
        print('Error1')
        return None, None, None

    if (cel_size / ratio) * ratio != cel_size:
        print('Error2')
        return None, None, None

    tp_new_pos = delaunay.ChangeOfBasis(pos, Q)
    # print('temp_postion1:\n', tp_new_pos)

    new_site = []
    tolerance = 0.00001
    
    # note: index maybe out of range if it isn't intialized with len(cel_size)
    for i in range(cel_size):
        new_site.append(i)     
        
    for cnt in range(100):
        for i in range(cel_size):
            for j in range(cel_size):
                if num[i] == num[j]:
                    if VecIsOverlap(np.array(tp_new_pos[i]), np.array(tp_new_pos[j]), prim_latt, tolerance):
                        
                        # change the equiv to the first showed one
                        if new_site[j] == j:
                            new_site[i] = j    
                            break
    #print('new site:\n', new_site)

    new_num = []
    for i in set(new_site):
        new_num.append(num[i])          # test: ['Zr', 'O', 'O']
    # print('check!', new_num)

    tp_p = []
    for i in new_site:
        tp_p.append(tp_new_pos[i])      # test: atom sites according to the new_site ['Zr', 'O', 'O']
    # print('temp_postion2:\n', tp_p)

    # boundary trim
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

        # average of the overlapped two atoms
        ave_pos = cel_size / len(new_posi)
        for j in range(len(new_posi)):
            #tp2 = np.zeros(3)
            for k in range(3):
                new_posi[j][k] = new_posi[j][k] / ave_pos
            
            # Don't use "tp2 = new_posi[j]"
            tp2 = new_posi[j].copy() 
               
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
            
            # Primitive cell flag
            is_prim = 1 
            # print('primitive lattice:\n', latt)  # delaunay of primitive lattice vec
            # print('primitive numbers:\n', num)
            # print('primitive cell positions:\n', pos)  # atom postion
            # return latt, num, pos
           
            # Delaunay Transformation
            flag, reduc_b, delauP = delaunay.Delaunay(latt, -1)
            if flag:
                new_pos = delaunay.ChangeOfBasis(pos, np.linalg.inv(delauP))
                
                # delaunay of primitive lattice vec
                print('primitive lattice (delaunay):\n', reduc_b)     
                
                # atomic postions
                print('primitive cell positions:\n', [list(i) for i in new_pos])  
                
                return reduc_b, num, new_pos, is_prim
        else:
            print('Input cell is not primitive! Changing...')
            
            # Convetional cell flag
            is_prim = 0  
            ori_tr = [np.array([1.000000000, 0.0000000000, 0.0000000000]), \
                      np.array([0.000000000, 1.0000000000, 0.0000000000]), \
                      np.array([0.000000000, 0.0000000000, 1.0000000000])]
            
            # delete the first vec in tr (ori-ori)
            prim_vec = tr[1:] + ori_tr      
            print('all primitive translations:\n', prim_vec)
            
            # delaunay of primitive lattice vec
            # min_prim_latt : np.zeros((3, 3))
            flag3, min_prim_latt = check_volume1(prim_vec, latt)         
            if flag3 == 'Found':
                # print('final primitive lattice\n', min_prim_latt)
                # prim_latt, new_num, new_pos = AtomPosConv2Prim(min_prim_latt, latt, pos, num)
                # return prim_latt, new_num, new_pos
                # print('new numbers:\n', new_num)
                # print('new positions:\n', new_pos)
                
                # Delaunay Transformation
                flag, reduc_b, delaup = delaunay.Delaunay(min_prim_latt, -1)
                if flag:
                    
                    # delaunay of primitive lattice vec
                    print('primitive lattice (delaunay):\n', reduc_b)
                    
                    # atomic postions
                    reduc_latt, new_num, new_pos = AtomPosConv2Prim(reduc_b, latt, pos, num)
                    print('primitive cell positions:\n', [list(i) for i in new_pos])  
                    
                    return reduc_latt, new_num, new_pos, is_prim


def test_delaunay():
    flag, reduc_b, delauP = Delaunay(latt, -1)
    print(reduc_b)
    print(delauP)
