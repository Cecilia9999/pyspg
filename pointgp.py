#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Created on 2019-03-10
Update  on 2020-06-23
Author: Cecilia9999
GitHub: https://github.com/Cecilia9999/
"""

'''
    This file is for finding point group
'''

import numpy as np
import primitive
import delaunay


def AllRotOperation():
    """
        Find all point operations
    """
    rot = [[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1],
           [-1, 0, 0],
           [0, -1, 0],      # 5 
           [0, 0, -1],
           [0, 1, 1],
           [1, 0, 1],
           [1, 1, 0],
           [0, -1, -1],     # 10
           [-1, 0, -1],
           [-1, -1, 0],
           [0, 1, -1],
           [-1, 0, 1],
           [1, -1, 0],      # 15
           [0, -1, 1],
           [1, 0, -1],
           [-1, 1, 0],
           [1, 1, 1],
           [-1, -1, -1],    # 20
           [-1, 1, 1],
           [1, -1, 1],
           [1, 1, -1],
           [1, -1, -1],
           [-1, 1, -1],     # 25 
           [-1, -1, 1],
           ]

    rotations = []

    for i in range(26):
        for j in range(26):
            for k in range(26):
                tp_rot = np.zeros((3, 3))
                tp_rot[0] = np.array(rot[i])
                tp_rot[1] = np.array(rot[j])
                tp_rot[2] = np.array(rot[k])
                if abs(np.linalg.det(tp_rot)) == 1:
                    rotations.append(tp_rot)
    
    # print('len rot', len(rotations))
    return rotations


def GetMetricMat(latt):
    """
        Get metric of lattice
    """
    metric = np.dot(latt, np.transpose(latt))
    return metric


def MetricCondition(metric, metric_rot):
    """
        Metric Matrix should satisfy G conditions
    """
    tolerance = 0.00001
    
    # Length part
    len_met = np.zeros(3)        
    len_met_rot = np.zeros(3)
    
    for i in range(3):
        len_met[i] = np.sqrt(metric[i, i])
        len_met_rot[i] = np.sqrt(metric_rot[i, i])
        # print(len_met[i], len_met_rot[i])
        
        if abs(len_met_rot[i] - len_met[i]) > tolerance:
            return 0

    # Angle part
    tp_set = ((0, 1), (0, 2), (1, 2))
    for s in tp_set:
        i = s[0]
        j = s[1]
        dangel = np.arccos(metric_rot[i, j] / len_met_rot[i] / len_met_rot[j]) - np.arccos(metric[i, j] / len_met[i] / len_met[j])
        # print(dangel)
        len_ave = (len_met_rot[i] + len_met[i]) * (len_met_rot[j] + len_met[j]) / 4.0
        # print('1:', np.square(np.sin(dangel)) * len_ave)
        # print('sin theta: ', np.square(np.sin(dangel)))
        # print('s2: ', np.square(np.sin(dangel)) * len_ave)
        
        # Not ensure whether, after delaunay changing, this condition should be added or not
        if np.square(np.sin(dangel)) > 0.000000000001:   
            if (np.square(np.sin(dangel)) * len_ave) > 0.0000000001:    # tolerance * tolerance:
                return 0
                
    # we can also add receive an angle tolerance, and use it to check directly
    return 1
    
 
def GetRotOfLowerSymm(new_latt, old_latt, old_rot):
    """
        get rotation of primitive lattice
        rot of prim have lower symmetry than that of delaunay
            W: rot of old latt    
            W': rot of new latt    
            Q: inv(new latt) * old latt
            W' = QWP
    """
    
    # Q = np.transpose(np.dot(old_latt, np.linalg.inv(new_latt)))
    P = np.dot(new_latt, np.linalg.inv(old_latt))   # // horizontal vec
    tolerance = abs(np.linalg.det(P)) / 10

    new_rot = []
    for i in old_rot:
        # tp_rot = np.dot(P, i)
        # tp_rot1 = np.dot(tp_rot, np.linalg.inv(P))  
        # new_rot W'=QWP
        tp_rot = np.dot(np.linalg.inv(np.transpose(P)), i)
        
        # standard horizontal
        tp_rot1 = np.dot(tp_rot, np.transpose(P))     
        tp_rot2 = tp_rot1.copy()
        
        # find matrix whose elements are all integer
        for j in range(3):
            for k in range(3):
                if tp_rot2[j, k] < 0.0:
                    tp_rot2[j, k] = int(tp_rot2[j, k] - 0.5)
                else:
                    tp_rot2[j, k] = int(tp_rot2[j, k] + 0.5)
        # print('1:', tp_rot1)
        # print('2:', tp_rot2)
        
        tp_flag = True
        for j in range(3):
            for k in range(3):
                if abs(tp_rot2[j, k] - tp_rot1[j, k]) > tolerance:
                    tp_flag = False
                    break
            if not tp_flag:
                break

        if tp_flag:
            if abs(np.linalg.det(tp_rot2)) == 1:    # tp_rot1 or tp_rot2
                new_rot.append(tp_rot2)
    
    # print('check new len: ', len(new_rot))
    
    if len(new_rot) == 0:
        print('Error! Cnnot find any new rotations of primitive lattice!')
        print('Return old rotations of minimum lattice...')
        return old_rot
    else:
        return new_rot



def SymRotation(prim_latt):
    """
        Get point group, G before rotation  
    """
    # min_latt, num, pos = primitive.GetPrimLattice()
    # min_latt should be the delaunay lattice, that is "delaunay.Delaunay(prim_latt)"
    flag, min_latt, delauP = delaunay.Delaunay(prim_latt, -1)
    
    # min_latt = prim_latt.copy()
    metric = GetMetricMat(min_latt)
    
    # G after rotation
    rotations = AllRotOperation()
    latt_oper = []
    new_rot = []
    count = 0
    for i in rotations:

        rot_latt = np.dot(np.transpose(i), min_latt)
        metric_rot = GetMetricMat(rot_latt)
        if MetricCondition(metric, metric_rot):
            if count < 48:
                latt_oper.append(i)
                count += 1
                # print('check here count: ', count)
            
            elif count == 48:
                print('Too many symmetric operations! Stop...')
                break
    # print('count', count)
    
    if count < 49:
        new_rot = GetRotOfLowerSymm(prim_latt, min_latt, latt_oper)

    # print(len(latt_oper), len(new_rot))
    # print(len(new_rot))  
    # All symmetry rotation parts
    return new_rot


def GetRotTable(allrot):
    """
        Get point goup type via LookUpTable
    """
    # { "det": [trace, operation] }
    LookUpTable1 = {            
        "-1": [[-2, -6],
               [-1, -4],
               [0, -3],
               [1, -2],
               [-3, -1]],
        "1": [[3, 1],
              [-1, 2],
              [0, 3],
              [1, 4],
              [2, 6]]
    }

    # identification of point group
    LookUpTable2 = [
        [  # /* 0 */
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            "     ",
            "   ",
            "HOLOHEDRY_NONE",
            "LAUE_NONE"
        ],
        [  # /* 1 */
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            "1    ",
            "C1 ",
            "TRICLI",
            "LAUE1"
        ],
        [  # /* 2 */
            [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
            "-1   ",
            "Ci ",
            "TRICLI",
            "LAUE1"
        ],
        [  # /* 3 */
            [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
            "2    ",
            "C2 ",
            "MONOCLI",
            "LAUE2M"
        ],
        [  # /* 4 */
            [0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
            "m    ",
            "Cs ",
            "MONOCLI",
            "LAUE2M"
        ],
        [  # /* 5 */
            [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
            "2/m  ",
            "C2h",
            "MONOCLI",
            "LAUE2M"
        ],
        [  # /* 6 */
            [0, 0, 0, 0, 0, 1, 3, 0, 0, 0],
            "222  ",
            "D2 ",
            "ORTHO",
            "LAUEMMM"
        ],
        [  # /* 7 */
            [0, 0, 0, 2, 0, 1, 1, 0, 0, 0],
            "mm2  ",
            "C2v",
            "ORTHO",
            "LAUEMMM"
        ],
        [  # /* 8 */
            [0, 0, 0, 3, 1, 1, 3, 0, 0, 0],
            "mmm  ",
            "D2h",
            "ORTHO",
            "LAUEMMM"
        ],
        [  # /* 9 */
            [0, 0, 0, 0, 0, 1, 1, 0, 2, 0],
            "4    ",
            "C4 ",
            "TETRA",
            "LAUE4M"
        ],
        [  # /* 10 */
            [0, 2, 0, 0, 0, 1, 1, 0, 0, 0],
            "-4   ",
            "S4 ",
            "TETRA",
            "LAUE4M"
        ],
        [  # /* 11 */
            [0, 2, 0, 1, 1, 1, 1, 0, 2, 0],
            "4/m  ",
            "C4h",
            "TETRA"
            "LAUE4M"
        ],
        [  # /* 12 */
            [0, 0, 0, 0, 0, 1, 5, 0, 2, 0],
            "422  ",
            "D4 ",
            "TETRA",
            "LAUE4MMM"
        ],
        [  # /* 13 */
            [0, 0, 0, 4, 0, 1, 1, 0, 2, 0],
            "4mm  ",
            "C4v",
            "TETRA",
            "LAUE4MMM"
        ],
        [  # /* 14 */
            [0, 2, 0, 2, 0, 1, 3, 0, 0, 0],
            "-42m ",
            "D2d",
            "TETRA",
            "LAUE4MMM",
        ],
        [  # /* 15 */
            [0, 2, 0, 5, 1, 1, 5, 0, 2, 0],
            "4/mmm",
            "D4h",
            "TETRA",
            "LAUE4MMM"
        ],
        [  # /* 16 */
            [0, 0, 0, 0, 0, 1, 0, 2, 0, 0],
            "3    ",
            "C3 ",
            "TRIGO",
            "LAUE3"
        ],
        [  # /* 17 */
            [0, 0, 2, 0, 1, 1, 0, 2, 0, 0],
            "-3   ",
            "C3i",
            "TRIGO",
            "LAUE3"
        ],
        [  # /* 18 */
            [0, 0, 0, 0, 0, 1, 3, 2, 0, 0],
            "32   ",
            "D3 ",
            "TRIGO",
            "LAUE3M"
        ],
        [  # /* 19 */
            [0, 0, 0, 3, 0, 1, 0, 2, 0, 0],
            "3m   ",
            "C3v",
            "TRIGO",
            "LAUE3M"
        ],
        [  # /* 20 */
            [0, 0, 2, 3, 1, 1, 3, 2, 0, 0],
            "-3m  ",
            "D3d",
            "TRIGO",
            "LAUE3M"
        ],
        [  # /* 21 */
            [0, 0, 0, 0, 0, 1, 1, 2, 0, 2],
            "6    ",
            "C6 ",
            "HEXA",
            "LAUE6M"
        ],
        [  # /* 22 */
            [2, 0, 0, 1, 0, 1, 0, 2, 0, 0],
            "-6   ",
            "C3h",
            "HEXA",
            "LAUE6M"
        ],
        [  # /* 23 */
            [2, 0, 2, 1, 1, 1, 1, 2, 0, 2],
            "6/m  ",
            "C6h",
            "HEXA",
            "LAUE6M"
        ],
        [  # /* 24 */
            [0, 0, 0, 0, 0, 1, 7, 2, 0, 2],
            "622  ",
            "D6 ",
            "HEXA",
            "LAUE6MMM"
        ],
        [  # /* 25 */
            [0, 0, 0, 6, 0, 1, 1, 2, 0, 2],
            "6mm  ",
            "C6v",
            "HEXA",
            "LAUE6MMM"
        ],
        [  # /* 26 */
            [2, 0, 0, 4, 0, 1, 3, 2, 0, 0],
            "-6m2 ",
            "D3h",
            "HEXA",
            "LAUE6MMM"
        ],
        [  # /* 27 */
            [2, 0, 2, 7, 1, 1, 7, 2, 0, 2],
            "6/mmm",
            "D6h",
            "HEXA",
            "LAUE6MMM"
        ],
        [  # /* 28 */
            [0, 0, 0, 0, 0, 1, 3, 8, 0, 0],
            "23   ",
            "T  ",
            "CUBIC",
            "LAUEM3"
        ],
        [  # /* 29 */
            [0, 0, 8, 3, 1, 1, 3, 8, 0, 0],
            "m-3  ",
            "Th ",
            "CUBIC",
            "LAUEM3"
        ],
        [  # /* 30 */
            [0, 0, 0, 0, 0, 1, 9, 8, 6, 0],
            "432  ",
            "O  ",
            "CUBIC",
            "LAUEM3M"
        ],
        [  # /* 31 */
            [0, 6, 0, 6, 0, 1, 3, 8, 0, 0],
            "-43m ",
            "Td ",
            "CUBIC",
            "LAUEM3M"
        ],
        [  # /* 32 */
            [0, 6, 8, 9, 1, 1, 9, 8, 6, 0],
            "m-3m ",
            "Oh ",
            "CUBIC",
            "LAUEM3M"
        ]
    ]

    # table = {'-6': 0, '-4': 0, '-3': 0, '-2': 0, '-1': 0, '1': 0, '2': 0, '3': 0, '4': 0, '6': 0}
    table = {}
    for m in LookUpTable1.values():
        for j in m:
            table[str(j[1])] = 0

    for i in allrot:
        tp_det = np.linalg.det(i)
        tp_tra = np.trace(i)  # // trace of matrix
        for tr in LookUpTable1[str(int(tp_det))]:
            if tr[0] == int(tp_tra):
                table[str(tr[1])] += 1  # table = {'-6': 1, '-4': 0, ....}
    # print(table)

    pgnum = 0
    pgtable = None
    symbol = None
    schoenf = None
    holo = None
    Laue = None
    for count, i in enumerate(LookUpTable2):
        if i[0] == list(table.values()):
            pgnum = count
            pgtable = i[0]
            symbol = i[1]
            schoenf = i[2]
            holo = i[3]
            Laue = i[4]

    pgt_type = [pgnum, pgtable, symbol, schoenf, holo, Laue]

    return pgt_type


def GetRotAxes(pgt_type, new_allrot):
    rot_axes = [[1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 0],
                [0, 1, -1],
                [-1, 0, 1],
                [1, -1, 0],
                [1, 1, 1],      # 10
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1],
                [0, 1, 2],
                [2, 0, 1],
                [1, 2, 0],
                [0, 2, 1],
                [1, 0, 2],
                [2, 1, 0],
                [0, -1, 2],     # 20
                [2, 0, -1],
                [-1, 2, 0],
                [0, -2, 1],
                [1, 0, -2],
                [-2, 1, 0],
                [2, 1, 1],
                [1, 2, 1],
                [1, 1, 2],
                [2, -1, -1],
                [-1, 2, -1],    # 30
                [-1, -1, 2],
                [2, 1, -1],
                [-1, 2, 1],
                [1, -1, 2],
                [2, -1, 1],    
                [1, 2, -1],
                [-1, 1, 2],
                [3, 1, 2],
                [2, 3, 1],
                [1, 2, 3],      # 40
                [3, 2, 1],
                [1, 3, 2],
                [2, 1, 3],
                [3, -1, 2],
                [2, 3, -1],  
                [-1, 2, 3],
                [3, -2, 1],
                [1, 3, -2],
                [-2, 1, 3],
                [3, -1, -2],    # 50
                [-2, 3, -1],
                [-1, -2, 3],
                [3, -2, -1],
                [-1, 3, -2],
                [-2, -1, 3],   
                [3, 1, -2],
                [-2, 3, 1],
                [1, -2, 3],
                [3, 2, -1],
                [-1, 3, 2],     # 60
                [2, -1, 3],
                [1, 1, 3],
                [-1, 1, 3],
                [1, -1, 3],
                [-1, -1, 3],  
                [1, 3, 1],
                [-1, 3, 1],
                [1, 3, -1],
                [-1, 3, -1],
                [3, 1, 1],      # 70
                [3, 1, -1],
                [3, -1, 1],
                [3, -1, -1],
                ]

    axes = np.zeros((3, 3))
    prop_rot = np.zeros((3, 3))
    sum_prop = np.zeros((3, 3))
    inv_rot = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])

    if pgt_type[5] == 'LAUE1':
        # print('here 1')
        axes[0] = rot_axes[0]  
        axes[1] = rot_axes[1]
        axes[2] = rot_axes[2]

    elif pgt_type[5] == 'LAUE2M':
        # print('here 2')
        # get first axis as "b"
        for i in new_allrot:
            prop_rot = np.linalg.det(i) * i

            # search npri = 2 for "b"
            if np.trace(prop_rot) == -1:
                # find e which satisfies W e = e
                for ax in rot_axes:
                    if (np.dot(prop_rot, np.array(ax)) == np.array(ax)).all():
                        axes[1] = ax
                        break
                if not (axes[1] == np.zeros(3)).all():
                    print('axes2 is found')
                    break

        # get second and third axis  S * e2(e3) = 0 as "a" and "c", respectively
        sum_prop = prop_rot + np.dot(prop_rot, prop_rot)
        ortho_candi = []
        for ax in rot_axes:
            if (np.dot(sum_prop, np.array(ax)) == np.zeros(3)).all():
                ortho_candi.append(ax)

        if len(ortho_candi) == 0:
            print('Cannot find another two axes!')
            return None

        # print(ortho_candi)
        axes[0] = ortho_candi[0]
        axes[2] = ortho_candi[1]
        init_det = np.linalg.det(axes)
        tp_axes = axes.copy()

        for i in range(len(ortho_candi)):
            tp_axes[0] = ortho_candi[i]
            for j in range(len(ortho_candi)):
                if ortho_candi[i] != ortho_candi[j]:
                    tp_axes[2] = ortho_candi[j]
                    
                    # print('before change:\n', axes[0], axes[2])
                    if abs(np.linalg.det(tp_axes)) < abs(init_det):
                        init_det = np.linalg.det(tp_axes)
                        axes[0] = tp_axes[0]
                        axes[2] = tp_axes[2]
                        # print('after change:\n', axes[0], axes[2])
                        
        # print('axes', axes)

        flag = True
        for i in range(3):
            if (axes[i] == np.zeros(3)).all():
                flag = False
        print('All axes be found? ', flag)
        if not flag:
            print('Error! Cannot find completed axes!')
            return None

        tp_a = axes.copy()
        
        # when det < 0, only change sec and thr axis
        if np.linalg.det(axes) < 0:  
            tp_a[0] = axes[0]
            axes[0] = axes[2]
            axes[2] = tp_a[0]

        print('axes', axes)

    elif pgt_type[5] == 'LAUEMMM' or pgt_type[5] == 'LAUEM3' or pgt_type[5] == 'LAUEM3M':
        # print('here 3')
        order = None
        if pgt_type[5] == 'LAUEMMM' or pgt_type[5] == 'LAUEM3':
            order = 2
        elif pgt_type[5] == 'LAUEM3M':
            order = 4

        ax_cnt = 0
        for i in new_allrot:
            # get proper rot
            prop_rot = np.linalg.det(i) * i   
            # print('poper_rot ide', prop_rot)
            
            if (np.trace(prop_rot) == -1 and order == 2) or (np.trace(prop_rot) == 1 and order == 4):
                if (prop_rot == np.identity(3)).all():
                    continue
                if ax_cnt < 3:
                    tp_ax = np.zeros(3)
                    for ax in rot_axes:
                        # print('1:', np.dot(prop_rot, np.array(ax)))
                        # print('2: ', np.array(ax))
                        if (np.dot(prop_rot, np.array(ax)) == np.array(ax)).all():
                            tp_ax = ax.copy()
                            # print('yes', tp_ax)
                            break
                    if not ((tp_ax == axes[0]).all() or (tp_ax == axes[1]).all() or (tp_ax == axes[2]).all()):
                        # print('here', ax_cnt)
                        axes[ax_cnt] = tp_ax
                        ax_cnt += 1
                        # print(ax_cnt)

        # check proper rot is right or not (satisfy rotation order)
        tp_prop = prop_rot.copy()   # S = W
        for i in range(order - 1):
            tp_prop = np.dot(tp_prop, prop_rot)  # tp = WW   tp = WWW   tp = WWWWW
        if not (tp_prop == np.identity(3)).all():
            print('Error! Please check n-fold axis again!')

        flag = True
        for i in range(3):
            if (axes[i] == np.zeros(3)).all():
                flag = False
                break
        print('All axes be found? ', flag)
        if not flag:
            print('Error! Cannot find completed axes!')
            return None

        a_index = []
        for i in range(3):
            for cnt, rt in enumerate(rot_axes):
                if (np.array(rt) == axes[i]).all():
                    a_index.append([cnt, rt])
        # print(a_index)

        sort_axe = sorted(a_index, key=lambda setx: setx[0])
        # print('sort axe: ', sort_axe)

        for i in range(3):
            axes[i] = sort_axe[i][1]
        # print(axes)

        tp_a = axes.copy()
        
        #  when det < 0, only change sec and thr axis
        if np.linalg.det(axes) < 0:  
            tp_a[0] = axes[1]
            axes[1] = axes[2]
            axes[2] = tp_a[0]

        print('axes', axes)

    elif pgt_type[5] == 'LAUE4M' or pgt_type[5] == 'LAUE4MMM' or pgt_type[5] == 'LAUE3' or pgt_type[5] == 'LAUE3M' or pgt_type[5] == 'LAUE6M' or pgt_type[5] == 'LAUE6MMM':
        # print('here 4')
        order = None
        if pgt_type[5] == 'LAUE4M' or pgt_type[5] == 'LAUE4MMM':
            order = 4
        elif pgt_type[5] == 'LAUE3' or pgt_type[5] == 'LAUE3M' or pgt_type[5] == 'LAUE6M' or pgt_type[5] == 'LAUE6MMM':
            order = 3

        # find first axis as "c"
        for i in new_allrot:
            # get proper rot
            prop_rot = np.linalg.det(i) * i     
            
            if (np.trace(prop_rot) == 1 and order == 4) or (np.trace(prop_rot) == 0 and order == 3):
                for ax in rot_axes:
                    if (np.dot(prop_rot, np.array(ax)) == np.array(ax)).all():
                        axes[2] = ax
                        break
                if not (axes[2] == np.zeros(3)).all():
                    print('axes3 is found')
                    break

        # find second axis S * e2 = 0 as "a"

        tp_prop = prop_rot.copy()   # S = W1
        sum_prop = prop_rot.copy()
        for i in range(order - 1):
            # tp = WW   tp = WWW   tp = WWWWW
            # S = W + W*W  S = W + WW + WWW  S = W + WW + WWW + WWWW
            tp_prop = np.dot(tp_prop, prop_rot)  
            sum_prop += tp_prop     
       
        if not (tp_prop == np.identity(3)).all():
            print('Error! Please check n-fold axis again!')

        ortho_candi = []
        for ax in rot_axes:
            if (np.dot(sum_prop, np.array(ax)) == np.zeros(3)).all():
                ortho_candi.append(ax)

        if len(ortho_candi) == 0:
            print('Cannot find another two axes!')
            return None

        cnnn = 0
        for i in range(len(ortho_candi)):
            # set axes1 in ortho_candidates as second axis
            axes[0] = ortho_candi[i]
            
            # find third axis  e3 = W * e2
            axes[1] = np.dot(prop_rot, axes[0])
            
            is_found = 0
            for j in range(len(ortho_candi)):
                if (axes[1] == np.array(ortho_candi[j])).all() or (axes[1] == -np.array(ortho_candi[j])).all():
                    # print('axes1 is found')
                    # print('axes2 is found')
                    is_found = 1
                    break
            if not is_found:
                continue
            if abs(np.linalg.det(axes)) < 4:
                cnnn += 1
                # print('cnnn', cnnn)
                break

        tp_a = axes.copy()
        
        # when det < 0, only change sec and thr axis
        if np.linalg.det(axes) < 0:  
            tp_a[0] = axes[0]
            axes[0] = axes[1]
            axes[1] = tp_a[0]
            print('axes1 is found')
            print('axes2 is found')

        flag = True
        for i in range(3):
            if (axes[i] == np.zeros(3)).all():
                flag = False
        print('All axes be found? ', flag)
        if not flag:
            print('Error! Cannot find completed axes!')
            return None

        print('axes', axes)

    else:
        print('here error')
        return None

    return axes


def GetPointGroup(prim_latt, num, pos, dictp):
    """
        Get point group and its symmetry operations
    """
    # get all space group operations: rotation + translation
    # o_latt, o_pos, o_num, dictp = delaunay.StructRead()
    # prim_latt, num, pos = primitive.GetPrimLattice(o_latt, o_pos, o_num, dictp)
    allrot = SymRotation(prim_latt)
   
    # print('check2: ', len(allrot))
    symmetry = []
    new_allrot = []
    for tprot in allrot:
        trans = primitive.SymTrans(tprot, dictp, num, pos, prim_latt)
        # print('check trans:\n', len(trans))
        # print(trans == [])
        
        if trans != []:
            new_allrot.append(tprot)
            for i in range(len(trans)):
                # alltrans.append(trans[i])
                # new_allrot.append(tprot)  
                # save the rotation which the corresponding trans is not None!
                symmetry.append([tprot, trans[i]])

    # print('CH1: ', len(symmetry))

    pgt_type = GetRotTable(new_allrot)
    # print(pgt_type[5])
    axis = GetRotAxes(pgt_type, new_allrot)
    axis = np.array(axis, dtype=int)
    
    # print('axis', axis)
    return symmetry, pgt_type, axis