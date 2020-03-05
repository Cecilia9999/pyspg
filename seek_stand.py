#!/usr/bin/env python
# -*- coding=utf-8 -*-
'''
This code is for generating standardized & primitive cell
input: POSCAR
output: POSCAR-stand, POSCAR-prim

NOTE! This script is based on seekpath

--- syx
'''

import os.path
import seekpath
import spglib
import numpy as np
import pymatgen.symmetry.analyzer as psa

def get_atom(symbol):
    symbol_map = {
            "H": 1,
            "He": 2,
            "Li": 3,
            "Be": 4,
            "B": 5,
            "C": 6,
            "N": 7,
            "O": 8,
            "F": 9,
            "Ne": 10,
            "Na": 11,
            "Mg": 12,
            "Al": 13,
            "Si": 14,
            "P": 15,
            "S": 16,
            "Cl": 17,
            "Ar": 18,
            "K": 19,
            "Ca": 20,
            "Sc": 21,
            "Ti": 22,
            "V": 23,
            "Cr": 24,
            "Mn": 25,
            "Fe": 26,
            "Co": 27,
            "Ni": 28,
            "Cu": 29,
            "Zn": 30,
            "Ga": 31,
            "Ge": 32,
            "As": 33,
            "Se": 34,
            "Br": 35,
            "Kr": 36,
            "Rb": 37,
            "Sr": 38,
            "Y": 39,
            "Zr": 40,
            "Nb": 41,
            "Mo": 42,
            "Tc": 43,
            "Ru": 44,
            "Rh": 45,
            "Pd": 46,
            "Ag": 47,
            "Cd": 48,
            "In": 49,
            "Sn": 50,
            "Sb": 51,
            "Te": 52,
            "I": 53,
            "Xe": 54,
            "Cs": 55,
            "Ba": 56,
            "La": 57,
            "Ce": 58,
            "Pr": 59,
            "Nd": 60,
            "Pm": 61,
            "Sm": 62,
            "Eu": 63,
            "Gd": 64,
            "Tb": 65,
            "Dy": 66,
            "Ho": 67,
            "Er": 68,
            "Tm": 69,
            "Yb": 70,
            "Lu": 71,
            "Hf": 72,
            "Ta": 73,
            "W": 74,
            "Re": 75,
            "Os": 76,
            "Ir": 77,
            "Pt": 78,
            "Au": 79,
            "Hg": 80,
            "Tl": 81,
            "Pb": 82,
            "Bi": 83,
            "Po": 84,
            "At": 85,
            "Rn": 86,
            "Fr": 87,
            "Ra": 88,
            "Ac": 89,
            "Th": 90,
            "Pa": 91,
            "U": 92,
            "Np": 93,
            "Pu": 94,
            "Am": 95,
            "Cm": 96,
            "Bk": 97,
            "Cf": 98,
            "Es": 99,
            "Fm": 100,
            "Md": 101,
            "No": 102,
            "Lr": 103,
            "Rf": 104,
            "Db": 105,
            "Sg": 106,
            "Bh": 107,
            "Hs": 108,
            "Mt": 109,
            "Ds": 110,
            "Rg": 111,
            "Cn": 112,
            "Uut": 113,
            "Uuq": 114,
            "Uup": 115,
            "Uuh": 116,
            "Uus": 117,
            "Uuo": 118,
        }
    return symbol_map[symbol]

filename = input("Please input the POSCAR name: ")
filename = './testfile/' + filename
if os.path.isfile(filename):
    print("Using POSCAR to calculate...")
else:
    print("Error!\nPlease prepare a POSCAR!")

with open(filename, 'r') as POS:
    tmp = POS.readlines()
    tmp1 = tmp[1].split()
    latt_const = float(tmp1[0])
    #print(latt_const)
    tmp3 = []
    for i in tmp[2:5]:
        tmp2 = []
        for j in i.split():
            tmp2.append(float(j) * latt_const)
        tmp3.append(tmp2)
    # print(tmp3)
    lattice = np.array(tmp3)
    np.set_printoptions(formatter={'float': '{:0.10f}'.format})
    print('lattice:')
    print(lattice)
    species = tmp[5]
    tmp7 = species.split()  #[Si, O]
    tmp4 = [int(i) for i in tmp[6].split()]  # [4, 8]
    dictp = dict(zip(tmp7, tmp4))
    print(dictp)
    numbers = []
    for i in dictp.keys():
        for k in range(dictp[i]):
            numbers.append(get_atom(i))  #【40, 40, 40, 40, 8, 8, 8, 8, 8, 8, 8, 8】
    print('numbers:')
    print(numbers)

    tmp5 = tmp[7].split()
    # print(tmp5)
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

    cell = (lattice, positions, numbers)
    print('space group:\n' + spglib.get_spacegroup(cell))
    path = seekpath.get_path(cell)
    #print(path)
    print(spglib.get_symmetry_dataset(cell))
    print(spglib.delaunay_reduce(lattice))
    #print(spglib.find_primitive(cell))


output = filename + "-stand"
with open(output, 'w') as ST:
    ST.write("POSCAR-STAN\n" + "1.0\n")
    for i in path['conv_lattice']:
        for j in i:
            ST.write(str('{:.10f}'.format(j)) + ' ')
        ST.write('\n')

    ST.write(species)

    spe1 = {}
    cnt1 = {}
    conv = path['conv_positions'].tolist()
    m = 0
    for i in path['conv_types']:
        if i not in cnt1.keys():
            cnt1[i] = 1
            spe1[i] = [conv[m], ]
            m = m + 1
        else:
            cnt1[i] = cnt1[i] + 1
            spe1[i].append(conv[m])
            m = m + 1
    #print(spe1)

    for i in cnt1.values():
        ST.write(" " + str(i) + " ")
    ST.write('\n')

    ST.write("Direct\n")

    for i in spe1.values():
        for j in i:
            for k in j:
                ST.write(str('{:.12f}'.format(k)) + ' ')
            ST.write('\n')

    ST.close()

    '''
    for i in path['conv_positions']:
        for j in i:
            ST.write(str('{:.12f}'.format(j)) + ' ')
        ST.write('\n')
    ST.close()
    '''
output2 = filename + "-prim"
with open(output2, 'w') as PR:
    PR.write("POSCAR-PRIM\n" + "1.0\n")
    for i in path['primitive_lattice']:
        for j in i:
            PR.write(str('{:.10f}'.format(j)) + ' ')
        PR.write('\n')

    PR.write(species)

    spe2 = {}
    cnt2 = {}
    pri = path['primitive_positions'].tolist()
    m = 0
    for i in path['primitive_types']:
        if i not in cnt2.keys():
            cnt2[i] = 1
            spe2[i] = [pri[m], ]
            m = m + 1
        else:
            cnt2[i] = cnt2[i] + 1
            spe2[i].append(pri[m])
            m = m + 1
    #print(spe2)

    for i in cnt2.values():
        PR.write(" " + str(i) + " ")
    PR.write('\n')

    PR.write("Direct\n")

    for i in spe2.values():
        for j in i:
            for k in j:
                if k > 1.0:
                    k = k - 1.0
                elif k < 0:
                    k = k + 1
                PR.write(str('{:.12f}'.format(k)) + ' ')
            PR.write('\n')

    PR.close()