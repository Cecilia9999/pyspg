#!/usr/bin/env python
# -*- coding=utf-8 -*-
'''
This code is for generating standardized & primitive cell
input: POSCAR
output: POSCAR-stand, POSCAR-prim

NOTE! This script is based on Pymatgen(different from seekpath)

--- syx
'''

import os
import pymatgen
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

stru = Structure.from_file("POSCART2")
spg = SpacegroupAnalyzer(stru, symprec=0.01, angle_tolerance=5)
#fp = spg.find_primitive()
fp = spg.get_primitive_standard_structure()
stan = spg.get_conventional_standard_structure()

with open("POSCAR-stand", 'w') as ST:
    ST.write("POSCAR-STAN\n" + "1.0\n")
    ST.write(str(stan.lattice) + '\n')

    # count {species: nums}
    spe = {}
    for i in stan.species:
        spe[i] = spe.get(i, 0) + 1
    #print(spe)
    for i in spe.keys():
        ST.write(str(i) + " ")
    ST.write('\n')

    for i in spe.values():
        ST.write(str(i) + " ")
    ST.write('\n')

    ST.write("Direct\n")

    [rows, cols] = stan.frac_coords.shape
    for i in range(rows):
        for j in range(cols):
            ST.write(str('{:.12f}'.format(stan.frac_coords[i, j])) + "  ")
        ST.write('\n')
    ST.close()

with open("POSCAR-prim", 'w') as PR:
    PR.write("POSCAR-PRIM\n" + "1.0\n")
    PR.write(str(fp.lattice) + '\n')

    # count {species: nums}
    spe = {}
    for i in fp.species:
        spe[i] = spe.get(i, 0) + 1

    for i in spe.keys():
        PR.write(str(i) + " ")
    PR.write('\n')

    for i in spe.values():
        PR.write(str(i) + " ")
    PR.write('\n')

    PR.write("Direct\n")

    [rows, cols] = fp.frac_coords.shape
    for i in range(rows):
        for j in range(cols):
            PR.write(str('{:.12f}'.format(fp.frac_coords[i, j])) + "  ")
        PR.write('\n')
    PR.close()
