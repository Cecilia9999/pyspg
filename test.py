#!/usr/bin/env python
# -*- coding=utf-8 -*-


import numpy as np
import delaunay
import niggli
import csv
'''
c2 = np.array([[0.0000000000, 0.0000000000, 1.0000000000],
 [1.0000000000, 0.0000000000, 0.0000000000],
 [0.0000000000, 1.0000000000, 0.0000000000]])

c1 = np.array([[0, 0, 1],
       [1, 0, 0],
       [0, 1, 0]])
if (c1 == c2).all():
    print('true')

a = np .array([[1, 0, 1],
     [1, 2, 1],
    [-1, 0, 1]])


b = np.array([[0.0000000000, 2.5185500000, 2.5185500000],
 [2.5185500000, 0.0000000000, 2.5185500000],
 [2.5185500000, 2.5185500000, 0.0000000000]])

b1 = np.array([[2.5185500000, -2.5185500000, 0.0000000000],
 [-2.5185500000, -0.0000000000, -2.5185500000],
 [2.5185500000, 2.5185500000, 0.0000000000]])

t1 = np.transpose(np.dot(b1, np.linalg.inv(b)))   # P vertical
q1 = np.linalg.inv(t1)

m = np.array([[1.0000000000, 1.0000000000, -1.0000000000],
 [0.0000000000, 2.0000000000, 0.0000000000],
 [1.0000000000, 1.0000000000, 1.0000000000]])

print(np.cos(180))
'''

l1 = [[7.17851431, 0, 0],  # a
           [0, 3.99943947, 0],  # b
           [0, 0, 8.57154746]]


l2 = [[5.0759761474456697, 5.0759761474456697, 0],  # a
           [-2.8280307701821314, 2.8280307701821314, 0],  # b
           [0, 0, 8.57154746]]

l3 = [[8.57154746, 0, 0],  # a                  # // swap a - c
           [0, -3.99943947, 0],  # b
           [0, 0, 7.17851431]]  # c


pos = [[0.0, 0.84688439, 0.1203133],
          [0.0, 0.65311561, 0.6203133],
          [0.0, 0.34688439, 0.3796867],
          [0.0, 0.15311561, 0.8796867],
          [0.5, 0.34688439, 0.1203133],
          [0.5, 0.15311561, 0.6203133],
          [0.5, 0.84688439, 0.3796867],
          [0.5, 0.65311561, 0.8796867]]

pos2 = [[0.1203133, 0.84688439, 0.0],            # // swap a - c
          [0.6203133, 0.65311561, 0.0],
          [0.3796867, 0.34688439, 0.0],
          [0.8796867, 0.15311561, 0.0],
          [0.1203133, 0.34688439, 0.5],
          [0.6203133, 0.15311561, 0.5],
          [0.3796867, 0.84688439, 0.5],
          [0.8796867, 0.65311561, 0.5]]


p = np.transpose(np.dot(l3, np.linalg.inv(l1)))
print(p)

new_pos = []
for i in pos:
    new_pos.append(np.dot(i, np.linalg.inv(p)))

for i in new_pos:
    print(list(i))


'''
p = np.transpose(np.dot(l2, np.linalg.inv(l1)))
print(p)

new_pos = []
for i in pos:
    new_pos.append(np.dot(np.linalg.inv(p), i))
#print(new_pos)
for i in new_pos:
    print(list(i))

Pc = [[0.5, 0.5, 0],
      [-0.5, 0.5, 0],
      [0, 0, 1]]


'''




b = np.array([[6.6518325806, 0.0000000000, 0.0000000000],
              [3.3259162903, 5.7606559965, 0.0000000000],
              [3.3259162903, 1.9202186655, 5.4311985589]])    # input prim

b1 = np.array([[0.0000000000, 4.7035559250, 4.7035559250],
               [4.7035559250, 0.0000000000, 4.7035559250],
               [4.7035559250, 4.7035559250, 0.0000000000]])    # standard prim

b2 = np.array([[-3.3259162903, 5.7606559965, 0.0000000000],
               [-3.3259162903, -5.7606559965, -0.0000000000],   # delaunay
               [3.3259162903, 1.9202186655, 5.4311985589]])

b3 = np.array([[-0.0000000000, 7.6808746620, 5.4311985589],    # bravais
                [-6.6518325806, -3.8404373310, 5.4311985589],
                [6.6518325806, -3.8404373310, 5.4311985589]])

t_mat = np.array([[0.5, 0.5, 0.0],
                   [0, 0.5, 0.5],
                   [0.5, 0, 0.5]])

print('b--b2:\n', np.transpose(np.dot(b2, np.linalg.inv(b))))   # input --> delaunay
print('b2--b3:\n', np.transpose(np.dot(b3, np.linalg.inv(b2))))   # delaunay --> bravais
print('b--b3:\n', np.transpose(np.dot(b3, np.linalg.inv(b))))
c1 = np.linalg.inv(np.transpose(np.dot(b2, np.linalg.inv(b))))
mm1 = np.dot(np.transpose(c1), b3)
print('check ', mm1)
mm2 = np.array([[ 1, 1, -1],
                [ 0,  2,  0],
                [ 1,  1,  1]])
print(np.dot(mm1, mm2))
print('a:\n', np.dot(np.transpose(t_mat), b))

reduc_b = niggli.Niggli(b3)
print('de:', reduc_b)


print(np.linalg.det(b)/np.linalg.det(b1))
F_mat = np.array ([[0, 1. / 2, 1. / 2],
                   [1. / 2, 0, 1. / 2],
                   [1. / 2, 1. / 2, 0]])



c = np.transpose(np.dot(b3, np.linalg.inv(b)))   # P vertical
d = np.dot(np.linalg.inv(c), F_mat)

print(c, d)


t0 = np.array([[1.0000000000, 1.0000000000, -1.0000000000],
               [0.0000000000, 2.0000000000, 0.0000000000],
               [1.0000000000, 1.0000000000, 1.0000000000]])
print(np.linalg.inv(t0))

t1 = np.array([[6.6518325806, 0.0000000000, 0.0000000000],
 [3.3259162903, 5.7606559965, 0.0000000000],
 [3.3259162903, 1.9202186655, 5.4311985589]])

t2 = np.array([[0.5, 0.5, 0.0],
               [0.0, 0.5, 0.5],
               [0.5, 0.0, 0.5]])

print(np.dot(np.transpose(np.linalg.inv(t2)), t1))
'''
for i in range (3):
    for j in range (3):
        if c[i, j] < 0.0:
            c[i, j] = int (c[i, j] - 0.5)
        else:
            c[i, j] = int (c[i, j] + 0.5)
d = np.linalg.inv(c)
print(c)
print(d)



with open("HallDeltaTrans", 'r') as H:
    tmp = H.readlines()
    #print(tmp[1])
    new = []
    for t in tmp[1:]:
        #print(t)
        e_t = t.split()
        if e_t == ['};']:
            break
        #print(e_t)
        if e_t[3].startswith('('):
            if e_t[3] == '(':
                e_t2 = e_t[4:]
            else:
                e_t2 = e_t[3:]
        elif e_t[4].startswith('('):
            if e_t[4] == '(':
                e_t2 = e_t[5:]
            else:
                e_t2 = e_t[4:]
        #print(e_t2)
        num = e_t2[0].strip(")")  # // hall number
        if num.startswith('('):
            num = num.strip('(')
        #print(num)
        sym = []
        #sym.append(e_t2[2])
        for i in e_t2[2:]:
            t = []
            if i == '[':
                continue
            if i == '*/':
                break

            if len(i) >= 4:
                t = i.split(',')
                for m in t:
                    if m != '':
                        if m[-1] == ']':
                            sym.append(int(m.strip(']')))
                        else:
                            sym.append(int(m))
            else:
                if i[-1] == ']':
                    sym.append(int(i.strip(']')))
                else:
                    sym.append(int(i.strip(',')))

        #print(sym)
        new.append([int(num), sym])

    a = []
    for i in new:
        if len(i[1]) != 12:
            a.append(1)
    print(len(a))

    with open('newhall', 'w') as W:
        for i in new:
            W.write(str(i) + ',' + '\n')
        W.close()
    H.close()

    
    cn = 1 
    for i in range(len(new)):
        tp_vec = []
        if i[0] == cn:
            tp_vec.append(cn, i[1])
        elif i[0] != cn:
            cn += 1
    print(new)
    '''






'''
rot = [[1, 0, 0],
       [0, 1, 0],
       [0, 0, 1],
       [-1, 0, 0],
       [0, -1, 0],  # /* 5 */
       [0, 0, -1],
       [0, 1, 1],
       [1, 0, 1],
       [1, 1, 0],
       [0, -1, -1],  # /* 10 */
       [-1, 0, -1],
       [-1, -1, 0],
       [0, 1, -1],
       [-1, 0, 1],
       [1, -1, 0],  # /* 15 */
       [0, -1, 1],
       [1, 0, -1],
       [-1, 1, 0],
       [1, 1, 1],
       [-1, -1, -1],  # /* 20 */
       [-1, 1, 1],
       [1, -1, 1],
       [1, 1, -1],
       [1, -1, -1],
       [-1, 1, -1],  # /* 25 */
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
                rotations.append([i, j, k, tp_rot])
print(rotations)

with open("rotfile", 'w') as ROT:
    ROT.write ("ALL OPERATIONS:\n")
    for i in rotations:
        ROT.write(str(i[3]) + '\n')
    ROT.close()




a = [1, 2, ]
for i in range(2, 17):
    a.append(3)
for i in range(17, 56):
    a.append(4)
for i in range(56, 107):
    a.append(5)
for i in range(107, 124):
    a.append(6)
for i in range(124, 227):
    a.append(7)
for i in range(227, 349):
    a.append(8)
for i in range(349, 355):
    a.append(9)
for i in range(355, 357):
    a.append(10)
for i in range(357, 366):
    a.append(11)
for i in range(366, 376):
    a.append(12)
for i in range(376, 388):
    a.append(13)
for i in range(388, 400):
    a.append(14)
for i in range(400, 430):
    a.append(15)
for i in range(430, 435):
    a.append(16)
for i in range(435, 438):
    a.append(17)
for i in range(438, 446):
    a.append(18)
for i in range(446, 454):
    a.append(19)
for i in range(454, 462):
    a.append(20)
for i in range(462, 468):
    a.append(21)
a.append(22)
for i in range(2):
    a.append(23)
for i in range(471, 477):
    a.append(24)
for i in range(477, 481):
    a.append(25)
for i in range(481, 485):
    a.append(26)
for i in range(485, 489):
    a.append(27)
for i in range(489, 494):
    a.append(28)
for i in range(494, 503):
    a.append(29)
for i in range(503, 511):
    a.append(30)
for i in range(511, 517):
    a.append(31)
for i in range(517, 531):
    a.append(32)

print(len(a))
with open('spg2.csv', 'rt') as csvfile:
    rows = csv.reader(csvfile)
    with open('newspg.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for i, row in enumerate(rows):
           row.append(a[i])
           writer.writerow(row)
    #for i in range(1, 531):

    f.close()
    csvfile.close()


a = np.array([[0.3, 0, 1], [1, 1.5, 2], [2, 2, 1]])
b = np.array([[1, 2, 1], [0, 1, 2], [1, 0, 0]])
c = np.array([1, 2, 1])
d = np.array([1.5, 2, 0.5])
m = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        m[i, j] = -a[i, j]

def SS(a):
  a = a + 1
  return a

a = 4
a = SS(a)
print(a)


print(m)
print(a)

a = np.zeros((3, 3))
b = np.array([[1, 1, 1], [-1, 0, 1], [0, 2, 1]])

a[0] = [1, 1, 1]
print(-a[0])
b[0] = [0, 0, 3]
print(a[0])




t = []
for i in range(-3, 4):
    for j in range(-3, 4):
        for k in range(-3, 4):
            t.append([i, j, k])
for i in t:
    if i == [0, 0, 0]:
        del i
print(len(t))
ct = t.copy()
m = []

for i in ct:
    for j in t:
        if (np.array(i) == -np.array(j)).all():
            print(i, j)
            ct.remove(j)
            t.remove(j)

print(len(ct))
print(ct)

for c in range(171):
    for i in ct:
        count2 = 0
        count3 = 0
        for j in i:
            if j == 2 or j == -2:
                count2 += 1
            elif j == 3 or j == -3:
                count3 += 1
        if count2 > 1 or count3 > 1:
            ct.remove(i)

for c in range(103):
    for i in ct:
        count0 = 0
        for j in i:
            if j == 3 or j == -3:
                count0 += 1
            if count0 > 0:
                if 0 in i:
                    ct.remove(i)
                    break

for c in range(103):
    for i in ct:
        count0 = 0
        count1 = 0
        for j in i:
            if j == 0:
                count1 += 1
            if count1 == 2:
                if 2 in i or -2 in i:
                    ct.remove(i)
                    break

print(len(ct))
print(ct)

print(np.dot(a, c))
print(np.dot(a, np.transpose(c)))

for i in range(3):
    for j in range(3):
        if a[i][j] < 0.0:
            a[i][j] = int(a[i][j] - 0.5)
        else:
            a[i][j] = int(a[i][j] + 0.5)
print(a)

for i in range(3):
    for j in range(3):
        b[i][j] = np.round(b[i][j])

print(b)

print(int(0.5), int(0.7), int(-0.4), int(-0.5))
'''
