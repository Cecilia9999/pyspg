# pyspg

A python version of spglib (Astushi Togo), which is used for finding symmetry of crystal.  



## Functions

Main functions of pyspg are: 

- Change input cell to standard (conventional/primitive) cell

- Find point group/ space goup of input cell

- Other symmetry information of input cell



## Modulus

- delaunay.py 	:  Delauney transformation (find smallest cell)

- niggli.py		:  Niggli transformation (find smallest cell)

- primitive.py	:  Find primitive cell

- pointgp.py	:  find point group of input cell

- halldata.py 	:  Data of Hall symmetry operations

- hallsymb.py 	:  Find Hall symbol

- matchHall.py	:  Match Hall symbol

- spgroup.py	:  Find space group of input cell

- Inputcell.py	:  Handle input "POSCAR" format file and change it to standard/primitive cell

- change of b	:  Change axis for some certain type (refer: spglib package, and update it)

Some database and documents are in ./supports/



## Reference

spglib: https://github.com/spglib/spglib



### Author

Cecilia9999:  https://github.com/Cecilia9999/

