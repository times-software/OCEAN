#!/usr/bin/env python
# Copyright (C) 2014 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#


import sys
import os


old = sys.argv[1]

infile = open(old,'r')
inlines = infile.readlines()
infile.close()

natom = int( inlines[2][0:6] )
print 'number of atoms = ', natom

stepsX = int( inlines[3][0:6] )
stepsY = int( inlines[4][0:6] )
stepsZ = int( inlines[5][0:6] )
print 'dimensions of the charge density : ', stepsX, ' x ', stepsY, ' x ', stepsZ

ex = [0]*4
ey = [0]*4
ez = [0]*4

i = 3
while i < 6 :

	tmp = inlines[i].split()
	steps = tmp[0]
	ex[i-2] = tmp[1]
	ey[i-2] = tmp[2]
	ez[i-2] = tmp[3]
	
	i = i + 1

#print steps
print 'the lattice vectors are:'
print '  e_x = ', ex[1], ex[2], ex[3]
print '  e_y = ', ey[1], ey[2], ey[3]
print '  e_z = ', ez[1], ez[2], ez[3]

size = stepsX * stepsY * stepsZ


n = 1
#while n < 4 :

#	inlines = inlines1
offset = natom + 6

entry = [0]*size

k = 0
j = 0
while j < ( len(inlines) - offset) :

	tmp = inlines[j+offset].split()

	i = 0
	while i < len(tmp) :
		entry[k] = float( tmp[i] )
		k = k + 1
		i = i + 1

	j = j + 1


sum0 = 0
i = 0
while i < size :
	sum0 = sum0 + entry[i]
	i = i + 1

#print sum0
#print sum1, sum2

# just added 4/13/12, double check this
#
sum0 = sum0 / 4

#density = [[[0]*stepsZ]*stepsY]*stepsX
density = [0]*size

i  = 0
ix = 0
while ix < stepsX :
	iy = 0
	while iy < stepsY :
		iz = 0
		while iz < stepsZ :
			j = iz*(stepsX*stepsY) + iy*(stepsX) + ix
			density[j] = entry[i]
			i = i + 1
			iz = iz + 1
		iy = iy + 1
	ix = ix + 1


outlines = [''] * (size+1)
outlines[0] = '    i1    i2    i3      data \n'


i  = 0
iz = 0
while iz < stepsZ :
        iy = 0
        while iy < stepsY :
                ix = 0
                while ix < stepsX :
			j = iz*(stepsX*stepsY) + iy*(stepsX) + ix
			outlines[i+1] = '     ' + str(ix+1) + '     ' + str(iy+1) + '     ' + str(iz+1) + '     ' + str(density[j]) + '\n'
                        i = i + 1
                        ix = ix + 1
                iy = iy + 1
        iz = iz + 1


outfile = open('rhoofr','w')
outfile.writelines(outlines)
outfile.close()


