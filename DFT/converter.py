#!/usr/bin/env python

import sys
import os


old = sys.argv[1]

infile = open(old,'r')
inlines = infile.readlines()
infile.close()

#print inlines[2][0:6]
natom = int( inlines[2][0:6] )
print 'number of atoms = ', natom

stepsX = int( inlines[3][0:6] )
stepsY = int( inlines[4][0:6] )
stepsZ = int( inlines[5][0:6] )
print 'dimensions of the charge density : ', stepsX, ' x ', stepsY, ' x ', stepsZ

#outlines = inlines[0:natom+6]

ex = [0]*4
ey = [0]*4
ez = [0]*4

i = 3
while i < 6 :
	
        j = 0
        c = inlines[i][j]

        while c == ' ' :
                j = j + 1
                c = inlines[i][j]

	tmp = ''
        while c != ' ' :
                tmp = tmp + c
                j = j + 1
                c = inlines[i][j]
	steps = int( tmp )

        while c == ' ' :
                j = j + 1
                c = inlines[i][j]

	tmp = ''
        while c != ' ' :
                tmp = tmp + c
                j = j + 1
                c = inlines[i][j]
	ex[i-2] = float( tmp )

        while c == ' ' :
                j = j + 1
                c = inlines[i][j]

        tmp = ''
        while c != ' ' :
                tmp = tmp + c
                j = j + 1
                c = inlines[i][j]
        ey[i-2] = float( tmp )

        while c == ' ' :
                j = j + 1
                c = inlines[i][j]

        tmp = ''
        while c != '\n' :
                tmp = tmp + c
                j = j + 1
                c = inlines[i][j]
        ez[i-2] = float( tmp )

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

	i = 0
	c = ''
	while c != '\n' :

       		c = inlines[j+offset][i]

        	while c == ' ' :
       	        	i = i + 1
                	c = inlines[j+offset][i]

	        tmp = ''
       		while ( c != ' ' and c != '\n' ) :
               		tmp = tmp + c
	                i = i + 1
       		        c = inlines[j+offset][i]

		entry[k] = float( tmp )
#		print '   k = ', k, '  density = ', entry[k]
		k = k + 1

	j = j + 1

#	n = n + 1


#sum1 = 0
#sum2 = 0
#sum3 = 0
sum0 = 0
i = 0
while i < size :
#	sum1 = sum1 + entry1[i]
#	sum2 = sum2 + entry2[i]
#	sum3 = sum3 + entry3[i]
	sum0 = sum0 + entry[i]
	i = i + 1

print sum0
#print sum1, sum2

# just added 4/13/12, double check this
#
sum0 = sum0 / 4

#density = [[[0]*stepsZ]*stepsY]*stepsX
density = [0]*size

i  = 0
ix = 0
while ix < stepsX :
#	print " ix = ", ix
	iy = 0
	while iy < stepsY :
		iz = 0
		while iz < stepsZ :
#			print '   i = ', i, '  density = ', entry[i]
#			density[ix][iy][iz] = entry[i]
			j = iz*(stepsX*stepsY) + iy*(stepsX) + ix
			density[j] = entry[i]
#			density[j] = entry[i] / 6.74833304
#                           division is for anstrom to bohr conversion (cubic)
#                        print '   j = ', j, '  density = ', density[j]
			i = i + 1
			iz = iz + 1
		iy = iy + 1
	ix = ix + 1



outlines = [''] * (size+1)
outlines[0] = '    i1    i2    i3      data \n'


i  = 0
iz = 0
while iz < stepsZ :
#	print " iz = ", iz
        iy = 0
        while iy < stepsY :
                ix = 0
                while ix < stepsX :
#                        print '   i = ', i, '  density = ', density[ix][iy][iz]
#                        outlines = outlines + ['']
			j = iz*(stepsX*stepsY) + iy*(stepsX) + ix
			outlines[i+1] = '     ' + str(ix+1) + '     ' + str(iy+1) + '     ' + str(iz+1) + '     ' + str(density[j]) + '\n'
                        i = i + 1
                        ix = ix + 1
                iy = iy + 1
        iz = iz + 1


outfile = open('rhoofr','w')
outfile.writelines(outlines)
outfile.close()


