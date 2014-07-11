#!/usr/bin/env python

import sys
import os


# read ibrav
infile = open('ibrav','r')
input = infile.readlines()
infile.close()

ibrav = input[0].rstrip()

# read rscale
infile = open('rscale','r')
input = infile.readlines()
infile.close()

acell = input[0].split()

# read rprim
infile = open('rprim','r')
rprim = infile.readlines()
infile.close()

outfile = open('cellinfo','w')
outlines = ['']*6

# case statement on ibrav

if   ibrav=='0':
#	print 'ibrav = 0'
	print 'nothing to be done for ibrav = 0'

elif ibrav=='1':
#	print 'ibrav = 1'
	outlines[0] = '   celldm(1) = ' + acell[0] + '\n'

elif ibrav=='2':
#	print 'ibrav = 2'
	outlines[0] = '   celldm(1) = ' + acell[0] + '\n'

elif ibrav=='3':
#	print 'ibrav = 3'
	outlines[0] = '   celldm(1) = ' + acell[0] + '\n'

elif ibrav=='4':
#	print 'ibrav = 4'
	outlines[0] = '   celldm(1) = ' + acell[0] + '\n'
        outlines[1] = '   celldm(3) = ' + str( float( acell[2] ) / float( acell[0] ) ) + '\n'

elif ibrav=='5':
#	print 'ibrav = 5'
	print 'ibrav = 5 not yet implemented'

elif ibrav=='6':
#	print 'ibrav = 6'
        outlines[0] = '   celldm(1) = ' + acell[0] + '\n'
       	outlines[1] = '   celldm(3) = ' + str( float( acell[2] ) / float( acell[0] ) ) + '\n'

elif ibrav=='7':
#	print 'ibrav = 7'
	print 'ibrav = 7 not yet implemented'

elif ibrav=='8':
#	print 'ibrav = 8'
	outlines[0] = '   celldm(1) = ' + acell[0] + '\n'
	outlines[1] = '   celldm(2) = ' + str( float( acell[1] ) / float( acell[0] ) ) + '\n'
	outlines[2] = '   celldm(3) = ' + str( float( acell[2] ) / float( acell[0] ) ) + '\n'

else:
	print 'ibrav > 8 not yet implemented'

outfile.writelines(outlines)
outfile.close()
