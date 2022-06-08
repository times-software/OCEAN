# Jvinson 2022
"""
Using ASE and SPGLIB, read a variety of input formats (whatever ASE read 
supports) and attempt to build an OCEAN input file.

Various defaults can be stored in the ocean.json file, a copy of which is 
co-located with this script. If ocean.json is present in the directory in 
which this script is run, that file will take precedent. 

No error checking is currently done.
"""
from ase.atoms import Atoms
from ase.io import read
from ase.units import Bohr
from spglib import standardize_cell
import os
import sys
import json


def write_ocean_in( filename: str, atoms: Atoms, input_data: dict ):

    filename = os.path.expanduser(filename)
    mode = 'w'
    fd = open(filename, mode)

    input_data_str = []
    for key in input_data:
      input_data_str.append( str(key) + ' { ' + str(input_data[key]) + ' }\n' )

    fd.write( ''.join(input_data_str ))

    if any(atoms.get_initial_magnetic_moments()):
        raise NameError( 'Spin=2 not implemented yet' )

    species = sorted(set(atoms.numbers))
    fd.write('znucl {{ {} }}\n'.format(' '.join(str(Z) for Z in species)))
    fd.write('typat')
    fd.write('{\n')
    types = []
    for Z in atoms.numbers:
        for n, Zs in enumerate(species):
            if Z == Zs:
                types.append(n + 1)
    n_entries_int = 20  # integer entries per line
    for n, type in enumerate(types):
        fd.write(' %d' % (type))
        if n > 1 and ((n % n_entries_int) == 1):
            fd.write('\n')
    fd.write(' }\n')

    atomic_positions_str = []
    for atom in atoms:
        atomic_positions_str.append( '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}\n'.format(
              coords=[atom.a, atom.b, atom.c] ))

    fd.write( 'xred {\n' )
    fd.write( ''.join(atomic_positions_str))
    fd.write( '}\n' )

    fd.write( 'acell {{ {acell[0]} {acell[0]} {acell[0]} }} \n'.format( acell=[1/Bohr] ) )

    fd.write( 'rprim {{ {cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
              '        {cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
              '        {cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}  }}\n'
                   ''.format(cell=atoms.cell))

    fd.close()


############

def main():

    print( "Starting ... " )
    
    oceanJSON = dict()

    jsonFile = os.path.join( os.getcwd(), 'ocean.json' )
    if os.path.exists( jsonFile ):
        print( "Found local ocean.json")
        with open( jsonFile, 'r' ) as fd:
            oceanJSON = json.load(fd)
    else:
        jsonFile = os.path.join( os.path.dirname(os.path.realpath(__file__)), 'ocean.json' )
        if os.path.exists( jsonFile ):
            print( "Found ocean.json with this script")
            with open( jsonFile, 'r' ) as fd:
                oceanJSON = json.load(fd)
        else:
            print( "No ocean.json found" )

    inputFile = input("Name of input file?\n")

    atoms = read( inputFile )
    


    makePrimitive = input("Reduce to a primitive cell? [y/n]\n")
    if makePrimitive == 'y' or makePrimitive == 'Y':
        print( "Reducing to primitive" )
        makePrimitive = True
    else:
        print( "Not reducing to primitive" )
        makePrimitive = False



    new_cell = standardize_cell( (atoms.cell, atoms.get_scaled_positions(), atoms.numbers ), to_primitive=makePrimitive, symprec=5e-3 )

    new_unit_cell, new_scaled_positions, new_numbers = new_cell
    new_atoms = Atoms(new_numbers, cell=new_unit_cell, scaled_positions=new_scaled_positions)



    edgeElement = input("Select element for XAS: " + ", ".join(set(new_atoms.symbols))  + "\n" ) 
    
    edgeElementNumber = 0 
    for el in set(new_atoms.symbols):
        if el == edgeElement:
            edgeElementNumber = Atoms( el ).numbers[0]
            print( "Selected edge element {0} atomic number {1}".format( el, edgeElementNumber) )
        else:
            print( "Didn't select {0}".format( el) )

    if edgeElementNumber != 0:
        nl = input("Specify N and L for the x-ray edge (e.g. '1 0' for K or '2 1' for L23)\n")

        oceanJSON["edges"] = "-{0} {1}".format( edgeElementNumber, nl )
        
    
    try:
        write_ocean_in( "ocean.in", new_atoms, input_data=oceanJSON )
    except:
        raise Exception("FAILED while trying to write ocean.in")


if __name__ == '__main__':
    main()
