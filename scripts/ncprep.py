#!/bin/env python

"""
Pull information from gaussian output (chk) and save as netcdf file

#TODO : make connections between gaussian file from cclib import to ase's ListOfAtoms

Author: Rich Alesi

"""
import os, sys
from Jacapo import Jacapo
#from ase import *
from cclib.parser import Gaussian
import logging
import numpy as np
from ase.atoms import Atoms
from ase.gui import view

def convertStr(s):
    """Convert string to either int or float."""
    try:
        ret = int(s)
    except ValueError:
        #Try float.
        try:
            ret = float(s)
        except ValueError:
            ret = s
    return ret

def gau_get(in_file):

    gout = Gaussian(in_file)
    gout.logger.setLevel(logging.ERROR)
    gdata = gout.parse()

    gatoms = Atoms(numbers = gdata.atomnos, positions = gdata.atomcoords[-1])
    return gatoms

def xyz_get(in_file):

    f = open(in_file)
    atoms = f.read().split('\n')
    # assumes a format of int float float float
    els = []
    pos = []
    for lines in atoms:
        sp_line = lines.split()
        print sp_line
        sp_line = [convertStr(a) for a in sp_line]
        print sp_line
        if len(sp_line) >3 and type(sp_line[0]) is int and type(sp_line[3]) is float:
            els.append(int(sp_line[0]))
            pos.append([float(a) for a in sp_line[1:4]])
    xyzatoms = Atoms(numbers = els, positions = pos)
    return xyzatoms
"""
main function that accepts command line arguments (files)
"""

if __name__ == "__main__":

    in_file = str(sys.argv[1])
    if os.path.exists(in_file):
        # Filename is of form mlc-func.basis.g03
        base, ext = os.path.splitext(in_file)
        headPath, tail = os.path.splitext(base)

        gauNC = '%s.nc' % (headPath)
        if ext == '.out' or ext == '.g03':
            gatoms = gau_get(in_file)
        elif ext == '.xyz' or ext == '.gjf' or ext == '.com':
            gatoms = xyz_get(in_file)
        else:
            print "File not recognized."
            exit()
    else:
        print "File not found."
        exit()

    # Determine distances between atoms to set cell size
    pos = gatoms.get_positions()
    nn = []

    for i, atom in enumerate(gatoms):
        for j in range(i, len(gatoms)):
            dist = np.ndarray.tolist(np.abs(pos[i]-pos[j]))
            nn.append(dist)

    # Find maximum distance in all dimensions and factor by 2
    xx, yy, zz = zip(*nn)
    cellSize = [2*(max(xx)+2), 2*(max(yy)+2), 2*(max(zz)+2)]

    gatoms.set_cell(cellSize)
    gatoms.center()

    # Now view to check cell
    #view(gatoms)

    # Now set up the Jacapo calculator for gatoms

    calc = Jacapo(atoms = gatoms,
            pw = 300,
            kpts = (1,1,1),
            nc = gauNC,
            ft = 0.01)

    gatoms.set_calculator(calc)

    view(gatoms)

    nBands = int(calc.get_valence() * 0.65 + 4)

    calc.set_nbands(nBands)

    calc.write_nc(nc = gauNC)

    print 'Conversion complete.'
