#!/bin/env python
# Make all these imports only apply to the sections that need them (make CLI faster)
import os, re
from ase.calculators.jacapo import Jacapo
from cclib.parser import Gaussian
#from pylab import *
import logging
import geom
from ase.gui.view import View
import numpy
from ase.atoms import Atoms

"""
Additional Functions Needed:
    Molecular mass
    ThetaROT
    ZPE
    Sum of Elect + ZPE
    Gibbs
    Forces
    CPU Time
    Atomic Charge
    Complexity (PW Cutoff or NBasis)
"""


def atom_nos(tp, func, txt=None):
    """ txt = comp, func, or ncn """
    if txt is not None:
        f = open('/home/walesi/data/gaussian/tests/%s.txt' % func, 'r')
        info = f.read().split('\r')
        dict = {}
        for line in info:
            entry = line.split('\t')
            dict["".join(entry[0:1])] = list(entry[-1])
        return dict

    else:
        return None

class Cwrap():
    """
    This class represents a universal wrapper for Gaussian and Dacapo Calculations.  With data from this class
    and by parsing information from the file name, we can compare data and get all pertinent data from python
    """

    def __init__(self, filename):
        """
        Constructor.

        filename -- filename to get info from, either Gaussian (out, g03) or netCDF file (nc)
        """
        self.file = filename
        if os.path.exists(filename):
            ext = filename.split('.')[-1]
            self.complex = filename.split('.')[-2]
            #print ext
            if ext == 'nc':
                self.type = 'nc'
                # read directly from nc to avoid calculation
                #self.data = netCDF(filename,'r')
                self.data = Jacapo.read_atoms(filename)
            elif ext == 'out' or ext == 'g03' or ext == 'log':
                gout = Gaussian(filename)
                gout.logger.setLevel(logging.ERROR)
                self.type = 'gau'
                self.data = gout.parse()
        else:
            raise Exception, 'File is not of identifiable format'

    def get_pos(self):
        """
        get the positions of all the atoms
        """
        if self.type == 'nc':
            positions = self.data.get_positions()
        elif self.type == 'gau':
            positions = self.data.atomcoords[-1]
        return positions

    def get_atm(self):
        if self.type == 'nc':
            atom = self.data.get_atomic_numbers()
        elif self.type == 'gau':
            atom = self.data.atomnos
        return atom

    def get_nrg(self):
        if self.type == 'nc':
            pot = self.data.get_potential_energy()
        elif self.type == 'gau':
            pot = self.data.scfenergies[-1]
        return pot

    def get_vib(self):
        if self.type == 'gau':
            # this correction comes from Lang 96 for all DFT calculations in Gaussian
            correction = 0.9613
            vibfrq = []
            for i in self.data.vibfreqs[-len(self.data.vibfreqs)/2:]:
                vibfrq.append(i*correction)
        return vibfrq

    def get_cal(self):
        if self.type == 'nc':
            comp = self.data.get_cutoff()
        elif self.type == 'gau':
            comp = self.data.nbasis
        return comp

    def get_chg(self):
        if self.type == 'gau':
            f = open(self.file)
            f = f.read()
            regexp = re.compile('(?:APT atomic charges:\s+\d+\s+)[0-9\.\-A-Za-z\s]*')
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
        else:
            charges = [0 for i in self.get_natoms()]
        return charges

    def get_frc(self):
        if self.type == 'nc':
            force = self.data.get_forces()
        if self.type == 'gau':
            f = open(self.file)
            f = f.read()
            regexp = re.compile('(?:Forces[A-Za-z\s\(\)\,\/]*[-]{67})[\s0-9\.\-]*')
            force = [i.split()[2:] for i in regexp.findall(f)[-1].split('\n')[3:-2]]
        return force

    def get_na(self):
        if self.type == 'nc':
            comp = self.data.get_number_of_atoms()
        elif self.type == 'gau':
            comp = self.data.natom
        return comp

    def view(self):
        if self.type == 'nc':
            View(self.data)
        elif self.type == 'gau':
            View(Atoms(numbers = self.data.atomnos, positions = self.data.atomcoords[-1]))

    def get_dis(self, a1, a2):
        return False

    def get_ang(self, a1, a2, a3):
        """
        calculate angle from three points
        """
        l1 = self.get_bond_length(a1,a2)
        l2 = self.get_bond_length(a2,a3)
        angle = numpy.arctan(l2/l1)
        return geom_.degrees(angle)

