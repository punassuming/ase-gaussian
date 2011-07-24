#!/usr/bin/env python

"""
Workflow

*create Atoms object*

from xyz
>>> atoms = xyz_to_atoms(xyz_file)

from logfile
>>> gc = Gaussian('water')
>>> atoms = gc.get_atoms()

*sp calculation*
>>> gc = Gaussian('single_point',job = {'SP':[]}, atoms=atoms)

*using header file*
>>> gc = Gaussian('using_header',header = '#P Opt=tight
                   Pop=(chelpg,hirschfield)'

from checkpoint ?
 TODO use izmat to obtain input string from checkpoint
"""

__docformat__ = 'restructuredtext'

import exceptions, glob, os, pickle, string
import numpy as np
import subprocess as sp
from cclib.parser import Gaussian as gparse
from ase import Atoms
import logging

try:
    import cPickle as pickle
except:
    import pickle

log = logging.getLogger('Gaussian')

def coords_from_checkpoint(filename):
    if os.path.exists(filename):
        pass

def coords_from_log(filename)
    if os.path.exists(filename):
        pass

def xyz_to_atoms(filename)
    if os.path.exists(filename):
        pass

def job_from_log(filename):
    pass

def read(basename)
    """ 
    Return atoms object from following inputs:

    chk, fchk, gjf, log, gc

    atoms = read('water.chk')
    """
    calc = Gaussian(basename)
    atoms = calc.get_atoms()
    return atoms

class Gaussian():
    """
    An ASE calculator associated with Gaussian
    """

    __name__ = 'Gaussian'
    __version__ = '0.01'

    # additional parameters are passed to kwargs for input file creation
    def __init__(self,
                 basename = 'g09_test',
                 output = None,
                 debug=logging.WARN,
                 **kwargs)
        
        default_jobs = {
            'Density':[],
            'Force':[],
            'Freq':[],
            'Geom':[],
            'Guess':[],
            'Massage':[],
            'NMR':[],
            'Opt':[],
            'Polar':[],
            'Pop':[],
            'SCRF':[],
            'SP':[],
            'Scan':[],
            'Volume':[]
        }

        default_kwrgs = {
            'atoms':None,
            'xc':'B3LYP',
            'basis':'6-31+',
            'job':{'SP':[]},
            'multiplicity':1,
            'charge':0,
            'mem':0.5,
            'header':None
        }

        """
        Parameters that are defined
        
            basename : string
             basename used for gaussian log and checkpoint files and binary
             extensions looked for:
                 gc              : gaussian ASE object as binary file
                 log | g03 | g09 : gaussian output log file
                 chk | fchk      : gaussian checkpoint file

            output : string
             different basename of output file to allow copying old job,
             will revert to basename if not chosen

            debug : integer
             logging debug level

            kwargs : **
             includes all major features of gaussian

                atoms : ASE.Atoms object
                 atoms object is attached to calculator to save information
                
                xc : text
                 checked against possible gassian parameters
                 setup currently only for DFT jobs:
                     ['B3LYP',

                basis : text
                 checked against possible gaussian parameters

                job : dict
                 contains optional keyword arguments to place in gaussian input,
                 keys are available as the Gaussian job options:

                    ['Density','Force','Freq','Geom','Guess','Massage','NMR',
                     'Opt','Polar','Pop','SCRF','SP','Scan','Volume']

                 ** Opt and Freq can be together, but not with SCRF **

                 e.g. a standard optimization will be {'opt':[]}
                 * a more complex string would be:
                    {'opt':['tight'],'freq':[],'pop':['Hirshfield','NBO']}

                multiplicity : integer
                
                charge : integer

                mem : float
                 value in GB to reserve for memory

                header : string 
                 plain text header to overwrite entire job string, 
                 ** disables the use of many methods **

                header : string 
                 plain text header to overwrite entire job string, 
                 ** disables the use of many methods **

            Instead of an input file, which is typically used for gaussian
            jobs, as ASE atoms object with calculator information is saved for
            easier job creation.  This information contains all the job
            submission information in addition to the atom positions.  In
            order to obtain the energies, the log and checkpoint files are
            parsed to determine if the job has already been run and has
            completed successfully.  If so, the energy is returned, otherwise,
            the input string is constructed based on the atoms object and the
            job options specified on init.  This string is passed to gaussian
            through stdin, preventing the need for temporary input files. 

            Examples
            
             requires water file to already exist
            >>> gau = Gaussian('water')

             restart previous optimization
            >>> gau = Gaussian('CO2', job = {'opt':['restart']})

        """

        self.debug = debug
        log.setLevel(debug)
        
        # set up parameters to look for
        #self.pars = Gaussian.default_kwargs.copy()
        #self.pars_uptodate = {}

        #log.debug(self.pars)
        
        #for key in self.pars:
            #self.pars_uptodate[key] = False


        # setup file associations, self.chk, self.log, self.gc
        self.set_files(basename)

        # default to False and switch if job complete
        self.ready = False


        # set up kwargs to change / add 
        self.kwargs = Gaussian.default_kwargs.copy()

        # if gc file exists, load kwargs from file
        if os.path.exists(self.gc):
            # load atoms and calculator information from binary file 
            self.load_gc(self, self.gc)

        self.command_kwargs = kwargs

        # default is overwritten by gc which is overwritten by 

    def load_gc(self, gcfile):
        """ loads atoms and pertinent settings to calculator """

            gc = open(gcfile, 'rb')

            self.gc_kwargs = pickle.load(gc)

            if 'atoms' not in old_kwargs:
                log.debug('no atoms object was found in binary file')
                gc.close()
                raise Exception, 'binary file exists but does not contain atoms'

    def save_gc(self, gcfile):
        """ save atoms, calc and pertinent settings to calculator """

        
    def set_files(self, basename):
        """ set up filenames for binary, log and checkpoint files """

        gc = basename + '.gc'

        # if not setup
        if not hasattr(self, 'gc'):
            self.gc = gc

        # if different than previous and doesnt exist, copy from old
        if gc != self.gc and not os.path.exist(gc):
            log.debug('copying %s to %s' % (self.gc, gc))

            if os.path.exists(self.gc):
                status = os.system('cp %s %s' % (self.gc, gc))
            self.gc = gc

        # if exists, don't set it
        elif os.path.exists(gc):
            # get atoms object from gc
            self.atoms = self.read_atoms(gc)
            self.gc = gc

        self.chk = basename + '.chk'

        for possible_log in ['.log','.g03','.g09']:
            if os.path.exists(basename + possible_log):
                self.log = basename + possible_log

        if not hasattr(self, 'log'):
            self.log = basename+'.g09'

    def get_atoms(self, gc):
        """ return the atoms object associated with a calculator """

        if hasattr(self, 'atoms'):
            if self.atoms is None: 
                return None
            atoms = self.atoms.copy()
            atoms.set_calculator(self)
        else:
            atoms = None
        return atoms

    def format_job_line(self):
        pass

    def format_input_string(self):
        pass



"""
TODO Additional Functions Needed:
    Molecular mass
    ThetaROT
    ZPE
    Sum of Elect + ZPE
    Gibbs
    CPU Time
"""

class Parser(gparse):
    """
    This class represents a wrapper for cclib's gaussian parser
    """

    def __init__(self, logfile):
        """
        filename -- filename to get info from Gaussian logfile 
        """
        if os.path.exists(logfile):
            super().__init__(logfile)
            self.parse()
            
            self.file = open(logfile,'rb').read()


    def get_pos(self):
        """ get the positions of all the atoms """
        return self.data.atomcoords[-1]

    def get_atom(self):
        return self.data.atomnos

    def get_natoms(self):
        return len(self.data.atomnos)

    def get_potential_energy(self):
        return self.data.scfenergies[-1]

    def get_frequencies(self):
        # this correction comes from Lang 96 for all DFT calculations in gparse
        correction = 0.9613
        vibfrq = []
        for i in self.data.vibfreqs[-len(self.data.vibfreqs)/2:]:
            vibfrq.append(i*correction)
        return vibfrq

    def get_basis(self):
        return self.data.nbasis

    def get_charge(self, method):
        #if method in self.charge_methods.keys():
            #regexp = re.compile(self.charge_methods[method])
        if hasattr(self,'get_chg_'+method):
            getchg = 'self.get_chg_%s()' % method
            return eval(getchg)

        def get_chg_apt(self):
            f = self.file
            regexp = re.compile('(?:APT atomic charges:\s+\d+\s+)[0-9\.\-A-Za-z\s]*')
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
            return charges

        def get_chg_chelpg(self):
            f = self.file
            #regexp = re.compile('(?:Charge=\s+[0-9\.]+\sDipole=+\s+[ 0-9\.\-]+Tot=[0-9\.]+\s+\d+\s+)[0-9\.\-A-Za-z\s]*')
            regexp = re.compile('(?:Charge=\s+[0-9\.\ \-A-Z a-z\=]*Tot=\s+)[0-9\.\-A-Za-z\s]*')
            #print regexp.findall(f)[-1].split('\n')[2:-3]
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-2]]
            return charges

        def get_chg_mulliken(self):
            f = self.file
            regexp = re.compile('(?:Mulliken atomic charges:\s+\d+\s+)[0-9\.\-A-Za-z\s]*')
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
            return charges

        def get_chg_nbo(self):
            f = self.file
            regexp = re.compile('(?:Core\s+Valence\s+Rydberg\s+Total\s+[-]+)[0-9\.\-A-Za-z\s]*')
            charges = [i.split()[2] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
            return charges

    def get_forces(self):
        f = self.file
        regexp = re.compile('(?:Forces[A-Za-z\s\(\)\,\/]*[-]{67})[\s0-9\.\-]*')
        forces = [i.split()[2:] for i in regexp.findall(f)[-1].split('\n')[3:-2]]
        return forces

    def make_atoms(self):
        from ase import Atoms
        self.atoms = Atoms(numbers = self.data.atomnos, positions = self.data.atomcoords[-1]))

def calculate_distance(self, a1, a2):
    return False

def calculate_angle(self, a1, a2, a3):
    """
    calculate angle from three points
    """
    l1 = self.get_bond_length(a1,a2)
    l2 = self.get_bond_length(a2,a3)
    angle = numpy.arctan(l2/l1)
    return geom.degrees(angle)

