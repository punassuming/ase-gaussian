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
>>> gc = Gaussian('single_point',route = {'SP':[]}, atoms=atoms)

*using header file*
>>> gc = Gaussian('using_header',header = '#P Opt=tight
                   Pop=(chelpg,hirschfield)')

from checkpoint ?
 TODO use izmat to obtain input string from checkpoint
"""


__docformat__ = 'restructuredtext'

import exceptions, re, os, pickle, string
import numpy as np
import subprocess as sp
import cclib
from ase import Atoms
import logging
import io

try:
    import cPickle as pickle
except:
    import pickle

log = logging.getLogger('Gaussian')

def route_from_log(filename):
    pass

def read(basename):
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

    default_kwargs = {
        'atoms':None,
        'xc':'B3LYP',
        'basis':'6-31+',
        'route':{'SP':[]},
        'multiplicity':1,
        'charge':0,
        'mem':0.5,
        'desc':""
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

                route : dict
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
                 plain text header to overwrite entire route string, 
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
            >>> gau = Gaussian('CO2', route = {'opt':['restart']})

"""
    # additional parameters are passed to kwargs for input file creation
    def __init__(self,
                 basename = 'g09_test',
                 output = None,
                 debug=logging.WARN,
                 status = 'New',
                 **kwargs):

        self.debug = debug
        log.setLevel(debug)
        
        # setup file associations, self.chk, self.log, self.gc
        self.set_files(basename)

        self.ready = False 
        self.status = status
        
        # set up kwargs base initial call of Gaussian object
        self.params = Gaussian.default_kwargs.copy()
        self.kwargs = kwargs

        # if gc file exists, load kwargs from file
        if os.path.exists(self.gc):
            # load gaussian object parameters from .gc binary file and should have to do nothing else
            #assuming we want all of the parameters that are currently written in .gc file
            #could possible update this to check if anythign needs to be updated 
            try:
                self.load_gc(self.gc)
                if 'atoms' in kwargs:
                    atoms = kwargs['atoms']
                    atoms.set_positions(self.atoms.get_positions())
                    atoms.calc = self
            except:
                raise Exception, 'should delete .gc file and run again '
            
        else:
            #read in kwargs from existing file 
            #else if .log file exists, then define gausian object from log file and make .gc file
            if os.path.exists(self.log):
                self.read_kwargs, stat = io.read_log(self.log)
                if stat == 'Finished':
                    self.status = stat
            #elif os.path.exists(self.chk):
            #elif os.path.exists(self.gjf):
            #if no files exist, that just use what was called with the ojbect 
            else:
                kwargs['atoms'].calc = self
                self.read_kwargs = kwargs

            #check to see if any inputed kwargs in object call are not in read_kwargs
            #will overide read in kwargs with object call kwargs 
            for key in kwargs:
                if key in self.read_kwargs and (self.read_kwargs[key] <> kwargs[key]):
                    self.read_kwargs[key] = kwargs[key]
                elif key in self.default_kwargs:
                    self.read_kwargs.update({key:kwargs[key]})
                else:
                    raise Exception, 'not a valid parameter'
            
            if len(self.read_kwargs)>0:
                if 'atoms' in self.read_kwargs:
                    self.set_atoms(self.read_kwargs['atoms'])
                    del self.read_kwargs['atoms']
                self.set(self.read_kwargs)


            self.save_gc(self.gc)

        print self.atoms.calc
        self.route = self.params['route']

        if self.status == 'Finished':
            self.cclib_parser = cclib.parser.Gaussian(self.log).parse()
            self.file = open(self.log, 'rb').read()
        
            
    def set(self, newkwargs):
        for key in newkwargs:
            if key not in self.default_kwargs:
                raise Exception, key + "not valid kwarg"
            
            if self.params[key]<> newkwargs[key]:
                self.params[key]=newkwargs[key]
                
    def set_atoms(self, atoms):
        if hasattr(self, 'atoms') and self.atoms is not None:
            if self.atoms_are_equal(atoms):
                return
        self.atoms=atoms.copy()
        self.atoms.calc = self
        self.params['atoms'] = self.atoms
                
    def atoms_are_equal(self, atoms):
        '''
        comparison of atoms to self.atoms using tolerances to account
        for float/double differences and float math.
        '''
        
        TOL = 1.0e-6 #angstroms

        a = self.atoms.arrays
        b = atoms.arrays

        #match number of atoms in cell
        lenmatch = len(atoms) == len(self.atoms)
        if lenmatch is not True:
            return False #the next two comparisons fail in this case.
        
        #match positions in cell
        posmatch = (abs(a['positions'] - b['positions']) <= TOL).all()
        #match cell
        cellmatch = (abs(self.atoms.get_cell()
                         - atoms.get_cell()) <= TOL).all()
        
        if lenmatch and posmatch and cellmatch:
            return True
        else:
            return False


    def set_status(self, file):
        pass
    
    def load_gc(self, gcfile):
        """ loads atoms and pertinent settings to calculator """

        gc = open(gcfile, 'rb')
        self.gc_kwargs = pickle.load(gc)

        if 'atoms' not in self.gc_kwargs:
            log.debug('no atoms object was found in binary file')
            gc.close()
            raise Exception, 'binary file exists but does not contain atoms'
        else:
            # read atoms only from gc file 
            self.params = self.gc_kwargs
            self.atoms = self.gc_kwargs['atoms'].copy()
            self.atoms.calc = self
            self.status =  self.gc_kwargs['atoms'].calc.status 

    def save_gc(self, gcfile):
        """ save atoms, calc and pertinent settings to calculator """
        gc = open(gcfile, 'wb')
	pickle.dump(self.params, gc)
        
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
                os.system('cp %s %s' % (self.gc, gc))
            self.gc = gc

        # if exists, don't set it
        elif os.path.exists(gc):
            # get atoms object from gc
            #self.atoms = self.read_atoms(gc)
            self.gc = gc

        self.chk = basename + '.chk'

        for possible_log in ['.log','.g03','.g09']:
            if os.path.exists(basename + possible_log):
                self.log = basename + possible_log

        if not hasattr(self, 'log'):
            self.log = basename+'.g09'
        
        if not hasattr(self, 'gjf'):
            self.gjf = basename+".gjf"

    def get_atoms(self):
        """ return the atoms object associated with a calculator """
        if hasattr(self, 'atoms'):
            if self.atoms is None: 
                return None
            atoms = self.atoms.copy()
            atoms.set_calculator(self)
        else:
            atoms = None
        return atoms

    def format_route(self):
        kwrd = ""
        for key in self.params['route']:
            options = ",".join(self.params['route'][key])
            kwrds = key + "=(" + options + ") "
        this_route = "#" + self.params['xc']+"/"+self.params['basis'] +" " + kwrds

        return this_route

    def format_input_string(self):
        pass

    #neeed to test this 
    def calculation_required(self, calc_def):
##        calc_key = calc_def.keys()[0]
##        route = self.kwargs['route']
        
##        if calc_key in route:# and self.kwargs['status']="Finsihed":
##            if route[calc_key].find(calc_def[calc_key][0]):
##                return False
##            else:
##                return True
##        elif calc_def.keys()[0]=='freq' and 'opt' in self.route:
##            if route['opt'].find('calcall')<>-1:
##                return False
##            else:
##                return True
##        else:
##                return True
        if self.status == 'Finished':
            pass
        elif self.status == 'New':
            return True
    
    def update_route(self, calc_def):
        calc_key = calc_def.keys()[0]

        if calc_key in self.kwargs['route']:
            self.kwargs['route'][calc_key]=self.kwargs['route'][calc_key]+calc_def[calc_key]
        else:
            self.kwargs['route'][calc_key]=calc_def[calc_key]

        self.save_gc
    
    def get_potential_energy(self, atoms=None):
        if self.calculation_required(atoms):
                #log.debug('calculation required for energy'
                self.calculate()
        try:
            e = self.cclib_parser.scfenergies[-1]
            return e
        except:
            raise RuntimeError('Error in cauclating the total energy')

    #should I return this as a string or a floating point number?
    def get_ZPE(self):

        if self.calculation_require({'freq': ""}):
            self.update_route({'freq': ""})
            self.write_gjf()
            self.calculate()
        
        exprZPE = re.compile('Zero-point correction.*')
        try:
            ZPE_ht = float(exprZPE.findall(self.file)[-1].split()[2])
            #convert to hartree
            ZPE = str(ZPE_ht * 27.211396132)
        except:
            ZPE = 'No ZPE calculated'
        return ZPE

    def get_homo_energy(self):
        try:
            i_homo=self.get_homo_index
            homo_energy= self.get_orbital_energies(i_homo)
        except:
            homo_energy = 'cannot determine homo energy'
        return homo_energy
    
    def get_lumo_energy(self):
        try:
            i_lumo = self.get_homo_index + 1
            lumo_energy = self.get_orbital_energies(i_lumo)
        except:
            lumo_energy = 'cannot determine lumo energy'
        return lumo_energy

    def get_homo_index(self):
        return self.cclib_parser.homos[0]

    def get_orbital_energies(self):
        return self.cclib_parser.moenergies[0]

    def get_charge_density(self, method):
        #if method in self.charge_methods.keys():
            #regexp = re.compile(self.charge_methods[method])
        if hasattr(self,'get_chg_'+method):
            getchg = 'self.get_chg_%s()' % method
            return eval(getchg)

    def get_chg_apt(self):
        f = self.file
        regexp = re.compile('(?:APT atomic charges:\s+\d+\s+)[0-9\.\-A-Za-z\s]*')

        try: 
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
        except:
            charges = ['apt charge not calculated']

        return charges

    def get_chg_chelpg(self):
        f = self.file
        #regexp = re.compile('(?:Charge=\s+[0-9\.]+\sDipole=+\s+[ 0-9\.\-]+Tot=[0-9\.]+\s+\d+\s+)[0-9\.\-A-Za-z\s]*')
        regexp = re.compile('(?:Charge=\s+[0-9\.\ \-A-Z a-z\=]*Tot=\s+)[0-9\.\-A-Za-z\s]*')
        #print regexp.findall(f)[-1].split('\n')[2:-3]
        try:
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-2]]
        except:
            charges ['chelpg charges not calculated']
        return charges

    def get_chg_mulliken(self):
        f = self.file
        regexp = re.compile('(?:Mulliken atomic charges:\s+\d+\s+)[0-9\.\-A-Za-z\s]*')
        try:
            charges = [i.split()[-1] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
        except:
            charges = ['muliken charge not calculated']
        return charges

    def get_chg_nbo(self):
        f = self.file
        regexp = re.compile('(?:Core\s+Valence\s+Rydberg\s+Total\s+[-]+)[0-9\.\-A-Za-z\s]*')
        try:
            charges = [i.split()[2] for i in regexp.findall(f)[-1].split('\n')[2:-1]]
        except:
            charges = ['nbo charge not calcualted']
        return charges

    def get_chg_hirshfeld(self):
        f = self.file
        regexp= re.compile('(?:Hirshfeld spin densities).*[\0-9\-A-Z]*')
        try:
            charges = [i.split()[3] for i in regexp.findall(f)[-1].split("\n")[2:-1]]
        except:
            charges = ['hirshfeld charges not calculated']
        return charges

    def get_chg_bader(self):
    
        # need to develop way to name bader files with header name, will prob have to move .dat files to specific files names once calcualted
        
        pass
    
    def get_forces(self):
        f = self.file
        regexp = re.compile('(?:Forces[A-Za-z\s\(\)\,\/]*[-]{67})[\s0-9\.\-]*')
        forces = [i.split()[2:] for i in regexp.findall(f)[-1].split('\n')[3:-2]]
        return forces

    def set_basis(self, new_basis):
        if new_basis <> self.params['basis']:
            self.set({basis: new_basis, status:'New'})
            self.write_gc
        else:
            pass
        
    def get_basis(self):
        return self.params['basis']

    def set_xc(self, new_xc):
        if new_xc <> self.params['xc']:
            self.set({xc: new_xc, status: 'New'})
            self.write_gc
        else:
            pass
        
    def get_xc(self):
        return self.params['xc']

    def get_desc(self):
        return self.params['desc']
    
    def set_charge(self, new_charge):
        if new_charge <> self.params['charge']:
            self.set({charge: new_charge, status: 'New'})
            self.write_gc
        else:
            pass
        
    def get_charge(self):
        return self.params['charge']
    
    def set_multiplicity(self, new_mult):
        if new_mult <> self.params['multiplicity']:
            self.set({multiplicity: new_mult, status: 'New'})
            self.write_gc
        else:
            pass
        
    def get_multiplicity(self):
        return self.params['multip;licity']
    
     
    def calculate(self):
        #make .gjf file if necessary
        self.write_gjf()
        #run the calculation

        #if terminate, update status
        #if positions have changes, update atoms positions
        #make cclib parser, and test_data object 
        pass

    def write_gjf(self):
        gjf = open(self.gjf, 'wb')
        gjf.writelines("%chk="+self.chk+"\n")
        gjf.writelines("%mem="+str(self.params['mem'])+"\n")
        gjf.writelines("%nproc = 1"+"\n")
        gjf.writelines(self.format_route()+"\n")
        gjf.writelines(self.params['desc']+"\n")
        gjf.writelines(str(self.params['charge'])+str(self.params['multiplicity'])+"\n")

        mol_nos=self.atoms.get_atomic_numbers()
        mol_coords = self.atoms.get_positions()
        
        for i in range(len(mol_nos)):
            j = mol_coords[i]
            gjf.writelines(str(mol_nos[i]) + ' ' + str(j[0]) + ' '+ str(j[1])+ ' '+ str(j[2])+"\n")

        #insert end of file write for custom spec, like PCM model
        
        gjf.write("\n\n\n")
        gjf.close()
        
class GaussianParser(cclib.parser.gaussianparser.Gaussian):
    """
    This class represents a wrapper for cclib's gaussian parser
        TODO Additional Functions Needed:
        Molecular mass
        ThetaROT
        ZPE
        Sum of Elect + ZPE
        Gibbs
        CPU Time
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

    def make_atoms(self):
        from ase import Atoms
        self.atoms = Atoms(numbers = self.data.atomnos, positions = self.data.atomcoords[-1])

def calculate_distance(self, a1, a2):
    return False

def calculate_angle(self, a1, a2, a3):
    """
    calculate angle from three points
    """
    l1 = self.get_bond_length(a1,a2)
    l2 = self.get_bond_length(a2,a3)
    angle = np.arctan(l2/l1)
    return geom.degrees(angle)

