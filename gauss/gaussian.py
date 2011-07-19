#!/usr/bin/env python

"""


"""

import exceptions, glob, os, pickle, string
import numpy as np
import subprocess as sp

import logging

log = logging.getLogger('Gaussian')

def from_checkpoint(filename):
    if os.path.exists(filename):
        pass

def from_log(filename)
    if os.path.exists(filename):
        pass


class Gaussian():
    """
    An ASE calculator associated with Gaussian
    """

    __name__ = 'Gaussian'
    __version__ = '0.01'

    # additional parameters are passed to kwargs for input file creation
    def __init__(self,
                 basename = None,
                 output = None,
                 debug=logging.WARN,
                 **kwargs)
        
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

                 note: Opt and Freq can be together, but not with SCRF

                 e.g. a standard optimization will be {'opt':[]}
                 * a more complex string would be:
                    {'opt':['tight'],'freq':[],'pop':['Hirshfield','NBO']}

                mem : float
                 value in GB to reserve for memory

                debug : logging level used for parser and job creation

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

        #Asign filenames derived from basename
        self.basename = basename
        self.log = basename + '.log'
        self.chk = basename + '.chk'
        self.pick = basename + '.gc'


