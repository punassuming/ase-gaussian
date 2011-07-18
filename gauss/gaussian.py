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

            output : string
             different name of output file if the user desires, 
             will revert to basename if not chosen

            debug : integer
             logging debug level

            kwargs : **
             contains many possible items, listed below

                xc : text
                 checked against possible gassian parameters

                basis : text
                 checked against possible gaussian parameters

                keywords : dict
                 contains optional keyword arguments to place in gaussian input

                mem : options to specify 
        
        
        
        
        Asign filenames derived from basename
        self.basename = basename
        self.log = basename + '.log'
        self.chk = basename + '.chk'
        self.pick = basename + '.gc'


