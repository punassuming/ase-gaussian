import os, re
from cclib.parser import Gaussian as cclibGaussian
import gaussian
from ase import atoms 

# I think we can read the chk, log and gjf files for kwargs to gaussian object
kwargs = {'atoms':None,
	'xc':'B3LYP',
	'basis':'6-31+',
	'route':{'SP':[]},
	'multiplicity':10,
	'charge':10,
	'mem':2,
	'desc':None
	}

def read_chk(filename):
    '''
    create gaussian object from object from checkpoint file
    '''
    if os.path.exists(filename):
        subprocess.call(['formchk ' + filename + " temp.fchk -3"],shell=True)
        text = open('temp.fchk', 'rb').read()
        return kwargs
    
    else:
        print '.chk file does not exst'


def read_log(filename):
    '''
    create gaussian object from gaussian log file
    '''
    if os.path.exists(filename):
		myfile = cclibGaussian(filename)
		try:
			data = myfile.parse()
			atom_numbers = data.atomnos
			atom_coords = data.atomcoords[-1]
			myatoms = atoms.Atoms(atom_numbers, atom_coords)
			kwargs['atoms']=myatoms
                        kwargs['charge']=data.charge
                        kwargs['multiplicity']=data.mult
		except:
			data = 'unable to open' + filename

                #need to write function to return status  
                kwargs['xc'], kwargs['basis'], kwargs['route'] = get_calcParams(filename)
                #actually write code to check status of calcualtion
                stat = 'Finished'
                return kwargs, stat
    else:
        print "log file does not exist" 

def get_calcParams(file):
    text = open(file, 'rb').read()

    basis = ""
    calcType = ""
    route = {}
    
    if(file[-4:]==".log"):
        expBasis = re.compile('Standard basis: .*')
        basis = expBasis.findall(text)[0].strip('Standard basis: ').split(' ')[0]
        expHeader = re.compile('(?:#).*[\d\0-9.\-A-Za-z\s]*')
        readHeader = expHeader.findall(text)[0]
        readHeader = readHeader[:readHeader.find("\n ---")].replace("\n","").replace(", ", ",").replace(' ,', ',')
        headerList = readHeader.lower().split(' ')

        for i in headerList:
            index = i.lower().find(basis.lower())
            if index <>-1:
                calcType = i[:index-1]
            else:
                if i.find("=")<>-1:
                    route.update({i[:i.find("=")].strip('#'): i[i.find("=")+1:].replace(')',"").replace("(","").strip().split(",")})
                else:
                    route.update({i.strip("#"):[]})
    elif(file[-4:]=='.chk'):
        basis = ""
        calcType=""
        route = {}

    elif(file[-4:]==".gjf"):
        basis = ""
        calcType =""
        route ={}

    return calcType, basis, route
        

