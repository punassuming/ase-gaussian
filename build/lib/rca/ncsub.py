#!/usr/bin/env python
"""
Generic Jacapo submission script

Author: Rich Alesi

"""
from Jacapo import *
from ase import *
import torque
import sys, os, re

# Accept only one input, the original molecule path (form: MLC-FUNC.nc)

fileName = str(sys.argv[1])

base, ext = os.path.splitext(fileName)
init =[]
energies =[]
cutOffs = [300, 350, 400, 450, 500]

for cutOff in cutOffs:
	ncFile = '%s.%i.nc' % (base, cutOff)

	# Check to see if file exists
	if os.path.exists(ncFile):
		nccalc = Jacapo(ncFile, pw=cutOff)
	else:
		nccalc = Jacapo(fileName, outnc=ncFile, pw=cutOff)

	ncatoms = nccalc.atoms

	ncatoms.set_calculator(nccalc)

	init.append(ncatoms.get_potential_energy)

	from ase.optimize.qn import *
	dyn = QuasiNewton(ncatoms)
	dyn.run(fmax=0.05)

	energies.append(ncatoms.get_potential_energy())
from pylab import *
plot(cutOffs,energies,'bo-')
plot(cutOffs,init,'ro')
xlabel('cutoff energies (eV)')
ylabel('energies (eV)')
savefig('%s.png' % (base))
