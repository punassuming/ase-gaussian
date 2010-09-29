#!/bin/env python

"""
Script to make CP ready gaussian start file from old output file

Necessary format:

# cmd line (same as old) + Counterpoise=2 and no Freq command

0,1
1 0.00 0.00 0.92 1
9 0.17 0.00 2.73 2
1 0.77 0.00 2.73 2
9 0.00 0.00 0.00 1
"""

import sys
import os
from cclib.parser import Gaussian
import logging

for arg in sys.argv[1:]:
	gau_file = str(arg)

	if os.path.exists(gau_file):
		newname, ext = os.path.splitext(gau_file)
		cp_file = '%s-cp.com' % (newname)
		chk_file = '%s-cp.chk' % (newname)
	else:
		print 'File not found.'
		exit()

	gau_out = Gaussian(gau_file)
	gau_out.logger.setLevel(logging.ERROR)
	gau_data = gau_out.parse()

	g = open(cp_file,'w')

	header = '%%chk=%s\n%%mem=1900MB\n%%nproc=1\n# opt=(calcall) b3lyp/cc-pVDZ geom=angle Counterpoise=2\n\nCounterpoise calculation\n\n' % (chk_file)
	g.write(header)
	g.write('%i,%i\n' % (gau_data.charge,gau_data.mult))

	for num in range(0,gau_data.natom):
		atm = gau_data.atomnos[num]
		x, y, z = list(gau_data.atomcoords[-1][num])
		g.write('%i %0.5f %0.5f %0.5f 1\n' % (atm, x, y, z))

	g.close()
