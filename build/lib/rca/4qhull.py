#!/usr/bin/env python

import glob
import rca

fname = 'pvdz-charge.txt'

f = open(fname, 'w')

for i in glob.glob('M*.g03'):
	typ, func = i.split('.')[0].split('-')
	mthp = rca.Cwrap(i)
	atomnos = mthp.get_atoms()
	atomcoords = mthp.get_pos()
	atomcharge = mthp.get_charge()

	f.write('%s\n' % func)
	f.write('%s\n' % str(atomnos))
	f.write('%s\n' % str(atomcoords))
	f.write('%s\n' % str(atomcharge))


f.close()
