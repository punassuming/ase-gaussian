#!/usr/bin/env python

import glob
import rca

fname = 'pvtz-vib.txt'

f = open(fname, 'w')

for i in glob.glob('*.g03'):
	x = i.split('.')[0].split('-')
	if len(x) == 1:
		typ = str(x)
		func = ''
	else:
		typ = x[0]
		func = x[1]

	mthp = rca.Cwrap(i)
	atomnrg = mthp.get_nrg()
	atomfreq = mthp.get_vib()

	f.write('%s\t%s\n' % (typ,func))
	f.write('%0.6f\n' % atomnrg)

	for i in atomfreq:
		for j in i:
			f.write('%0.3f ' % j)
		f.write(',')
	f.write('\n\n')


f.close()
