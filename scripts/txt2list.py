#!/usr/bin/env python

"""
Takes information output to txt file, processes it and feeds it into lists for
further use in
"""

import sys, re, os
import geom

def volume():
	fname = 'pvtz-charge.txt'

	g = open(fname, 'r')
	f = g.read()
	g.close()
	mthp = re.compile(r'[A-Zlr23]+\n')
	atoms = re.compile(r'\[[\s0-9]+\]')
	coords = re.compile(r'\[{2}[\[\]\s\.\-\+0-9\n]+\]{2}')

	mthp = [i.strip('\n') for i in mthp.findall(f)]

	atoms = [i.strip('][').split() for i in atoms.findall(f)]
	atm = []
	for atom in atoms:
		atm.append([int(i) for i in atom])
	coords = [re.sub('[\]\'\[]' ,'',i).split('\n') for i in coords.findall(f)]

	cds = []
	for items in coords:
		rw = []
		for rows in items:
			rw.append([float(i) for i in rows.split()])
		cds.append(rw)

	for n, mlc in enumerate(mthp):
		if mlc == 'C2H3':
			print cds[n]
		fname = 'tz-%s-vol.txt' % mlc

		if not os.path.exists(fname):

			indices = input('type in functional groups of %s:' % mlc)
			if isinstance(indices, tuple):
				indices = [i-1 for i in indices]
			else:
				indices -= 1

			if indices != '':
				geom.make_3d_space(atm[n],cds[n],indices=indices,fname=fname)
			else:
				geom.make_3d_space(atm[n],cds[n],fname=fname)

def charge(typ):
	fname = 'pvtz-charge.txt'

	g = open(fname, 'r')
	f = g.read()
	g.close()
	mthp = re.compile(r'[A-Zlr23]+\n')
	atoms = re.compile(r'\[[\s0-9]+\]')
	charge = re.compile(r'\[\'[\-0-9\.\'\,\ ]+')

	mthp = [i.strip('\n') for i in mthp.findall(f)]

	atoms = [i.strip('][').split() for i in atoms.findall(f)]
	atm = []
	for atom in atoms:
		atm.append([int(i) for i in atom])

	chg = []
	for ch in charge.findall(f):
		ch = re.sub(r'[\,\[\']','',ch)
		chg.append([float(i) for i in ch.split(' ')])

	fname = 'tz-%s-APT.txt' % typ
	f = open(fname,'w')
	for n, mlc in enumerate(mthp):
		cur_chg = chg[n]
		print cur_chg
		indices = input('%s : type in functional groups of %s:' % (typ, mlc))
		if isinstance(indices, tuple):
			indices = [i-1 for i in indices]
			sum = 0
			for i in indices:
				sum += cur_chg[i]
		else:
			indices -= 1
			sum = cur_chg[indices]
		f.write('%s : %0.3f\n' % (mlc, sum))
	f.close()

if __name__ == "__main__":
	volume()
	charge('func')
	charge('ncn')
