#!/bin/sh

for i in `find ./ -name *.nc`:
do
	/usr/matsim/pymatsim/Jacapo/tools/stripnetcdf $i
done
