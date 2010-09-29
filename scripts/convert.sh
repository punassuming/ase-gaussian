#!/bin/sh

for i in `find ./ -name '*.g03' -print`
do
	python ~/dev/dacapo/ncprep.py $i
done
