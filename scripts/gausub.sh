#!/bin/csh

#The following should contain your program and any arguments

set base = "BCAR-DBN.pvdz"
set inputfile  = "./$base.gjf"		# G03 Input file

set outputfile = "$base.log"		# G03 Output file

if ($?PBS_JOBID) then		# if this is a PBS job
	echo "Starting" $PBS_JOBID `date`
	echo "Initiated on `hostname`"
	echo ""
	cd "$PBS_O_WORKDIR" || exit 1		# connect to working directory of qsub
	mkdir /scratch/$PBS_JOBID
	/bin/cp `echo "$base.*"` "/scratch/$PBS_JOBID"
	cd "/scratch/$PBS_JOBID"
	setenv GAUSS_SCRDIR /scratch/$PBS_JOBID/g03_scr
	mkdir $GAUSS_SCRDIR
endif

g03 < $inputfile > $outputfile

set iter=1

while `grep -c 'Normal' $outputfile` != 2
	echo 'Restarting Optimization'
	sub ' opt ' ' opt=restart ' '$inputfile'
	sub ',connectivity' ',allcheck' '$inputfile'
	echo "iteration $iter" >> $PBS_O_WORKDIR/$base.out
	echo "steps completed `grep -c 'Steps' $outputfile`" >> $PBS_O_WORKDIR/$base.out
	cp *.log $PBS_O_WORKDIR
	iter=`expr $iter + 1`
	g03 < $inputfile > $outputfile
end

cp *.log $PBS_O_WORKDIR
fomchk *.chk
cp *.fchk $PBS_O_WORKDIR
cp *.chk $PBS_O_WORKDIR

echo "Removing scratch directory"
rm -R /scratch/$PBS_JOBID
