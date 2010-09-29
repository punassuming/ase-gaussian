 #!/bin/csh

#The following should contain your program and any arguments
set inputfile  = "./ZWIT-DBN.cp.6-31.gjf"		# G03 Input file

set outputfile = "ZWIT-DBN.cp.6-31.log"		# G03 Output file

if ($?PBS_JOBID) then		# if this is a PBS job
	echo "Starting" $PBS_JOBID `date`
	echo "Initiated on `hostname`"
	echo ""
	cd "$PBS_O_WORKDIR" || exit 1		# connect to working directory of qsub
endif

if ($?PBS_NODEFILE) then
	#count the number of processors assigned by PBS
	set NP = `wc -l < $PBS_NODEFILE `
endif

g03 < $inputfile >& $outputfile



# Cluster Verbose Batch Script
#PBS -N n2opt.verbose
#PBS -j oe
#PBS -m ae
#PBS -l cput=0:51:00
#PBS -l walltime=0:60:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/csh
set echo 
setenv
df
set WORK=$HOME/cc/g98/examples
cd $TMPDIR
pwd
set INPUT=n2opt.verbose.com
cp $WORK/$INPUT .
module load g98a_9
g98 < ./$INPUT
ls -al
cp * $WORK

# Define input and output files.

set input=
set output=`echo $input | sed "s/gjf/log/"`

# check to see if output already exists
if [ -f output]: then
    # determine if file is already optimized
    set jstage=`grep -c 'Normal' $output`
    if [$jstage -e 2]: then
        echo "Frequency Calculation Already Completed"
        echo "Ending Job"
        exit 0
    else if [$jstage -e 1]: then
        echo "Optimization Already Completed"
        echo "Continuing on to frequency iteration"
        # Change job line to restart frequency calculation
        sed "s/' freq '/' freq=restart '/g" "$input"
        # Remove all optimazation parameters from job line
        sed "s/opt[=0-9a-zA-Z,()]*//g"