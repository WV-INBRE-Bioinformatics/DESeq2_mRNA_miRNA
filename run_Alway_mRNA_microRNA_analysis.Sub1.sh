#!/bin/bash
# This script submits the R script (RSCRIPT) to the cluster
# Needs an argument (VAR) that represents a result variable

RSCRIPT=$1
VAR=$2
PWD=`pwd`

if [ -f $RSCRIPT ]; then
	source MODULES
	mkdir -p $PWD/$VAR

	echo "Working on $VAR"
	echo "#!/bin/bash" > $PWD/$VAR/qsub_script.$VAR.sh
	echo "cd $PWD/$VAR" >> $PWD/$VAR/qsub_script.$VAR.sh
	echo "Rscript $RSCRIPT $VAR" >> $PWD/$VAR/qsub_script.$VAR.sh
	qsub -q batch -V -N $VAR -d $PWD/$VAR $PWD/$VAR/qsub_script.$VAR.sh
else
	echo "Not working on $VAR. $RSCRIPT not found!"
fi
