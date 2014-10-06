#!/bin/bash
#
# Author: Heri
#
# Script de submitere a job-urilor pe fiecare coda, folosind compilatoare diferite
#

mprun.sh --job-name MyTestGcc-O --queue ibm-opteron.q \
	--modules "compilers/gcc-4.6.0" \
	--script exec_script.sh --show-qsub --show-script --batch-job
mprun.sh --job-name MyTestGcc-N --queue ibm-nehalem.q \
	--modules "compilers/gcc-4.6.0" \
	--script exec_script2.sh --show-qsub --show-script --batch-job
mprun.sh --job-name MpiTestGcc-Q --queue ibm-quad.q \
	--modules "compilers/gcc-4.6.0" \
	--script exec_script3.sh --show-qsub --show-script --batch-job
sleep 5
a=`qstat | wc -l`
echo nr de joburi este "$a"
while [[  "$a" -ne "0" ]]; do
	let a=`qstat | wc -l`
done
if [[ $a -eq 0 ]]; then
	echo apelez diff pe fisierele output
	diff optim.out blas.out > dif_opt_blas.out
	diff mana.out blas.out > dif_dmana_blas.out
	gnuplot file.gp
fi

