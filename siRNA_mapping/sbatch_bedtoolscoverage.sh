#!/bin/bash
#SBATCH -J bedtoolscov
#SBATCH -o bedtoolscov.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bedtoolscov.e%A.%a   # error file name
#SBATCH -a 1-90
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 1
###SBATCH --mem-per-cpu=1750
#SBATCH -p gcluster,p1     # queue (partition) -- normal, development, largemem, etc.
###SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin  # email me when the job starts
###SBATCH --mail-type=end    # email me when the job finishes

###gfffiles from 1-34
MFELINENUM=$1

PATH=/gpool/bin/bedtools2/bin:$PATH

ALIGNJOBFILES=alignmentfiles.txt
MFEJOBFILES=mfewindowfiles.txt
ALIGNFILE=$(head -n ${SLURM_ARRAY_TASK_ID} ${ALIGNJOBFILES} | tail -n 1)
MFEFILE=$(head -n ${MFELINENUM} ${MFEJOBFILES} | tail -n 1)
###[gprmnsSlt]*.final.sorted.bed.gz
COVOUTFILE=$(basename ${MFEFILE} .sorted.bed)_$(basename ${ALIGNFILE} .sorted.bed).coverage

echo "Begin on ${MFEFILE} and ${ALIGNFILE}"
bedtools coverage -a ${MFEFILE} -b ${ALIGNFILE} -d -sorted > covfiles/${COVOUTFILE}
echo "Finished"
cd covfiles
gzip ${COVOUTFILE}
