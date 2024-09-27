#!/bin/bash
#SBATCH -J genucounts
#SBATCH -o genucounts.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e genucounts.e%A.%a   # error file name
#SBATCH -a 7-90                #
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 1
###SBATCH --mem-per-cpu=1750
#SBATCH -p gcluster,p1     # queue (partition) -- normal, development, largemem, etc.
###SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin  # email me when the job starts
###SBATCH --mail-type=end    # email me when the job finishes

## SLURM ENV Variables
#  SLURM_CPUS_ON_NODE
#  SLURM_CPUS_PER_TASK
#  SLURM_ARRAY_TASK_ID

JOBFILE=ucounts.txt
i=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

awk 'NR%2 == 0' ${i} | sort | uniq -c | awk '{print ">"$2"_"$1"\n"$2}' > $(basename ${i} .fasta).ucounts.fasta
gzip ${i}
gzip $(basename ${i} .fasta).ucounts.fasta
