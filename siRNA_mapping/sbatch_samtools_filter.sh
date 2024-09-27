#!/bin/bash
#SBATCH -J samtoolsfilter
#SBATCH -o samtoolsfilter.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e samtoolsfilter.e%A.%a   # error file name
#SBATCH -a 1-90                #
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

###samtools 1.10
PATH=/gpool/bin/samtools-1.10/bin:$PATH

JOBFILE=bamfiles.txt
BAMFILE=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

samtools view -b -f 4 ${BAMFILE} > $(basename ${BAMFILE} .bam).unmapped.bam
samtools view -b -F 4 ${BAMFILE} > $(basename ${BAMFILE} .bam).mapped.bam
