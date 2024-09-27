#!/bin/bash
#SBATCH -J bioawk
#SBATCH -o bioawk.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bioawk.e%A.%a   # error file name
#SBATCH -a 1-28              #
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

###bioawk
PATH=/gpool/bin/bioawk:$PATH

JOBFILE=trimmedfiles.txt
i=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
OUTDIR=siRNA_filtered_reads

mkdir -p ${OUTDIR}
bioawk -cfastx '{if (length($seq) == 21) print $name"\n"$seq}' $i > ${OUTDIR}/$(basename $i .fastq.gz)_21.fasta
bioawk -cfastx '{if (length($seq) == 22) print $name"\n"$seq}' $i > ${OUTDIR}/$(basename $i .fastq.gz)_22.fasta
bioawk -cfastx '{if (length($seq) == 24) print $name"\n"$seq}' $i > ${OUTDIR}/$(basename $i .fastq.gz)_24.fasta
