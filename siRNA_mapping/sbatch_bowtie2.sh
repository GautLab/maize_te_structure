#!/bin/bash
#SBATCH -J btmatch
#SBATCH -o btmatch.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e btmatch.e%A.%a   # error file name
#SBATCH -a 1-90                #
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 16
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

###bowtie 2.4.1
PATH=/networkshare/bin/bowtie2-2.4.1-linux-x86_64:$PATH

JOBFILE=jobfile.txt
QRY=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
REF="Zea_mays.AGPv4.dna.toplevel.fa"
PREFIX=$(basename ${QRY} .fasta.gz)
PREFIX="maizeAGPv4.${PREFIX}"
OPTIONS="--threads ${SLURM_CPUS_PER_TASK} -a --score-min C,0,-1 -f"

bowtie2 $OPTIONS -x $(basename ${REF} .fa) -U ${QRY} -S ${PREFIX}.sam
