#!/bin/bash
#SBATCH -J btbuild
#SBATCH -o btbuild.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e btbuild.e%A.%a   # error file name
#SBATCH -a 1                #
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

REF="Zea_mays.AGPv4.dna.toplevel.fa"
OPTIONS="--threads ${SLURM_CPUS_PER_TASK}"

bowtie2-build ${OPTIONS} ${REF} $(basename ${REF} .fa)
