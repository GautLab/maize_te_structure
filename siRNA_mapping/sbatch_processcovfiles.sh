#!/bin/bash
#SBATCH -J proccovs
#SBATCH -o proccovs.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e proccovs.e%A.%a   # error file name
#SBATCH -a 1-3060               # 3060
###SBATCH -n 1                # Total number of mpi tasks requested
###SBATCH --ntasks-per-node 1
#SBATCH -c 1
###SBATCH --mem-per-cpu=1750
#SBATCH -p gcluster,p1     # queue (partition) -- normal, development, largemem, etc.
###SBATCH -t 48:00:00        # run time (hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin  # email me when the job starts
###SBATCH --mail-type=end    # email me when the job finishes

PATH=/gpool/bin/customscripts:${PATH}

JOBFILE=covjobfile.txt
#JOBFILE=failedjobs.txt
COVOUTFILE=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
#COVOUTFILE=${COVOUTFILE}.gz
OUTDIR=countfiles

processcovfiles.py -i ${COVOUTFILE} -o ${OUTDIR}
