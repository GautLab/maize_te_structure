#!/bin/bash
#SBATCH -J trimm
#SBATCH -o trimm.o%A.%a   # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e trimm.e%A.%a   # error file name
#SBATCH -a 1-28              #
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

#source /networkshare/.mybashrc
###java 1.8.221
PATH=/gpool/bin/jre1.8.0_221/bin:$PATH
export LD_LIBRARY_PATH=/gpool/bin/jre1.8.0_221/lib:$LD_LIBRARY_PATH
export JAVA_HOME=/gpool/bin/jre1.8.0_221/
export JRE_HOME=/gpool/bin/jre1.8.0_221/

#TRIMMOMATIC="/data/apps/trimmomatic/0.32/"
TRIMMOMATIC="/networkshare/bin/Trimmomatic-0.39"
WORKDIR="./"
READTYPE="SE"
JOBFILE=jobfile.txt
FASTQ=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)
SUFFIX=$(basename ${FASTQ} .fastq)

#TRIM="java -jar ${TRIMMOMATIC}/trimmomatic-0.32.jar"
TRIM="java -jar ${TRIMMOMATIC}/trimmomatic-0.39.jar"
READTYPE="SE"
OPTIONS1="${READTYPE} -threads ${SLURM_CPUS_PER_TASK} -phred33 -trimlog ${WORKDIR}${SUFFIX}.log"
INFILES="${WORKDIR}${SUFFIX}.fastq.gz"
OUTFILES="${WORKDIR}${SUFFIX}.trimm.fastq"
OPTIONS2="ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq-SE.fa:2:20:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:15 MINLEN:20 AVGQUAL:20"

$TRIM $OPTIONS1 $INFILES $OUTFILES $OPTIONS2

pigz -p ${SLURM_CPUS_PER_TASK} ${OUTFILES}
