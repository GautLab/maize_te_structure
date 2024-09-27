#!/bin/bash
#SBATCH -J bedtoolsb2b              # jobname
#SBATCH -o bedtoolsb2b.o%A.%a       # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bedtoolsb2b.e%A.%a       # error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-90                 # 12,14,15      # 1-11,16,17                    # start and stop of the array start-end
###SBATCH -n 1                  # -n, --ntasks=INT Maximum number of tasks. Use for requesting a whole node. env var SLURM_NTASKS
#SBATCH -c 1                    # -c, --cpus-per-task=INT The # of cpus/task. env var for threads is SLURM_CPUS_PER_TASK
#SBATCH -p p1,gcluster                 # queue (partition) -- normal, development, largemem, etc.
###SBATCH --mem-per-cpu=1750
###SBATCH -t 72:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
###SBATCH --mail-user=solarese@uci.edu
###SBATCH --mail-type=begin     # email me when the job starts
###SBATCH --mail-type=end       # email me when the job finishes

## SLURM ENV Variables
#  SLURM_CPUS_ON_NODE
#  SLURM_CPUS_PER_TASK
#  SLURM_ARRAY_TASK_ID

###bedtools 2.27.1
PATH=/gpool/bin/bedtools2-2.27.1/bin:$PATH
PATH=/gpool/bin/bedtools2-2.27.1/scripts:$PATH

JOBFILE=bamfiles.txt
BAM=$(basename $(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1) .bam).sorted.bam

bedtools bamtobed -i ${BAM} > $(basename ${BAM} .sorted.bam).bed
