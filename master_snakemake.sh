#/bin/bash
#PBS -N master_sm
#PBS -m abe
#PBS -e logs/${PBS_JOBID}_master_sm.err
#PBS -o logs/${PBS_JOBID}_master_sm.out
#PBS -M email@domain.edu
#PBS -l nodes=1:ppn=1,vmem=8gb,walltime=48:00:00
#PBS -d /location/of/Snakefile

TARGET=

snakemake --cluster 'qsub -m abe -M EMAIL@DOMAIN.edu -e logs/${{PBS_JOBID}}_{wildcards.sample}_{rule}.err -o logs/${{PBS_JOBID}}_{wildcards.sample}_{rule}.out -l nodes=1:ppn={params.ppn},vmem={params.vmem},walltime={params.walltime}' --jobs 100 --use-conda $TARGET
