#/bin/bash
#PBS -N master_sm
#PBS -m abe
#PBS -e logs/master-sm_${PBS_JOBID}.err
#PBS -o logs/master-sm_${PBS_JOBID}.out
#PBS -M drfultz@iu.edu
#PBS -l nodes=1:ppn=1,vmem=8gb,walltime=36:00:00
#PBS -d /N/slate/drfultz/working_tagged_longread_assembler/

TARGET=final_contigs/NOR2-tel_contig.fa

source ~/.bashrc
conda activate snakemake

snakemake --cluster 'qsub -m abe -M drfultz@iu.edu -e logs/{wildcards.contig}_{rule}_${{PBS_JOBID}}.err -o logs/{wildcards.contig}_{rule}_${{PBS_JOBID}}.out -l nodes=1:ppn={params.ppn},vmem={params.vmem},walltime={params.walltime}' --jobs 100 --use-conda $TARGET
