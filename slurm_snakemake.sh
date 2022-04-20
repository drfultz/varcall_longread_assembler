#!/bin/bash
#SBATCH -J sm_master
#SBATCH -p general
#SBATCH -o sm_master_%j.out
#SBATCH -e sm_master_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d@iu.edu
#SBATCH -D /N/slate/d/tagged_longread_assembler
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=36:00:00

TARGET=final_contigs/test-tel_contig.fa

source ~/.bashrc
conda activate snakemake

snakemake --cluster 'sbatch \
-J sm_sub \
-p general \
-o logs/{wildcards.contig}_{rule}_${{SLURM_JOBID}}.out \
-e logs/{wildcards.contig}_{rule}_${{SLURM_JOBID}}.err \
--mail-type=ALL \
--mail-user=d@iu.edu \
--nodes=1 \
--ntasks=1 \
--cpus-per-task={params.ppn} \
--mem-per-cpu={params.vmem} \
--time={params.walltime} \
' --jobs 100 --local-cores 1 --use-conda $TARGET
