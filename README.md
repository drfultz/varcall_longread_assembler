# tagged_longread_assembler
A snakemake pipeline for assembly of long DNA reads in difficult to assemble regions

**Introduction**\
Long-read sequencing (Oxford Nanopore MinION, Pacific Biosciences SMRT) has vastly improved the ability to assemble complete genomes.  However, even the longest reported reads from these methods (which extend into megabases) are insufficient to span repeat regions commonly found within eukaryotic genomes.  While there are often polymorphisms among repeat copies that could serve as landmarks for assembly, the inaccuracy of long-read methods limits their usefullness.  This poses a challenge for existing DNA assembly software, which relies primarily on per-base level of identity between two reads.  However, these regions may have structural variants or other landmarks across repeat copies that contain sufficient information, even within noisy reads, to allow assembly.  This pipeline uses primarily existing tools to generate and polish a consensus sequence using long reads from  difficult to assemble regions where unique features can be identified.

**Requirements**\
This pipeline is designed to be run by the workflow manager Snakemake, using its built-in interface with the package manager Conda.  The easiest way to achieve this is to install Conda, then use Conda to install Snakemake.  All other dependecies will be installed by Conda when Snakemake is run with the `--use-conda` option.

**Getting started**\
Clone this git repository.  Included in the repository is a test input file and a matching final output that is empty.\
Change directory to the repository folder.\
Run snakemake with the test file as the snakemake target and `--use-conda`.  This will automatically install all other dependencies and perform a relatively short test assembly.
```
git clone https://github.com/drfultz/tagged_longread_assembler.git
cd tagged_longread_assembler
snakemake -j 1 --use-conda polished/test_polished2X.fa
```
It will likely take several minutes for conda to install all of the dependencies on the initial run.  Upon completion, the file `polished/test_polished2X.fa` should be rewritten and contain a single fasta entry.

**Workflow**\
How it works

**Performing your own run**\
1) The first step for creating an assembly is to find reads from the region.  To do this, define sequence features (or patterns of features) that are unique to the region.  For example, one set of tandem repeats of the 45S ribosome genes in *Arabidopsis thaliana* are bordered on one end by the At4g00005 protein coding gene.  Reads containing this junction are obviously unique to this region, as this gene is unique.  Reads containing this feature can be retrieved by BLAST or minimap2 alignment to known sequence.

Less obvious, might be unique structural variants within repeats.  Alignments to a consensus repeat, dot plots against a consensus repeat, and/or automated feature annotations can be useful for identifying such variants.  For example, in a set of *A. thaliana* 45S repeats, annotating repeat reads revealed a region where two tandem repeat units were both truncated in the same manner.  While this type of truncation was seen to occur in several places, it appeared (from examining nearby features) that this tandem truncation occured in only one location in the genome.

2) 
