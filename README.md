# tagged_longread_assembler
A snakemake pipeline for assembly of long DNA reads in difficult to assemble regions

**Introduction**
Long-read sequencing (Oxford Nanopore MinION, Pacific Biosciences SMRT) has vastly improved the ability to assemble complete genomes.  However, even the longest reported reads from these methods (which extend into megabases) are insufficient to span repeat regions commonly found within eukaryotic genomes.  While there are often polymorphisms among repeat copies that could serve as landmarks for assembly, the inaccuracy of long-read methods limits their usefullness.  This poses a challenge for existing DNA assembly software, which relies primarily on per-base level of identity between two reads.  However, these regions may have structural variants or other landmarks across repeat copies that contain sufficient information, even within noisy reads, to allow assembly.  This pipeline uses primarily existing tools to generate and polish a consensus sequence using long reads from  difficult to assemble regions where unique features can be identified.

**Installation**
