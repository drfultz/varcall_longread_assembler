# tagged_longread_assembler
A snakemake pipeline for assembly of long DNA reads in difficult to assemble regions

**Introduction**\
Long-read sequencing (Oxford Nanopore MinION, Pacific Biosciences SMRT) has vastly improved the ability to assemble complete genomes.  However, even the longest reported reads from these methods (which extend into megabases) are insufficient to span repeat regions commonly found within eukaryotic genomes.  While there are often polymorphisms among repeat copies that could serve as landmarks for assembly, the inaccuracy of long-read methods limits their usefullness.  This poses a challenge for existing DNA assembly software, which relies primarily on per-base level of identity between two reads.  However, these regions may have structural variants or other landmarks across repeat copies that contain sufficient information, even within noisy reads, to allow assembly.  This pipeline uses primarily existing tools to generate and polish a consensus sequence using long reads from  difficult to assemble regions.  Unique features must be identified beforehand, and unique 'sequence tags' are inserted into the reads to act as landmarks for assembly.

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

**Preparing your own inputs**\
Preparing your inputs as fasta or fastq reads containing sequence tags is the most time consuming step and requires careful work.  Essentially this step is outlining where each read should be assembled by hand for the pipeline to perform later. A GUI program for working with DNA sequence is necessary for the following steps to be practical.  Geneious Prime works very well for this, but free alternatives like Genome Compiler, Unipro UGENE, or Benchling's sequence tools would likely be sufficient.

1\) The first step for creating an assembly is to find reads from the region.  To do this, define sequence features (or patterns of features) that are unique to the region.  For example, one set of tandem repeats of the 45S ribosome genes in *Arabidopsis thaliana* are bordered on one end by the At4g00005 protein coding gene.  Reads containing this junction are clearly unique to this region, as this gene is unique.  Reads containing this feature can be retrieved by BLAST or minimap2 alignment to the unique sequence.

Less obvious unique features might include unique structural variants within repeats.  Alignments to a consensus repeat, dot plots against a consensus repeat, and/or automated feature annotations can be useful for identifying such variants.  For example, in a set of *A. thaliana* 45S repeats, annotating repeat reads revealed a region where two tandem repeat units were both truncated in the same manner.  While this type of truncation was seen to occur in several places, it appeared (from examining nearby features) that this tandem truncation occured in only one location in the genome.

2\) Reverse-complement any reads that are in the reverse direction (relative to the desired contig).

3\) Now that reads are gathered and unique features are defined, plan where to insert the sequence tags.  (A set of randomly genereated 10kb tags is included in tags/10kb_tags.fasta)  The tags need not be limited to the location of unique features but should be placed relative to where unique features are located.  For example, one tag could be located directly upstream of a feature, while another could be located tens of kb downstream.

>Using the example of the 45S gene/At4g0005 gene above, tags were placed at the start of repeats.  This was done at every fourth repeat unit after the unique junction.  Tag A was placed at the 4th repeat, tag B at the 8th repeat, and tag C at the 12th repeat, etc.

Each tag should be inserted in the position of any read that represents a specific genomic coordinate.  The tags must also be placed in order from left to right according to the expected contig.  

4\) Using Geneious to perform the actual edits to the reads can be performed thus (I believe that Genome Compiler has very similar features):
- A) Open all of the reads in Geneious together (as a Geneious sequence list)
- B) Use the '**Annotate from database**' tool to mark a relatively invariable region of your repeat on all of the reads (this annotation need not be long, 200bp works).  If possible, also annotate your unique feature(s).
- C) Copy and paste each tag into the appropriate region on each read one tag at a time.  Place the tag precisely by clicking on the annotation, then pressing the left arrow key followed by right arrow key to position the cursor to exactly before the annotation.  Press **CTRL+V** to paste.  Once this is done for every tag, save the reads are ready.

5\) Export the Geneious sequence list as a fasta named `{sample_name}_tagreads.fa` and place it in the `tagged_longread_assembler/tagged_reads` folder.

6\) Make sure that the fasta file containg your sequence tags is listed as the `Tag_set` value in the `config.yaml` file.

7\) From the `tagged_longread_assembler` folder, run `snakemake -j 1 --use-conda {sample_name}_tagreads.fa`
