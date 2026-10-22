#!/usr/bin/env python3
import pysam

##### Concatenates the consensus sequences of the segments.

input_seq_fastas_list = list(snakemake.input.cons) # takes sorted order from snakemake
output_fasta_file = snakemake.output[0]

###concatenate the sequences from each file in order
in_seqs_list = []
#loop through the sequences
for input_seq_fasta in input_seq_fastas_list:
	in_seqs_list.append(list(pysam.FastxFile(input_seq_fasta))[0].sequence)

###write to file
with open(output_fasta_file, "w") as output_fasta:
	output_fasta.write(f">{snakemake.wildcards.contig}\n{''.join(in_seqs_list)}\n")
