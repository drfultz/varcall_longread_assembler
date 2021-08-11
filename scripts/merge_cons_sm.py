#!/usr/bin/env python3
import sys
import pysam

###attempts to find a tandem duplication of each input tag.
###If found, collapses into single tag.

input_tags_fa = snakemake.input.tag_fa
input_seq_fastas = snakemake.input.cons
output_fasta_file = snakemake.output[0]

###process and sort the input consensus filenames
input_seq_fastas_list = str(input_seq_fastas).split()
input_seq_fastas_list.sort()

###concatenate the sequences from each file in order
in_seq_joined = ""
#loop through the sequences except the final one, trimming the last 10000 bp of each
for input_seq_fasta in input_seq_fastas_list[:-1]:
	input_entry = list(pysam.FastxFile(input_seq_fasta))[0]
	in_seq_joined = in_seq_joined + input_entry.sequence[:-10000] #ought to change to remove based on length of used tags

#Add final sequence, since it needs to retain final tag
last_entry = list(pysam.FastxFile(input_seq_fastas_list[-1]))[0]
in_seq_joined = in_seq_joined + last_entry.sequence

#Probably not necessary now - was from using sequence matching in old script
in_seq_upper = in_seq_joined.upper() #changing all to uppercase eliminates issues where the read/tag junction is misaligned,
					#resulting in lowercase tag bases from cons command

###write to file
with open(output_fasta_file, "w") as output_fasta:
	output_fasta.write(f">{snakemake.wildcards.contig}\n{in_seq_upper}")
