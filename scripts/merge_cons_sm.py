#!/usr/bin/env python3
import sys
import pysam

#####Trims the 3' tags off of each consensus sequence and concatenates them.

input_tags_fa = snakemake.input.tag_fa
input_seq_fastas = snakemake.input.cons
output_fasta_file = snakemake.output[0]


###measure the length of the input tags.  This assumes all tags are of the same length.
with open(input_tags_fa, 'r') as input_tags_fa_obj:
	second_line = input_tags_fa_obj.readlines()[1]
	tag_length = len(second_line.split()[0])
	print(tag_length)

###process and sort the input consensus filenames
input_seq_fastas_list = str(input_seq_fastas).split()
input_seq_fastas_list.sort()

###concatenate the sequences from each file in order
in_seq_joined = ""
#loop through the sequences except the final one, trimming the last tag_length bases
for input_seq_fasta in input_seq_fastas_list[:-1]:
	input_entry = list(pysam.FastxFile(input_seq_fasta))[0]
	in_seq_joined = in_seq_joined + input_entry.sequence[:-tag_length]

#Add final sequence, since it needs to retain final tag
last_entry = list(pysam.FastxFile(input_seq_fastas_list[-1]))[0]
in_seq_joined = in_seq_joined + last_entry.sequence

#Probably not necessary now - was from using sequence matching in old script
in_seq_upper = in_seq_joined.upper() #changing all to uppercase eliminates issues where the read/tag junction is misaligned,
					#resulting in lowercase tag bases from cons command

###write to file
with open(output_fasta_file, "w") as output_fasta:
	output_fasta.write(f">{snakemake.wildcards.contig}\n{in_seq_upper}")
