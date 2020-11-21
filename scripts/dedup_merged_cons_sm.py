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
for input_seq_fasta in input_seq_fastas_list:
	input_entry = list(pysam.FastxFile(input_seq_fasta))[0]
	in_seq_joined = in_seq_joined + input_entry.sequence

in_seq_upper = in_seq_joined.upper() #changing all to uppercase eliminates issues where the read/tag junction is misaligned,
									#resulting in lowercase tag bases from cons command

###iterate the input tag sequences from the tags fasta
for tag in pysam.FastxFile(input_tags_fa):
	#if a perfect copy of the tag sequence is in the target sequence
	#look for a tandem duplication of that tag and replace it with a single copy
	tag_upper = tag.sequence.upper()
	if tag_upper in in_seq_upper:
		search_string = f"{tag_upper}{tag_upper}"
		if search_string in in_seq_upper:
			in_seq_upper = in_seq_upper.replace(search_string,tag_upper)
			#reports succseful deduplicaion.  This should occur where
			#segments of consensus sequence were merged.
			sys.stderr.write(f"Tag {tag.name} deduplicated.\n")
		else:
			#Reports a tag that is found, but not in a tandem duplication
			#For properly joined consensus sequence segements,
			#this should only occur at the first and last tag that
			#was present in the tagged reads
			sys.stderr.write(f"!!!Tag {tag.name} not joined!!!\n")

	else:
		#Reports tag not found. Should only occur at tags not included in
		#any read or not properly paired with an expected neighbor tag in the input reads
		sys.stderr.write(f"Tag {tag.name} not present.\n")

#write to file
with open(output_fasta_file, "w") as output_fasta:
	output_fasta.write(f">{snakemake.wildcards.sample}\n{in_seq_upper}")