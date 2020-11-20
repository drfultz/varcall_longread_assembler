#!/usr/bin/env python3
import sys
import pysam

###Finds, then removes tags

# input_tags_fa = sys.argv[1]
# input_tagged_fa = sys.argv[2]
# output_seq_file = sys.argv[2].split("_")[0] + "_contig.fa"

input_tags_fa = snakemake.input[0]
input_tagged_fa = snakemake.input[1]
output_seq_file = snakemake.output[0]

tagged_read = list(pysam.FastxFile(input_tagged_fa))[0]
tagged_seq_upper = tagged_read.sequence.upper()

###iterate the input tag sequences from the tags fasta
for tag in pysam.FastxFile(input_tags_fa):
	tag_upper = tag.sequence.upper()
	if tag_upper in tagged_seq_upper: #if a perfect copy of the tag sequence is in the target sequence
		tagged_seq_upper = tagged_seq_upper.replace(tag_upper,"")
		if tag_upper in tagged_seq_upper: #if a perfect copy of the tag sequence is STILL in the target sequence
			raise Exception(f"Tag {tag.name} occured more than once.  Please review the tag deduplication step")
		else:
			sys.stderr.write(f"Tag {tag.name} removed.\n")
	else:
		sys.stderr.write(f"Tag {tag.name} not found.\n")

with open(output_seq_file, "w") as output_seq:
	output_seq.write(f">{tagged_read.name}\n")
	output_seq.write(tagged_seq_upper)