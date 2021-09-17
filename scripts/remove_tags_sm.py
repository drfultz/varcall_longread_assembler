#!/usr/bin/env python3
import sys
import re
import pysam

###Finds, then removes tags

# input_tags_fa = sys.argv[1]
# input_tagged_fa = sys.argv[2]
# output_seq_file = sys.argv[2].split("_")[0] + "_contig.fa"

input_tags_fa = snakemake.input[0]
input_tagged_fa = snakemake.input[1]
output_seq_file = snakemake.output[0]

tagged_contig = list(pysam.FastxFile(input_tagged_fa))[0]
tagged_contig_upper = tagged_contig.sequence.upper()

###iterate the input tag sequences from the tags fasta
for tag in pysam.FastxFile(input_tags_fa):
	partial_tag_upper = tag.sequence.upper()[4:-4]
	
	#if partial copy of the inner tag sequence is in the target sequence
	if partial_tag_upper in tagged_contig_upper:
		
		#remove the partial tag and four bases on either side
		tagRE = r"\S{4}" + str(partial_tag_upper) + r"\S{4}"
		tagged_contig_upper = re.sub(tagRE, "", tagged_contig_upper, count=0)
		
		#if a partial copy of the tag sequence is STILL in the target sequence, exception
		if partial_tag_upper in tagged_contig_upper:
			raise Exception(f"Tag {tag.name} occured more than once.  Please review the tag deduplication step")
		else:
			sys.stderr.write(f"Tag {tag.name} removed.\n")
	else:
		sys.stderr.write(f"Tag {tag.name} not found.\n")

with open(output_seq_file, "w") as output_seq:
	output_seq.write(f">{tagged_contig.name}\n{tagged_contig_upper}\n")
