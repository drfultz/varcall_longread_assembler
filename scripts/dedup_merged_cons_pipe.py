#!/usr/bin/env python3
import sys
import re

import pysam

###attempts to find a tandem duplication of each input tag.
###If found, collapses into single tag.

input_tags_fa = sys.argv[1]
input_seq = sys.stdin.read()

in_seq_upper = input_seq.upper() #changing all to uppercase eliminates issues where the read/tag junction is misaligned,
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

#write to stdout for piping
sys.stdout.write(in_seq_upper)