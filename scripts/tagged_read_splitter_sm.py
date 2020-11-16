#!/usr/bin/env python3
import pysam
import sys
import subprocess
import re
import os

###First, pair the 'sequence tags' together as search patterns.
###Takes a first fastx and creates a search pattern from each successive pair of sequences.
###The pattern matches any sequence between an exact match of the each member of the pair.

###Then, look for that pattern 
###Then, a second fastx is iterated one time for each search pattern.
###The first occurence of the pattern per read (inclusive of the searched sequence)
###is saved into a fasta (one fasta per pattern).  No file is written if no match is found.

###This program is intended to find unique 'tags' artificially placed into noisy long reads
###at the same known positions within each read. The intent is to generate inputs for MSA bounded by
###these tags, which in turn will generate consensus sequences for segements of difficult to assemble
###regions of the genome

input_tags_fa = snakemake.input[0]
input_reads_fa = snakemake.input[1]

#input_tags_fa = sys.argv[1]
#input_reads_fa = sys.argv[2]

#gather sequence tags from the tags fasta
tags = []
for tag in pysam.FastxFile(input_tags_fa):
	tags.append(tag.sequence)


#iterate through each tag, paired with the subsequent tag
count = 0
for tagA,tagB in zip(tags[:-1:],tags[1::]):
	#print some tag identifiers
	print(f"\nTAG{count} to TAG{count+1}")
	print(f"{tagA[0:5]}   {tagB[0:5]}")
	
	#Define the search parameter as any alphanumeric between the exact tags (inclusive of tags)
	search_string = f"{tagA}\\w+{tagB}"
	content = ""

	#iterate through reads for this tag
	for read in pysam.FastxFile(input_reads_fa):
		print(read.name)
		try:
			subseq = re.search(search_string, read.sequence).group(0)
			#if match found, append a new fasta entry with the matched string
			print(len(subseq))
			content += (">" + read.name + "\n")
			content += (subseq + "\n")
		except AttributeError:
			print("not found")

	#if any matches are found for the pattern
	if content:
		#out_fa_file = "seg_" + str(count+1) + "to" + str(count) + ".fa"
		#make a directory using snakemake wildcard
		path = "segment_reads/" + snakemake.wildcards["sample"]
		os.makedirs(path, exist_ok=True)
		#write the matched fasta entries for this tag
		out_fa_file = f"{path}/seg_{str(count).zfill(2)}to{str(count+1).zfill(2)}_reads.fa"
		with open(out_fa_file, "w") as out_fa:
			out_fa.write(content)

	#increase the counter
	count += 1
