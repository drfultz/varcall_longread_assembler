#!/usr/bin/env python3
import pysam
import sys
import subprocess
import re
import os

#####First, pair the 'sequence tags' together as search patterns.
###Takes a 'tag' fastx and creates a search pattern from each successive pair of sequences.
###The pattern matches any sequence between an exact match of the first sequence and
###an exact match of the second sequence

#####Then, look for that pattern 
###A second fastx is iterated one time for each search pattern.
###The first occurence of the pattern per read (inclusive of the searched tag sequences)
###is saved into a fasta.  No file is written if no match is found.

###This program is intended to find unique 'tags' artificially placed into noisy long reads
###at the same known positions within each read. The intent is to generate MSA inputs bounded by
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


count = 0
in_interior = False
out_interior = False

#iterate through each tag, paired with the subsequent tag
for tagA,tagB in zip(tags[:-1:],tags[1::]):
	#print some tag identifiers
	print(f"TAG{str(count).zfill(3)} to TAG{str(count+1).zfill(3)}")
	#print(f"{tagA[0:5]}   {tagB[0:5]}")
	
	#Define the search parameter as any alphanumeric between the exact tags (inclusive of tags)
	search_string = f"{tagA}\\w+{tagB}"
	content = []

	#iterate through reads for this tag
	for read in pysam.FastxFile(input_reads_fa):
		try:
			subseq = re.search(search_string, read.sequence).group(0)
			#if match found, append a new fasta entry with the matched string
			print(f"\t{read.name[0:34]} - {len(subseq)}")
			content.append(f">{read.name}")
			content.append(subseq)
		except AttributeError:
			pass
			#print(f"{read.name[0:34]} - Not found")

	#if more than one matches are found for the pattern
	if len(content) > 3 and out_interior == False:
		#make a directory using snakemake wildcard
		path = "segment_reads/" + snakemake.wildcards["contig"]
		os.makedirs(path, exist_ok=True)
		#write the matched fasta entries for this tag
		out_fa_file = f"{path}/seg_{str(count).zfill(3)}to{str(count+1).zfill(3)}_reads.fa"
		with open(out_fa_file, "w") as out_fa:
			out_fa.write("\n".join(content) + "\n")
		in_interior = True

	elif len(content) < 3 and in_interior == True and out_interior == False:
		out_interior = True
		print("\n\t1st coverage drop\n")

	elif len(content) > 3 and out_interior == True:
		raise Exception("\nCoverage for a region between two tags dropped below 2 then went back up again, indicating a hole in this assembly.\n")

	#increase the counter
	count += 1
