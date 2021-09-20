#!/usr/bin/env python3
import sys
import re
import csv
import pysam
import glob

# features_file = sys.argv[1]
# input_reads_file = sys.argv[2]
# search_patterns_file = sys.argv[3]
# input_tags_file = sys.argv[4]
# tagged_reads_file = "%s_with_%s.fq" % (input_reads_file.split('.', 1)[0], search_patterns_file.split('.', 1)[0])
# used_readnames_file = "%s_usedreadnames.txt" % (input_reads_file.split('.', 1)[0])

features_file = snakemake.input[0]
input_reads_file = snakemake.input[1]
search_patterns_file = snakemake.input[2]
input_tags_file = snakemake.input[3]
used_readnames_file = snakemake.output[0]
tagged_reads_file = snakemake.output[1]


##### Parse the search patterns and distances
with open(search_patterns_file, 'r') as search_pattern_object:
	first_line = search_pattern_object.readline().strip()
	patterns = first_line.split('\t')
	second_line = search_pattern_object.readline().strip()
	spacings = second_line.split('\t')

if " " in first_line or " " in second_line:
	raise Exception("\nThere should not be space characters present in the _patterns.tsv file.\n")

if (len(patterns) - len(spacings)) != 0:
	print()
	raise Exception(f"\nThe number of search patterns {len(patterns)} (first line) should equal the number of spacings {len(spacings)} (second line).\n")
#print(patterns)
#print(spacings)

##### Create a fresh used_readnames_file for this contig
with open(used_readnames_file, 'w') as used_readnames_file_obj:
	pass

##### gather the used readnames
used_readnames = []
for other_used_readnames_file in glob.glob("used_readnames/*_usedreadnames.txt"):
	with open(other_used_readnames_file, 'r') as other_used_readnames:
		for line in other_used_readnames:
			used_readnames.append(line.strip())
#print(used_readnames)

##### Define tag insertion function for later
def tagit(tag, site, seq, qual):
	## insert tag into the sequence
	first_seq = seq[:site]
	second_seq = seq[site:]
	#print(f"  {len(seq)}")
	tagged_sequence = f"{first_seq}{tag.sequence}{second_seq}"
	#print(f"  s + {len(tag.sequence)} = {len(tagged_sequence)}")
	
	## insert dummy qualities into the qualities
	first_qual = qual[:site]
	second_qual = qual[site:]
	tag_qual = "I" * len(tag.sequence) # generate dummy quals for tag
	#print(f"  {len(qual)}")
	tagged_quality = f"{first_qual}{tag_qual}{second_qual}"
	#print(f"  q + {len(tag_qual)} = {len(tagged_quality)}")
	return tagged_sequence, tagged_quality

##### Open and loop through the features file for matching the patterns and inserting tags
with open(features_file, 'r') as feature_file_obj, open(tagged_reads_file,'w') as tagged_reads, open(used_readnames_file, 'w') as used_readnames_file_obj:

	### 'Mask' the readnames from the blacklist of used readnames
	feature_reader = csv.reader(feature_file_obj, delimiter='\t', quoting=csv.QUOTE_NONE)
	masked_feature_lines = []
	for feature_line in feature_reader:
		if feature_line[0] not in used_readnames:
			masked_feature_lines.append(feature_line)

	### loop through the features file with used readnames removed
	for mask_feature_line in masked_feature_lines:
		#print(f"\n\n!!!!!{mask_feature_line[0]}")
		offset = 0

		### loop through the patterns
		for pattern,spacing in zip(patterns,spacings):

			### increase the offset for later patterns
			offset += int(spacing)
			#print(f"searching for pattern with offset of {offset}")

			### find all occurences of the pattern
			matched_feature = re.findall(pattern, mask_feature_line[1])

			### if one match
			if len(matched_feature) == 1:
				print(f"read {mask_feature_line[0]} matched with pattern at offset {offset}")

				### loop through the fastq, looking for the first matching readname
				for input_read in pysam.FastxFile(input_reads_file):
					if mask_feature_line[0] == input_read.name: #find the read

						## Prepare stuff
						num_units = len(mask_feature_line[1].split(","))
						feature_starts = mask_feature_line[2].split(",")
						input_tags = list(pysam.FastxFile(input_tags_file))
						working_seq = input_read.sequence #!#!#! probably could just put the input_read object into the function
						working_qual = input_read.quality

						## Determine the number of units Left of the pattern (num_L_units)
						## in order to the set the 'frame' of tags to insert
						split_features = re.split(pattern, mask_feature_line[1]) #splits the features at the search pattern
						L_features = split_features[0]
						#print(L_features)
						num_L_units = len(L_features.split(",")) - 1

						#print(f"\tnumber units {num_units}")
						#print(f"\tleft units {num_L_units}")

						## Loop backwards (so as not to disturb downstream insertion sites)
						## through each units, adding tags
						for n in range(num_units-1, 0, -1):
							start_tag_num = n-num_L_units+10+offset
							start_site = int(feature_starts[n])
							#print(f"\tinserting tag #{start_tag_num}, named {input_tags[start_tag_num].name}, into pos {start_site}")
							working_seq, working_qual = tagit(input_tags[start_tag_num], start_site, working_seq, working_qual)

						## Write the tagged read output
						output_read = pysam.FastxRecord()
						output_read.name = f"{input_read.name}_tagged"
						output_read.sequence = working_seq
						output_read.quality = working_qual
						tagged_reads.write(str(output_read) + "\n")

						used_readnames_file_obj.write(input_read.name + "\n")

						## only first input read found
						break

				## if for loop completes without break from a match
				else: # no break
					print(f"\n\nfeature line {mask_feature_line[0]} not found in fastq\n")

				### break looping through further patterns after one is found
				break

			### if more than one match, exception
			elif len(matched_feature) > 1:
				print (f"\nPattern\n{pattern}\nwas present more than once in read {mask_feature_line[0]}\n")
				#raise Exception(f"\nPattern\n{pattern}\nwas present more than once in a single read\n")
