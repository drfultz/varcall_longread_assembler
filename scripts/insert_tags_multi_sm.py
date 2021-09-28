#!/usr/bin/env python3
import sys
import re
import csv
import pysam
import glob

# features_file = sys.argv[1]
# oriented_fastq_file = sys.argv[2]
# search_patterns_file = sys.argv[3]
# input_tags_file = sys.argv[4]
# tagged_reads_file = "%s_with_%s.fq" % (oriented_fastq_file.split('.', 1)[0], search_patterns_file.split('.', 1)[0])
# used_readnames_file = "%s_usedreadnames.txt" % (oriented_fastq_file.split('.', 1)[0])

features_file = snakemake.input[0]
oriented_fastq_file = snakemake.input[1]
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
if (len(patterns) - len(spacings)) != 0:
	raise Exception(f"The number of search patterns {len(patterns)} (first line) should equal the number of spacings {len(spacings)} (second line)")
#print(patterns)
#print(spacings)

##### Create a fresh used_readnames_file for this contig
with open(used_readnames_file, 'w') as used_readnames_file_obj:
	pass

##### gather the used readnames
used_readnames = set()
for other_used_readnames_file in glob.glob("used_readnames/*_usedreadnames.txt"):
	with open(other_used_readnames_file, 'r') as other_used_readnames:
		for line in other_used_readnames:
			used_readnames.add(line.strip())
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

##### Open and loop through the features file and fastq for matching the patterns and inserting tags
with open(features_file, 'r') as feature_file_obj, open(tagged_reads_file,'w') as tagged_reads, open(used_readnames_file, 'w') as used_readnames_file_obj:

	### open the features file and input tags
	feature_reader = csv.reader(feature_file_obj, delimiter='\t', quoting=csv.QUOTE_NONE)
	input_tags = list(pysam.FastxFile(input_tags_file))

	### loop through the features file and oriented fastq
	for feature_line,fastq_entry in zip(feature_reader,pysam.FastxFile(oriented_fastq_file)):
		offset = 0

		if feature_line[0] == fastq_entry.name:
			if feature_line[0] not in used_readnames:
				### loop through the patterns
				for pattern,spacing in zip(patterns,spacings):

					### increase the offset for later patterns
					offset += int(spacing)

					### find all occurences of the pattern
					matched_feature = re.findall(pattern, feature_line[1])

					### if one match
					if len(matched_feature) == 1:
						print(f"read {feature_line[0]} matched with pattern at unit {offset}")

						## Determine the number of units Left of the pattern (num_L_units)
						## in order to the set the 'frame' of tags to insert
						split_features = re.split(pattern, feature_line[1]) #splits the features at the search pattern
						L_features = split_features[0]
						#print(L_features)
						num_L_units = len(L_features.split(",")) - 1

						## Prepare stuff
						num_units = len(feature_line[1].split(","))
						feature_starts = feature_line[2].split(",")
						working_seq = fastq_entry.sequence
						working_qual = fastq_entry.quality

						## Loop backwards (so as not to disturb downstream insertion sites)
						## through each units, adding tags
						for n in range(num_units-1, 0, -1):
							start_tag_num = n-num_L_units+10+offset
							start_site = int(feature_starts[n])
							#print(f"\tinserting tag #{start_tag_num}, named {input_tags[start_tag_num].name}, into pos {start_site}")
							working_seq, working_qual = tagit(input_tags[start_tag_num], start_site, working_seq, working_qual)

						## Write the tagged read output
						output_read = pysam.FastxRecord()
						output_read.name = f"{fastq_entry.name}_tagged"
						output_read.sequence = working_seq
						output_read.quality = working_qual
						tagged_reads.write(str(output_read) + "\n")

						used_readnames_file_obj.write(fastq_entry.name + "\n")

						### break looping through further patterns after one is found
						break

					### if more than one match, make note
					elif len(matched_feature) > 1:
						print (f"\nPattern\n{pattern}\nwas present more than once in read {feature_line[0]}\n")
						#raise Exception(f"\nPattern\n{pattern}\nwas present more than once in a single read\n")
		else:
			raise Exception(f"\n\tThe .feat file readname {feature_line[0]} does not match the paired entry in the oriented fastq {fastq_entry.name}.\n")
