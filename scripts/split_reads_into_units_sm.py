#!/usr/bin/env python3
import sys
import re
import csv
import pysam
import glob
import statistics
import os
from collections import defaultdict

# features_file = sys.argv[1]
# oriented_fastq_file = sys.argv[2]
# search_anchors_file = sys.argv[3]
# input_tags_file = sys.argv[4]
# tagged_reads_file = "%s_with_%s.fq" % (oriented_fastq_file.split('.', 1)[0], search_anchors_file.split('.', 1)[0])
# used_readnames_file = "%s_usedreadnames.txt" % (oriented_fastq_file.split('.', 1)[0])

features_file = snakemake.input[0]
oriented_fastq_file = snakemake.input[1]
search_anchors_file = snakemake.input[2]
q_cutoff = snakemake.params.Segment_q_cutoff
max_depth = snakemake.params.Max_depth
read_info_file = snakemake.output[0]
metadata_file = snakemake.output[2]

#####
##### Setup
#####

##### Parse the search anchors and distances
with open(search_anchors_file, 'r') as search_anchor_object:
	first_line = search_anchor_object.readline().strip()
	anchors = first_line.split('\t')
	second_line = search_anchor_object.readline().strip()
	spacings = second_line.split('\t')
if (len(anchors) - len(spacings)) != 0:
	raise Exception(f"The number of search anchors {len(anchors)} (first line) should equal"
		f" the number of spacings {len(spacings)} (second line)")

### Create a fresh read_info_file for this contig to clear it of used reads
with open(read_info_file, 'w') as read_info_file_obj:
	pass

### gather the other used readnames
used_readnames = set()
for other_used_readnames_file in glob.glob("split_read_info/*_read_info.tsv"):
	with open(other_used_readnames_file, 'r') as other_used_readnames:
		for line in other_used_readnames:
			used_readnames.add(line.strip().split("\t")[0])


#####
##### Use anchor matching to bin segments from different reads into the same unit
#####

def find_unit_startpos(readname, feature_calls, anchors, spacings):
	"""Return the contig unit number of the read's first unit, using the first anchor that matches exactly once."""
	offset = 0
	### loop through the anchors
	for anchor, spacing in zip(anchors, spacings):
		### increase the offset for later anchors
		offset += int(spacing)
		### count all occurences of the anchor
		n_matches = len(re.findall(anchor, feature_calls))

		### If one match, determine the number of units left of the anchor (num_L_units) in order
		### to the set the 'frame' of the read for binning
		if n_matches == 1:
			print(f"read {readname} matched with anchor at unit {offset}")
			left_calls = re.split(anchor, feature_calls)[0]
			start_unit = offset - (len(left_calls.split(",")) - 1)
			if start_unit < 1:
				raise Exception(f"\n\tread {readname} extends further than the first unit number at the 5'\n")
			return start_unit

		### Output warning if more than one match
		if n_matches > 1:
			print(f"\nanchor\n{anchor}\nwas present more than once in read {readname}\n")

	### Return none if no proper matches found
	return None

##### Loop through the features file and fastq for matching the anchors and binning read segments
fastq_seq_dict = defaultdict(list)
fastq_len_dict = defaultdict(list)
with open(features_file, 'r') as feature_file_obj, open(read_info_file, 'w') as read_info_file_obj:

	### open the features file
	feature_reader = csv.reader(feature_file_obj, delimiter='\t', quoting=csv.QUOTE_NONE)

	### loop through each DNA sequence and feature calls to find anchor matches
	for feature_line,fastq_entry in zip(feature_reader,pysam.FastxFile(oriented_fastq_file), strict=True):
		### Set initial variables
		cur_readname = feature_line[0]
		cur_featcalls = feature_line[1]

		### Check if the readname is usable
		if cur_readname != fastq_entry.name:
			raise Exception(f"\n\tThe .feat file readname {cur_readname} does not match the paired entry"
				f" in the oriented fastq {fastq_entry.name}.\n")
		if cur_readname in used_readnames:
			continue

		##### Find the position of the first unit in the read by searching for all anchors in it
		start_unit = find_unit_startpos(cur_readname, cur_featcalls, anchors, spacings)
		if start_unit is None:
			continue

		## Prepare stuff
		num_units = len(cur_featcalls.split(","))
		feature_starts = feature_line[2].split(",")
		feature_starts.append(len(fastq_entry.sequence))
		read_dir = '+'

		## Write read info
		feature_starts_string = ",".join([str(x) for x in feature_starts])
		read_info_file_obj.write(f"{cur_readname}\t{start_unit}\t{num_units}\t{read_dir}\t{feature_starts_string}\n")

		## split the read into segments binned into each unit number of the assembling contig
		for n in range(num_units):
			## set numbering and nt positions
			unit_num = n + start_unit
			start_site = int(feature_starts[n])
			end_site = int(feature_starts[n+1])

			## create segment fastq object
			current_seg_read = pysam.FastxRecord()
			current_seg_read.name = f"{fastq_entry.name}_segment{n}_unit{unit_num}"
			current_seg_read.sequence = fastq_entry.sequence[start_site:end_site]
			current_seg_read.quality = fastq_entry.quality[start_site:end_site]

			## add segment fastq into dictionary
			fastq_seq_dict[unit_num].append(current_seg_read)
			fastq_len_dict[unit_num].append(end_site - start_site)

		# ## Loop backwards (so as not to disturb downstream insertion sites)
		# ## through each units, adding tags
		# for n in range(num_units-1, 0, -1):
		#   start_tag_num = n-num_L_units+10+offset
		#   start_site = int(feature_starts[n])
		#   #print(f"\tinserting tag #{start_tag_num}, named {input_tags[start_tag_num].name}, into pos {start_site}")
		#   working_seq, working_qual = tagit(input_tags[start_tag_num], start_site, working_seq, working_qual)

		# ## Write the tagged read output
		# output_read = pysam.FastxRecord()
		# output_read.name = f"{fastq_entry.name}_tagged"
		# output_read.sequence = working_seq
		# output_read.quality = working_qual
		# tagged_reads.write(str(output_read) + "\n")


##### Make the output directory
path = snakemake.output[1]
os.makedirs(path, exist_ok=True)

if not fastq_seq_dict:
	raise Exception(f"No reads matched any anchor for contig {snakemake.wildcards.contig}")


#####
##### Filter the sets of read segments for each unit
#####

def length_filter(seg_reads, seg_lengths, unit_label):
	"""For units deeper than 4, keep reads within 0.67-1.5x the median length. If 2 or fewer survive,
	return all reads unfiltered."""
	if len(seg_reads) <= 4:
		# return all reads
		return seg_reads

	### Calculate the median length and keep reads that are close
	median_len = statistics.median(seg_lengths)
	kept = [r for r, L in zip(seg_reads, seg_lengths) if median_len * 0.67 < L < median_len * 1.5]

	if len(kept) > 2:
		# return filtered reads
		return kept
	else:
		# Note high variability and return all
		print(f"\nThe segment reads at unit number {unit_label} were frequently dissimilar to the median length"
			  f" of {median_len}.\nsegment lengths: {seg_lengths}")
		return seg_reads

def quality_filter(seg_reads, q_cutoff, max_depth):
	"""Keep reads with arithmetic mean Phred > q_cutoff and at most max_depth best scoring reads.
	If 4 or fewer pass, return all input reads."""
	passing = []
	### Loop reads and gather simple mean q for those that pass cutoff
	for read in seg_reads:
		mean_q = sum(read.quality.encode()) / len(read.quality) - 33
		if mean_q > q_cutoff:
			passing.append((mean_q, read))

	### Return filtered or unfiltered reads based on number that pass
	if len(passing) > max_depth:
		passing.sort(key=lambda x: x[0], reverse=True)
		# return max depth best reads
		return [read for _, read in passing[:max_depth]], len(passing)
	elif len(passing) > 4:
		# return all q filtered reads
		return [read for _, read in passing], len(passing)
	else:
		# do not filter reads if less than 5
		return seg_reads, len(passing)

##### Loop segments to filter reads and write them out
entered_goodcoverage = False
exited_goodcoverage = False
metadata_lol = []

first_unit, last_unit = min(fastq_seq_dict), max(fastq_seq_dict)
width = len(str(last_unit))

for unit_int in range(first_unit, last_unit + 1):
	unit_num = str(unit_int).zfill(width)
	unit_seg_reads = fastq_seq_dict.get(unit_int, [])
	unit_lengths = fastq_len_dict.get(unit_int, [])
	unit_depth = len(unit_seg_reads)

	### Once coverage of a unit is higher than 1, start filtering and writing multi-fasta files.
	### (Coverage of 1 or 0 is never written.)
	if unit_depth > 1 and not exited_goodcoverage:
		entered_goodcoverage = True
		
		### Filter by length
		temp_seg_reads = length_filter(unit_seg_reads, unit_lengths, unit_num)

		### Filter by quality
		final_seg_reads, n_qfilt = quality_filter(temp_seg_reads, q_cutoff, max_depth)

		### Write retained segment reads as multi-fasta
		with open(f"{path}/seg_{unit_num}_reads.fa", "w") as output:
			for final_seg_read in final_seg_reads:
				output.write(f">{final_seg_read.name}\n")
				output.write(f"{final_seg_read.sequence}\n")

		## save the metadata for this segment
		metadata_lol.append([unit_num, str(unit_depth), str(len(temp_seg_reads)), str(n_qfilt), str(len(final_seg_reads))])


	### If coverage drops to 1 after going above, set exited_goodcoverage. This should only occur at the 3'end of a contig.
	elif unit_depth <= 1 and entered_goodcoverage == True and not exited_goodcoverage:
		exited_goodcoverage = True

	### If coverage has dropped, but goes back up to higher than 1, error. This is often the result of incorrect inputs.
	elif unit_depth > 1 and exited_goodcoverage == True:
		raise Exception("\n\tCoverage for the units dropped below 2 read segments then went back up again,"
			" indicating a hole in this contig.\n")

##### Write metadata output
with open(metadata_file, 'w') as metadata_file_obj:
	metadata_file_obj.write("unit#\tinitial_depth\tpost_length_filter\tq_filtered\tused_depth\n")
	for line in metadata_lol:
		metadata_file_obj.write("\t".join(line) + "\n")
