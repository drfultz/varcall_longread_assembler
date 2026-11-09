#!/usr/bin/env python3
import sys
import os
import array
import pysam
from collections import defaultdict
import re

#####
#####

#in_bam = sys.argv[1]
#in_positions = sys.argv[2]

in_bam = snakemake.input[0]
in_positions = snakemake.input[1]
in_metadata = snakemake.input[2]

##### Parse the splitting positions
with open(in_positions, "r") as in_positions_obj:
	read_info_dict = {}
	for in_line in in_positions_obj:
		fields = in_line.strip().split("\t")
		read_info_dict[fields[0]] = [int(fields[1]),int(fields[2]),fields[3],fields[4].split(',')]

##### Read units from read_splitter
with open(in_metadata) as metadata_obj:
	next(metadata_obj) #skip header
	# Gather the 0 padded unit numbers as int keys but store the padded values
	unit_names_dict = {int(line.split("\t")[0]): line.split("\t")[0] for line in metadata_obj}

##### Make the output directory
path = snakemake.output[0]
os.makedirs(path, exist_ok=True)



def parse_mm_modes(mm_tag):
	"""Map (base, strand, mod_code) -> mode suffix ('', '.', or '?') for every section of an MM tag."""
	modes = {}
	for section in mm_tag.rstrip(";").split(";"):
		header = section.split(",", 1)[0]
		base, strand, codes, mode = re.fullmatch(r"([ACGTUN])([+-])([a-z]+|\d+)([.?]?)", header).groups()
		code_list = [int(codes)] if codes.isdigit() else list(codes)
		for code in code_list:
			modes[(base, 0 if strand == "+" else 1, code)] = mode
	return modes

def encode_mm_ml(new_mods, modes, sequence):
	"""Encode segment-local (0-based) modification calls as MM/ML values,
	keeping every section of the original tag and its mode."""
	mm_parts = []
	ml_values = []
	for (base, strand, mod_code), mode in modes.items():
		pos_prob_list = sorted(new_mods.get((base, strand, mod_code), []))
		mod_pos_set = {pos for pos, _ in pos_prob_list}
		# pysam reports positions of the base as written in MM, for either strand
		deltas, skip = [], 0
		for i, b in enumerate(sequence):
			if b != base:
				continue
			if i in mod_pos_set:
				deltas.append(skip)
				skip = 0
			else:
				skip += 1
		if len(deltas) != len(pos_prob_list):
			raise ValueError(f"{len(pos_prob_list)} {base}{mod_code} calls but only {len(deltas)} matched a {base}")
		strand_char = "+" if strand == 0 else "-"
		mm_parts.append(",".join([f"{base}{strand_char}{mod_code}{mode}"] + [str(d) for d in deltas]))
		ml_values.extend(prob for _, prob in pos_prob_list)
	mm_tag = ";".join(mm_parts) + ";" if mm_parts else ""
	return mm_tag, array.array("B", ml_values)


##### Iterate the reads for splitting
with pysam.AlignmentFile(in_bam, "rb", check_sq=False) as in_bam_obj:

	out_reads_dict = defaultdict(list)
	num_split_success = 0
	num_wrong_len = 0
	
	for read in in_bam_obj:
		##### Check each read first
		### Safeguard against mapped reads in bam
		if not read.is_unmapped:
			raise ValueError(f"{in_bam} contains aligned records; an unaligned BAM is required")

		### Find the split positions
		if read.query_name in read_info_dict:
			print(f"splitting bam read {read.query_name}")

			read_info_line = read_info_dict[read.query_name]
			read_positions = [int(x) for x in read_info_line[3]]

			### Check that read length from bam matches that of the fastq - they must be the same for splitting to work
			if read.query_length != read_positions[-1]:
				num_wrong_len += 1
				print(f"skipping {read.query_name} - read length mismatch between BAM and fastq", file=sys.stderr)
				continue

			num_split_success += 1
			mods = read.modified_bases
			modes = parse_mm_modes(read.get_tag("MM")) if read.has_tag("MM") else {}

			# Grab read orientation and length
			read_dir = read_info_line[2]
			if read_dir not in ('+','-'):
				raise ValueError(f"\nField 4 of the read info should be either + or - , "
					f"but was {read_dir}\n")
			read_len = read_positions[-1]

			for n in range(read_info_line[1]):
				unit_int = read_info_line[0] + n
				
				# If the unit is in the metadata from read_splitter, match it's 0 padded value
				if unit_int in unit_names_dict:
					current_unit = unit_names_dict[unit_int]
				else: # else skip it
					continue

				# flip coordinates if read was RC during .feat calling
				if read_dir == '+':
					start_site = read_positions[n]
					end_site = read_positions[n+1]
				else:
					start_site = read_len - read_positions[n+1]
					end_site = read_len - read_positions[n]

				current_seg_read = pysam.AlignedSegment()
				current_seg_read.query_name = f"{read.query_name}_segment{n}_unit{current_unit}"
				current_seg_read.query_sequence = read.query_sequence[start_site:end_site]
				current_seg_read.query_qualities = read.query_qualities[start_site:end_site]
				current_seg_read.flag = 4

				if modes:
					new_mod_field = {k: [(p - start_site, q) for p, q in v if start_site <= p < end_site]
									 for k, v in mods.items()}
					MM_tag, ML_tag = encode_mm_ml(new_mod_field, modes, current_seg_read.query_sequence)
					current_seg_read.set_tag("MM", MM_tag)
					current_seg_read.set_tag("ML", ML_tag)

				out_reads_dict[current_unit].append(current_seg_read)

	### Sanity check for having bam reads successfully pass through splitting
	num_unused_splitfq = len(read_info_dict) - num_split_success - num_wrong_len
	print(f"\n{num_split_success} bam reads split successfully\n"
		f"{num_wrong_len} bam reads mismatched their length from the fastq\n"
		f"{num_unused_splitfq} fastq reads (from read_info) did not have matching bam reads", file=sys.stderr)
	if num_split_success == 0:
		raise RuntimeError("No reads successfully split")
	elif num_split_success < 0.11*num_wrong_len:
		print("\n!!!!! More than 90 percent of matched bam reads did not match their length from the fastq."
			" Check that the uBAM reads are from the same basecalling as the fastq !!!!!\n", file=sys.stderr)


	##### Output one unaligned bam per unit from the read_splitter - produces empty bam if no bam segments were found
	for unit in unit_names_dict.values():
		with pysam.AlignmentFile(f"{path}/seg_{unit}_bamreads.bam", "wb", header=in_bam_obj.header) as out_bam_obj:
			for unit_read in out_reads_dict.get(unit, []):
				out_bam_obj.write(unit_read)
