#!/usr/bin/env python3
import sys
import csv
import re
import itertools

###Takes a .feat file and counts the number of occurences of SV sequence 'kmers' using
###only SVs defined in the second input, at a single length of SV kmer defined by third input.
###Can be used on long read .feat files and assembly .feat files.

if len(sys.argv)!=4:
	print("\tERROR: please input a .feat file, a file with target SVs, and the length of SV kmers")
	sys.exit(2)

feat_file = sys.argv[1]
target_SVs_file = sys.argv[2]
max_SVkmer_len = int(sys.argv[3])

###read in the target SVs that will be used to generate SV kmers
with open(target_SVs_file, 'r') as target_SVs_obj:
	target_SVs = []
	for line in csv.reader(target_SVs_obj, delimiter = "\t"):
		variants = line[1].split(',')
		l = []
		for variant in variants:
			l.append(f"{line[0]}:{variant}")
		target_SVs.append(l)

###Generate all combos
all_combos = list(itertools.product(*target_SVs))

###read in the feat file
with open(feat_file, 'r') as feat_obj:
	feat_fields = []
	for line in csv.reader(feat_obj, delimiter = "\t"):
		feat_fields.append(line[1])

###find the SV kmers
for int_combo in itertools.product(*[all_combos for i in range(max_SVkmer_len)]):
	pre_re_string = r""
	unit_pieces = []
	for unit in int_combo:
		unit_piece = "[^,]*".join(unit)
		unit_pieces.append(unit_piece)

	pre_re_string = "[^,]*,[^,]*".join(unit_pieces)

	re_string = r"(?=(" + pre_re_string + r"))"
	#print(re_string)

	count = 0
	for feat_field in feat_fields:
		#print(len(re.findall(re_string, feat_field)))
		count += len(re.findall(re_string, feat_field))

	print(str(int_combo) + "\t" + str(count))
