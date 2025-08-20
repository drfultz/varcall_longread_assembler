#!/usr/bin/env python3
import sys
import csv
from collections import defaultdict
import tagged_longread_mod as mod

###Takes a .feat file and counts the number of occurences of SV sequence '3mers' using
###only 3mers that exist in .
###Can be used on long read .feat files and assembly .feat files.

if len(sys.argv)!=3:
	print("\n\tERROR: please input two .feat files. Prints to stdout\n",file=sys. stderr)
	sys.exit(2)

feat_file1 = sys.argv[1]
feat_file2 = sys.argv[2]

###read in the feat file
with open(feat_file1, 'r') as feat_obj1, open(feat_file2, 'r') as feat_obj2:
	feats_lol1 = mod.read_in_feat_aslol(feat_obj1)
	feats_lol2 = mod.read_in_feat_aslol(feat_obj2)

def count_3mers(in_lol):
	three_mer_counts = defaultdict(int)
	for feat_list in in_lol:
		for i in range(len(feat_list) - 2):
			# Create a tuple of the current 3-mer
			three_mer = (feat_list[i], feat_list[i + 1], feat_list[i + 2])
			# Increment the count for this 3-mer
			three_mer_counts[three_mer] += 1

	return three_mer_counts


###Count em
counts1 = count_3mers(feats_lol1)
counts2 = count_3mers(feats_lol2)

###Output matched lists
# make a set of all keys present in both
all_keys = set(counts1.keys()).union(set(counts2.keys()))

# Print the counts for each 3-mer
for three_mer in all_keys:
	count1 = counts1.get(three_mer, 0)  # Get count from dict1 or 0 if not present
	count2 = counts2.get(three_mer, 0)  # Get count from dict2 or 0 if not present
	print(f"{three_mer}\t{count1}\t{count2}")

