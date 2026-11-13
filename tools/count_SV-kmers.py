#!/usr/bin/env python3
import sys
from collections import defaultdict
import tagged_longread_mod as mod

### Takes two .feat files and counts the number of occurences of SV sequence 'kmers' using
### a user defined k and only kmers present in at least one of the files.
### Can be used on long read .feat files and assembly .feat files, and is intended to
### compare the two types.

if len(sys.argv)!=4:
	print("\n\tERROR: please input two .feat files and a k value for kmer length. Prints to stdout\n",file=sys.stderr)
	sys.exit(2)

feat_file1 = sys.argv[1]
feat_file2 = sys.argv[2]
k = int(sys.argv[3])

###read in the feat file
with open(feat_file1, 'r') as feat_obj1, open(feat_file2, 'r') as feat_obj2:
	feats_lol1 = mod.read_in_feat_aslol(feat_obj1)
	feats_lol2 = mod.read_in_feat_aslol(feat_obj2)


def count_kmers(in_lol, k):
	kmer_counts = defaultdict(int)
	for feat_list in in_lol:
		for i in range(len(feat_list) - k + 1):
			win_kmer = tuple(feat_list[i:i + k])
			kmer_counts[win_kmer] += 1
	return kmer_counts

###Count em
counts1 = count_kmers(feats_lol1, k)
counts2 = count_kmers(feats_lol2, k)


###Output matched lists
# make a set of all keys present in either
all_keys = counts1.keys() | counts2.keys()

# Print the counts for each kmer, sorted primarily by counts1, highest first
for kmer in sorted(all_keys, key=lambda km: (-counts1.get(km, 0), -counts2.get(km, 0), km)):
	count1 = counts1.get(kmer, 0)  # Get count from dict1 or 0 if not present
	count2 = counts2.get(kmer, 0)  # Get count from dict2 or 0 if not present
	print("\t".join(kmer) + f"\t{count1}\t{count2}")