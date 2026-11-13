#!/usr/bin/env python3
import sys
import csv
from collections import Counter

##### Suggests unique patterns based on expected coverage only.
##### Counts every window ('kmer') of 1 to max kmer length units in a .feat file and reports
##### those found at about 1X the coverage expected for a kmer of that length.
##### Output patterns are exact matches of every variant call in each unit of the kmer.
##### Typically, longer kmers that contain a shorter one already reported are left out.

if len(sys.argv) != 4:
	print("\tERROR: please input a .feat file, an estimated coverage, and a max kmer length. Writes to stdout.")
	sys.exit(2)

feat_file = sys.argv[1]
est_cov = float(sys.argv[2])
max_kmer_len = int(sys.argv[3])


##### read in the feat file (variant calls only), parsed into a list of lists of full units
feats = []
total_units = 0
with open(feat_file, 'r') as feat_obj:
	for tsv_line in csv.reader(feat_obj, delimiter = "\t"):
		# ignore "%" markers on calls
		units = tsv_line[1].replace("%","").split(",")

		# count the two partial end units as one unit's worth of sequence
		total_units += len(units) - 1

		# keep only full units by ignoring both likely partial units at ends
		feats.append(units[1:-1])


##### Count the occurrences of every kmer, and the number of windows for each kmer length
kmer_counts = Counter()
num_windows = Counter()
for units in feats:
	### loop through each kmer length
	for kmer_len in range(1, max_kmer_len + 1):
		## slide a window of that many units along the read
		for start in range(len(units) - kmer_len + 1):
			kmer = tuple(units[start:(start + kmer_len)])
			kmer_counts[kmer] += 1
			num_windows[kmer_len] += 1


##### Calculate the expected coverage of a unique kmer at each length
baseline_expected = total_units / est_cov
exp_cov = {}
for kmer_len in num_windows:
	exp_cov[kmer_len] = num_windows[kmer_len] / baseline_expected

def cov_ratio(kmer):
	"""Return a kmer's count relative to the expected coverage of a unique kmer of its length."""
	return kmer_counts[kmer] / exp_cov[len(kmer)]


##### Output suggested unique kmers, sorted by length
# print header
print("pattern\tkmer_len\tcount\texp_cov")
### Loop kmers, sorted by length
for kmer in sorted(kmer_counts, key=lambda kmer: (len(kmer), kmer)):

	### skip kmers that are not at about 1X of expected coverage
	if not .7 <= cov_ratio(kmer) <= 1.3:
		continue

	### keep only the shortest unique kmers
	# skip unless the kmers one unit shorter at either end are both repeated
	if len(kmer) > 1:
		if cov_ratio(kmer[:-1]) <= 1.3 or cov_ratio(kmer[1:]) <= 1.3:
			continue

	print(f"""{",".join(kmer)}\t{len(kmer)}\t{kmer_counts[kmer]}\t{exp_cov[len(kmer)]:.1f}""")