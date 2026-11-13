#!/usr/bin/env python3
import sys
import re
import math
from collections import defaultdict
import tagged_longread_mod as mod

##### Estimates how close to the expected coverage each pattern in an _anchors.tsv file is,
##### normalized based on the distribution of unit 'windows' of a given length.


if len(sys.argv) != 4:
	print("\tERROR: please input a .feat file, an anchors file, and an estimated coverage value. Writes to stdout.")
	sys.exit(2)
	
feat_file = sys.argv[1]
patterns_file = sys.argv[2]
est_cov = int(sys.argv[3])


### read the feat file into memory
with open(feat_file, 'r') as feat_obj:
	feat_fields = [",".join(units) for units in mod.read_in_feat_aslol(feat_obj)]
	
### read search patterns into memory
with open (patterns_file) as patterns_obj:
	patterns = [pattern.strip() for pattern in patterns_obj.readline().split("\t")]

### Make distribution of unit length counts
distribution = defaultdict(int)
for features in feat_fields:
	### Count delimiter occurence, plus one for final unit
	feat_length = features.count(",") + 1

	### Loop the possible window sizes a read can contain
	for x in range(1, feat_length + 1):
		# Add a count equal to the total number of windows possible at that size
		distribution[x] += feat_length - x + 1

norm_factor = distribution[1] / est_cov

### loop through patterns to get counts of pattern occurences
for pattern in patterns:
	matches = 0 
	for features in feat_fields:
		## Count pattern occurences
		matches += sum(1 for _ in re.finditer(pattern, features))

	pat_len = re.sub(r"\[[^\]]*\]", "", pattern).count(",") + 1 # remove bracketed text then count remaining commas
	norm_FC = math.log2((matches + 0.0000001)/ (distribution[pat_len] / norm_factor))

	print(f"{pattern}\t{matches}\t{norm_FC}")
	