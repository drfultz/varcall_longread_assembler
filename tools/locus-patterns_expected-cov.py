#!/usr/bin/env python3
import sys
import csv
import re
import math
import tagged_longread_mod as mod

##### Estimates how close to the expected coverage each pattern in a _pattern.tsv file is,
##### normalized based on the distribution of unit 'windows' of a given length.


if len(sys.argv) != 4:
	print("\tERROR: please input a .feat file, a patterns file, and an estimated coverage value. Writes to stdout.")
	sys.exit(2)
	
feat_file = sys.argv[1]
patterns_file = sys.argv[2]
est_cov = int(sys.argv[3])


### read the feat file into memory
with open(feat_file, 'r') as feat_obj:
	feat_fields = mod.read_in_feat(feat_obj)

### read search patterns into memory
with open (patterns_file) as patterns_obj:
	patterns = [pattern.strip() for pattern in patterns_obj.readline().split("\t")]

### Make distribution of unit length counts
distribution = {}
for features in feat_fields:
	### Count delimiter occurence, plus one for final unit
	count = features.count(",") + 1

	num_wins = 1
	for x in range(count,0,-1): 

		if x in distribution:
			distribution[x] += num_wins
		else:
			distribution[x] = num_wins
		num_wins += 1

norm_factor = distribution[1] / est_cov

### loop through patterns to get counts of pattern occurences
for pattern in patterns:
	###
	count = 0 
	for features in feat_fields:
		## Find pattern by re.split
		resplit_feature_entry = re.split(r'%s' % pattern, features)
		## Add 1 for each occurence of the pattern in the line (multiple occurences count more than once)
		count += len(resplit_feature_entry) - 1

	pat_len = pattern.count("+,") + len(re.findall(r'\w,', pattern)) + 1
	norm_FC = math.log2((count + 0.0000001)/ (distribution[pat_len] / norm_factor))

	print(f"{pattern}\t{count}\t{norm_FC}")
	
