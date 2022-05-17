#!/usr/bin/env python3
import sys
import csv
import re

##### Using a locus patterns file (_patterns.tsv), align all varcall reads at each pattern
##### along the entire contig.

if len(sys.argv)!=3:
	print("\tERROR: please input a .feat file and a locus patterns file. Prints to STDOUT.")
	sys.exit(2)

feat_file = sys.argv[1]
patterns_file = sys.argv[2]

##### read in the feat file
with open(feat_file, 'r') as feat_obj:
	feat_fields = []
	switch = False
	for tsv_line in csv.reader(feat_obj, delimiter = "\t"):

		split_units = tsv_line[1].split(",")
		trimmed_units = []
		### loop through each unit
		for u, unit in enumerate(split_units):
			## if first or last unit
			if u == 0 or u == (len(split_units)-1):
				# keep unit if any call is not Unmapped
				for subunit in unit.split("&"):
					if subunit[-8:] != "Unmapped":
						trimmed_units.append(unit)
						break
				# else: # no break
				# 	print(f"all unmapped {u}")
			## keep all internal units
			else:
				trimmed_units.append(unit)

		# join trimmed units back together, ignore "^" markers on calls, and put into list
		feat_fields.append(",".join(trimmed_units).replace("^",""))

		# count number of features for first line
		if switch == False:
			num_features_present = len(tsv_line[1].split(",")[0].split("&"))
			switch = True

##### read in patterns file
with open(patterns_file, 'r') as patterns_obj:
	patterns_line = patterns_obj.readline().strip().split("\t")
	spacings_line = patterns_obj.readline().strip().split("\t")

print(patterns_line)
print(spacings_line)

def align_varcall_table(search_pattern, varcalls, num_feat, spacing):
	### make varcall alignment
	split_features = []
	L_offset = 0
	R_offset = 0

	### gather any lines that match pattern
	for line in varcalls:
		## Find pattern by re.split
		resplit_feature_entry = re.split(r'%s' % search_pattern, line)

		## If pattern is found
		if len(resplit_feature_entry) > 1:

			## Find the number of units left of the pattern
			L_units = resplit_feature_entry[0]
			num_L_units = len(L_units.split(",")) - 1
			## Keep track of the largest L_offset
			if num_L_units > L_offset:
				L_offset = num_L_units

			## Find the number of units right of the pattern
			R_units = resplit_feature_entry[1]
			num_R_units = len(R_units.split(","))
			## Keep track of the largest R_offset
			if num_R_units > R_offset:
				R_offset = num_R_units

			## split the features entry into units and calls and put into lol
			split_features.append([line.replace("&",",").split(","), num_L_units, num_R_units])

	### calculate spacing offset for whole contig
	extra_spacing = spacing - L_offset
	if extra_spacing < 0:
		extra_spacing = 0

	### now that max offsets found, go back and add spaces to make all sublists equal length
	aligned_features = []
	for entry in split_features:
		L_spaces = L_offset - entry[1]
		for x in range((L_spaces+extra_spacing)*num_feat):
			entry[0].insert(0, "")

		R_spaces = R_offset - entry[2]
		for y in range(R_spaces*num_feat):
			entry[0].append("")

		aligned_features.append(entry[0])

	return aligned_features


##### Output alignment file
abs_spacing = 0
for pattern,rel_spacing in zip(patterns_line,spacings_line):
	abs_spacing += int(rel_spacing)
	print("\t"*abs_spacing*num_features_present + pattern)
	for out_line in align_varcall_table(pattern, feat_fields, num_features_present, abs_spacing):
		print("\t".join(out_line))
