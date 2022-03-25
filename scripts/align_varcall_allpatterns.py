#!/usr/bin/env python3
import sys
import re

if len(sys.argv)!=4:
	print("\tERROR: please input a .feat file, a contig pattern file, and how many features for scaffold reads")
	sys.exit(2)

features_file = sys.argv[1]
pattern_file = sys.argv[2]
num_desired_features = int(sys.argv[3])

### read search patterns into memory
with open (pattern_file) as patterns_obj:
	patterns = [pattern.strip() for pattern in patterns_obj.readline().split("\t")]

### read features file into memory
with open(features_file,'r') as features_obj:
	features_lines = []
	for line in features_obj:
		features_lines.append(line.strip().split("\t"))

### Define a function that outputs one unit of variant calls, limited to a certain number of features,
### per entry in a list for use later in outputting one unit per line.
def trim_entry():
	trimmed_entry = [entry[0],0]
	Lempty_features = ["-\t"*num_desired_features]*L_spaces
	select_features = []
	Rempty_features = ["-\t"*num_desired_features]*R_spaces
	for unit in entry[1]:
		select_features.append("\t".join(unit[0:num_desired_features]) + "\t")
	trimmed_entry[1] = Lempty_features + select_features + Rempty_features
	return trimmed_entry

### loop through search patterns and output a block of aligned reads for each one
for search_pattern in patterns:
	### reset values
	features_parsed = []
	L_offset = 0
	R_offset = 0
	switch = False

	### loop through the features lines to find the given pattern and place the values in to features_parsed
	for split_line in features_lines:
		## Find pattern by re.split
		resplit_feature_entry = re.split(r'%s' % search_pattern, split_line[1])
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

			## split the features entry into units
			split_feature_entry = split_line[1].split(",")

			## split the units into individual variant calls
			features_lol = []
			for unit in split_feature_entry:
				features_lol.append(unit.split("&"))

				if switch is False:
					num_features_present = len(unit.split("&"))
					switch = True
			features_parsed.append([split_line[0], features_lol, num_L_units, num_R_units])

	### output search pattern and coverage as first line
	print(f"{search_pattern}\t\t\t{len(features_parsed)}")

	### Make preparations for finding and processing reads that extend the furthest 5' from the pattern
	longest_left_length = 10000
	longest_right_length = 10000
	longest_left_found = False
	longest_right_found = False

	trimmed_longest1 = ["Pattern not found","-"]
	trimmed_longest2 = ["Pattern not found","-"]
	trimmed_longest3 = ["Pattern not found","-"]
	trimmed_longest4 = ["Pattern not found","-"]
	
	### Output appropriatly spaced feature calls and save the longest feature lines for next step
	for entry in features_parsed:
		## Add blanks to both sides of feature lines to length of longest lines
		L_spaces = L_offset - entry[2]
		L_spaced_entries = L_spaces*num_features_present*"-\t"

		R_spaces = R_offset - entry[3]
		R_spaced_entries = R_spaces*num_features_present*"-\t"

		features_string = ""
		for x in entry[1]:
			for y in x:
				features_string = features_string + y + "\t"

		print(f"{entry[0]}\t{L_spaced_entries}\t{features_string}{R_spaced_entries}")

		## Find two reads that stretch the most 5' to pattern
		#if number of L_spaces is 0, read is a longest read but only set as first entry if there isn't one already
		if L_spaces == 0 and longest_left_found == False:
			trimmed_longest1 = trim_entry()
			longest_left_found = True

		elif L_spaces < longest_left_length:
			longest_left_length = L_spaces
			trimmed_longest2 = trim_entry()

		if R_spaces == 0 and longest_right_found == False:
			trimmed_longest3 = trim_entry()
			longest_right_found = True
	
		elif R_spaces < longest_right_length:
			longest_right_length = R_spaces
			trimmed_longest4 = trim_entry()
			

	### Output the unit by unit transposed longest reads for use in manual scaffolding
	print("\n")
	print("\t" + trimmed_longest1[0] + "\t"*num_desired_features + trimmed_longest2[0] + "\t"*num_desired_features + trimmed_longest3[0] + "\t"*num_desired_features + trimmed_longest4[0])
	for one,two,three,four in zip(trimmed_longest1[1],trimmed_longest2[1],trimmed_longest3[1],trimmed_longest4[1]):
		print("\t" + one + two + three + four)
	print("\n")
