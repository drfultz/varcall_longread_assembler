#!/usr/bin/env python3
import sys
import re

if len(sys.argv)!=4:
	print("\tERROR: please input a .feat file, a regular expression search string, and how many features for scaffold reads")
	sys.exit(2)
	
features_file = sys.argv[1]
search_pattern = sys.argv[2]
num_desired_features = int(sys.argv[3])

##### Parse the input features file and measure the distances on either side of the search pattern
##### for every feature line
with open(features_file,'r') as features_file_object:
	features_parsed = []
	L_offset = 0
	R_offset = 0
	switch = False

	for line in features_file_object:
		split_line = line.split("\t")
		readname = split_line[0]

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

			features_parsed.append([readname, features_lol, num_L_units, num_R_units])

print(f"coverage = {len(features_parsed)}",file=sys.stderr)

##### Make preperations for finding and processing reads that extend the furthest 5' from the pattern
## Define a function that outputs one unit of variant calls, limited to a certain number of features,
## per entry in a list for use later in outputting one unit per line.
def trim_entry():
	trimmed_entry = [entry[0],0]
	Lempty_features = ["-\t"*num_desired_features]*L_spaces
	select_features = []
	Rempty_features = ["-\t"*num_desired_features]*R_spaces
	for unit in entry[1]:
		select_features.append("\t".join(unit[0:num_desired_features]) + "\t")
	trimmed_entry[1] = Lempty_features + select_features + Rempty_features
	return trimmed_entry

longest_length = 10000
longest_found = False


##output search pattern as first line
print(search_pattern)

##### Output appropriatly spaced feature calls and save the longest feature lines for next step
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
	if L_spaces == 0 and longest_found == False:
		trimmed_longest = trim_entry()
		longest_found = True

	elif L_spaces < longest_length:
		longest_length = L_spaces
		trimmed_longest2 = trim_entry()

##### Output the  unit by unit transposed longest reads for use in manual scaffolding
print("\n")
print("\t" + trimmed_longest[0] + "\t"*num_desired_features + trimmed_longest2[0])
for one,two in zip(trimmed_longest[1],trimmed_longest2[1]):
	print("\t" + one + two)
