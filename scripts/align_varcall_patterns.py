#!/usr/bin/env python3
import sys
import re



if len(sys.argv)!=3:
	print("\tERROR: please input a .feat file and a regular expression search string")
	sys.exit(2)
	
features_file = sys.argv[1]
search_pattern = sys.argv[2]
output_file = "%s_aligned.tsv" % (features_file.split('.', 1)[0])

with open(features_file,'r') as features_file_object:
	features_parsed = []
	offset = 0
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
			#print(num_L_units)

			## Keep track of the largest offset
			if num_L_units > offset:
				offset = num_L_units

			## split the features entry into units
			split_feature_entry = split_line[1].split(",")
			
			## split the units into individual variant calls
			features_lol = []
			for unit in split_feature_entry:
				#print(unit.split("&"))
				features_lol.append(unit.split("&"))

				if switch is False:
					num_called_features = len(unit.split("&")) 
					switch = True
			#print(features_lol)

			features_parsed.append([readname, features_lol, num_L_units])

#with open(output_file, 'w') as output:
for entry in features_parsed:
	spaces = offset - entry[2]
	spaced_entries = spaces*num_called_features*"NA\t"

	features_string = ""
	for x in entry[1]:
		for y in x:
			features_string = features_string + "\t" + y

	print(f"{entry[0]}\t{spaced_entries}{features_string}")
	