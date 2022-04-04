#!/usr/bin/env python3
import sys
import csv
import re
import itertools

##### Starting from a seed pattern, attempt to find additional contiguous unique patterns

if len(sys.argv)!=5:
	print("\tERROR: please input a .feat file, an estimated coverage, a seed pattern, and a max window size")
	sys.exit(2)

feat_file = sys.argv[1]
est_cov = int(sys.argv[2])
current_pattern = sys.argv[3]
max_SVkmer_len = int(sys.argv[4])

patterns_out = feat_file.split(".")[0] + "_patterns.tsv"

### read in the feat file
with open(feat_file, 'r') as feat_obj:
	feat_fields = []
	switch = False
	for line in csv.reader(feat_obj, delimiter = "\t"):
		# ignore "^" markers on calls
		feat_fields.append(line[1].replace("^",""))
		if switch == False:
			num_features_present = len(line[1].split(",")[0].split("&"))
			switch = True

##### Make distribution of unit length counts and normalize
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
			##maybe use this part to trim any leading or trailing unmapped calls

	### calculate spacing offset for whole contig
	extra_spacing = spacing - L_offset
	if extra_spacing < 0:
		extra_spacing = 0

	### now that max offsets found, go back and add spaces to make all sublists equal length
	aligned_features = []
	for entry in split_features:
		L_spaces = L_offset - entry[1]
		for x in range((L_spaces+extra_spacing)*num_feat):
			entry[0].insert(0, "-")

		R_spaces = R_offset - entry[2]
		for y in range(R_spaces*num_feat):
			entry[0].append("-")

		aligned_features.append(entry[0])

	return aligned_features

##### align patterns and find new ones until can't find unique pattern or consensus is ambiguous
contiguous = True
patterns_line = [current_pattern]
spacing_line = ["0"]
absolute_spacing = 40 #!!!!! need to fix this
while contiguous == True:
	### align current pattern
	current_table = align_varcall_table(current_pattern, feat_fields, num_features_present, absolute_spacing)
	tp_current_table = list(map(list, zip(*current_table)))

	### evaluate and make consensus
	## for column in original table (using transposed table)
	consensus_calls = []
	for col in tp_current_table:
		col_total = 0
		calls = {}
		## go through each field to tally calls present
		for field in col:
			## if not a -
			if field != "-":
				## add to total
				col_total += 1

				## count occurences of every call
				if field in calls:
					calls[field] += 1
				else:
					calls[field] = 1

		## make a consensus
		if col_total >= 4:
			highest_count = 0
			for call in calls:
				if calls[call] > highest_count:
					consensus_call = call
					highest_count = calls[call]
		else:
			consensus_call = "-"
		consensus_calls.append(consensus_call)

		## if sufficient depth, look for bad columns
		if col_total >= 10:
			## tally mismatches
			col_mismatches = 0
			for field in col:
				if field != "-" and field[-8:] != "Unmapped":
					if field != consensus_call:
						col_mismatches += 1
			## reject if any single column shows a divided pattern
			if col_mismatches/col_total > .3:
				contiguous = False
				print("bad column")
				print(col)

	### print alignment table and consensus regardless of quality
	print("\n" + current_pattern)
	for aligned_line in current_table:
		print("\t".join(aligned_line))
	#print("\t".join(consensus_calls))

	### stop if columns bad
	if contiguous == False:
		break


	##### look for next search pattern
	#####
	### convert consensus back to string for re searching
	cons_string = ""
	for i, cons_call in enumerate(consensus_calls, start=1):
		if cons_call != "-":	
			if i % num_features_present != 0:
				cons_string += cons_call + "&"
			else:
				cons_string += cons_call + ","
	#print(cons_string)

	## Find pattern by re.split
	resplit_cons = re.split(r'(%s)' % current_pattern, cons_string[:-1]) #cons string has trailing character

	## If pattern is found, parse it
	if len(resplit_cons) < 3:
		raise Exception("\npattern not found in consensus.\n")
	else:
		R_set = resplit_cons[2].split(",")
	#print("\t".join(R_set))

	### Greedy search for next pattern
	## loop sizes first
	for size in range(1,max_SVkmer_len+1):
		## loop through right most patterns, back to left
		for w in range(1,len(R_set[size:])+1):
			#print(f"{-w-size}\t{-w}")
			current_set = R_set[(-w-size):-w]

			## make string
			pre_re_string = r""
			pre_re_string += ",".join(current_set)

			## using lookahead with internal capture to not consume overlapping matches
			re_string = r"(?=(" + pre_re_string + r"))"
			#print(re_string)

			## count occurences
			count = 0
			for feat_field in feat_fields:
				count += len(re.findall(re_string, feat_field))

			## if expected coverage found, keep as next pattern
			current_exp_cov = distribution[len(current_set)] / norm_factor
			#print(f"{current_exp_cov}\t{count}")
			if 0.7*current_exp_cov <= count <= 1.3*current_exp_cov:
				#
				relative_spacing = resplit_cons[1].count(",") + len(R_set) - w - size
				absolute_spacing += relative_spacing
				print(f"\n{relative_spacing}\t{current_exp_cov}\t{count}")
				current_pattern = pre_re_string

				### generate outputs
				patterns_line.append(current_pattern)
				spacing_line.append(str(relative_spacing))
				break

		else: # no break at this size
			# Continue the outer loop
			continue
		# inner loop was broken, break the outer.
		break

	else: # no break in finding new pattern
		break # break out of building patterns

##### Output patterns file
with open(patterns_out, 'w') as patterns_obj:
	patterns_obj.write("\t".join(patterns_line) + "\n")
	patterns_obj.write("\t".join(spacing_line) + "\n")


# #####
# # intially set current as filtered set and real as empty
# current_combos = filt_combos
# real_combos = [[]]
# ### for each length up to max length
# for size in range(1,max_SVkmer_len+1):

# 	### set values
# 	working_combos = []
# 	exp_cov = distribution[size] / norm_factor
# 	low_cov = 0.7*exp_cov
# 	high_cov = 1.3*exp_cov

# 	print(f"\nsize {size}\t exp {exp_cov}")

# 	### Look for multi unit combo that is in expected range
# 	##loop through existing multi combos
# 	for seed_combo in current_combos:
# 		##loop through combos to add
# 		for add_combo in real_combos:
# 			## add combo to seed
# 			multi_combo = seed_combo + add_combo

# 			## make string
# 			pre_re_string = r""
# 			pre_re_string += "[^,]*,[^,]*".join(multi_combo)

# 			## using lookahead with internal capture to not consume overlapping matches, I think?
# 			re_string = r"(?=(" + pre_re_string + r"))"
# 			#print(re_string)

# 			## count occurences
# 			count = 0
# 			for feat_field in feat_fields:
# 				count += len(re.findall(re_string, feat_field))

# 			### if real
# 			if count > low_cov:
# 				## if within unique pattern range
# 				if count < high_cov:
# 					print("\n" + pre_re_string + "\t" + str(count))

# 					### run function
# 					aligned_table = align_varcall_table(pre_re_string,feat_fields)
# 					tp_aligned_table = list(map(list, zip(*aligned_table)))

# 					### evaluate and make consensus
# 					## for column in original table (using transposed table)
# 					consensus_calls = []
# 					for col in tp_aligned_table:
# 						col_total = 0
# 						calls = {}
# 						## go through each field to tally calls present
# 						for field in col:
# 							## if not a -
# 							if field != "-":
# 								## add to total
# 								col_total += 1

# 								## count occurences of every call
# 								if field in calls:
# 									calls[field] += 1
# 								else:
# 									calls[field] = 1

# 						## test col for matching
# 						if col_total >= 4:
# 							highest_count = 0
# 							for call in calls:
# 								if calls[call] > highest_count:
# 									consensus_call = call
# 									highest_count = calls[call]
# 						else:
# 							consensus_call = "-"
# 						consensus_calls.append(consensus_call)

# 						## tally mismatches
# 						mismatches = 0
# 						for field in col:
# 							if field != "-" and field[-8:] != "Unmapped":
# 								if field != consensus_call:
# 									mismatches += 1


# 					# if table is properly matching
# 					#if threshold:
# 					# for aligned_line in aligned_table:
# 					# 	print("\t".join(aligned_line))

# 					print("\t".join(consensus_calls))

# 				## else include in real combos for further use
# 				else:
# 					working_combos.append(multi_combo)
# 					#print("added for next round")

# 	### write over initial combo set with smaller passing set
# 	print(working_combos)
# 	if size == 1:
# 		real_combos = working_combos
# 	current_combos = working_combos
