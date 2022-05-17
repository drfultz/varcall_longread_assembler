#!/usr/bin/env python3
import sys
import csv
import re
import itertools

#####

if len(sys.argv)!=5:
	print("\tERROR: please input a .feat file, a file with target SVs, an estimated coverage, and a max window size")
	sys.exit(2)

feat_file = sys.argv[1]
target_SVs_file = sys.argv[2]
est_cov = int(sys.argv[3])
max_SVkmer_len = int(sys.argv[4])

##### read in files
##read in the target SVs that will be used to generate SV kmers
with open(target_SVs_file, 'r') as target_SVs_obj:
	target_SVs = []
	for line in csv.reader(target_SVs_obj, delimiter = "\t"):
		variants = line[1].split(',')
		l = [""]
		for variant in variants:
			l.append(f"{line[0]}:{variant}")
		target_SVs.append(l)

##read in the feat file
with open(feat_file, 'r') as feat_obj:
	feat_fields = []
	for line in csv.reader(feat_obj, delimiter = "\t"):
		# ignore "^" markers on calls
		feat_fields.append(line[1].replace("^",""))

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

##### Generate all single unit SV combinations and filter to those specifying 2-3 features
### Generate all possible combos
all_combos = list(itertools.product(*target_SVs))

### Filter the SV combos to those with only 2 specified
filt_combos = []
# For each combo
for combo in all_combos:
	spec_count = 0
	filt_comboL = []
	# loop through the positions and count and keep only those that aren't empty strings
	for pos in combo:
		if pos != "":
			spec_count += 1
			filt_comboL.append(pos)
	# if between 2 and 2 SVs specified, keep that as part of filt combos
	if 1 < spec_count < 3:
		string_of2 = "[^,]*".join(filt_comboL)
		#only keep patterns specifying IGS as well ### needs to be made customizable
		if "B:" in string_of2:
			filt_combos.append([string_of2])
			#print(string_of2, file=sys.stderr)


def align_varcall_table(search_pattern, varcalls):
	### make varcall alignment
	split_features = []
	L_offset = 0
	R_offset = 0
	switch = False

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

			## count the number of calls per unit
			if switch is False:
				num_features_present = len(line.split(",")[0].split("&")) 
				switch = True

			## split the features entry into units and calls and put into lol
			split_features.append([line.replace("&",",").split(","), num_L_units, num_R_units])

	### now that max offsets found, go back and add spaces to make all sublists equal length
	aligned_features = []
	for entry in split_features:
		L_spaces = L_offset - entry[1]
		for x in range(L_spaces*num_features_present):
			entry[0].insert(0, "-")

		R_spaces = R_offset - entry[2]
		for y in range(R_spaces*num_features_present):
			entry[0].append("-")

		aligned_features.append(entry[0])

	return aligned_features


#####
# intially set current as filtered set and real as empty
current_combos = filt_combos
real_combos = [[]]
### for each length up to max length
for size in range(1,max_SVkmer_len+1):

	### set values
	working_combos = []
	exp_cov = distribution[size] / norm_factor
	low_cov = 0.7*exp_cov
	high_cov = 1.3*exp_cov

	print(f"\nsize {size}\t exp {exp_cov}")

	### Look for multi unit combo that is in expected range
	##loop through existing multi combos
	for seed_combo in current_combos:
		##loop through combos to add
		for add_combo in real_combos:
			## add combo to seed
			multi_combo = seed_combo + add_combo

			## make string
			pre_re_string = r""
			pre_re_string += "[^,]*,[^,]*".join(multi_combo)

			## using lookahead with internal capture to not consume overlapping matches, I think?
			re_string = r"(?=(" + pre_re_string + r"))"
			#print(re_string)

			## count occurences
			count = 0
			for feat_field in feat_fields:
				count += len(re.findall(re_string, feat_field))

			### if real
			if count > low_cov:
				## if within unique pattern range
				if count < high_cov:
					print("\n" + pre_re_string + "\t" + str(count))

					### run function
					aligned_table = align_varcall_table(pre_re_string,feat_fields)
					tp_aligned_table = list(map(list, zip(*aligned_table)))

					### evaluate and make consensus
					## for column in original table (using transposed table)
					consensus_calls = []
					for col in tp_aligned_table:
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

						## test col for matching
						if col_total >= 4:
							highest_count = 0
							for call in calls:
								if calls[call] > highest_count:
									consensus_call = call
									highest_count = calls[call]
						else:
							consensus_call = "-"
						consensus_calls.append(consensus_call)

						## tally mismatches
						mismatches = 0
						for field in col:
							if field != "-" and field[-8:] != "Unmapped":
								if field != consensus_call:
									mismatches += 1


					# if table is properly matching
					#if threshold:
					# for aligned_line in aligned_table:
					# 	print("\t".join(aligned_line))

					print("\t".join(consensus_calls))

				## else include in real combos for further use
				else:
					working_combos.append(multi_combo)
					#print("added for next round")

	### write over initial combo set with smaller passing set
	print(working_combos)
	if size == 1:
		real_combos = working_combos
	current_combos = working_combos
