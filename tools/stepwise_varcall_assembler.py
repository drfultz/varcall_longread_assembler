#!/usr/bin/env python3
import sys
import re
import os
import itertools
import hashlib
import tagged_longread_mod as mod

##### Starting from a seed anchor pattern, attempt to find additional contiguous anchor patterns
##### Outputs a suggested _anchors.tsv and an SV align with a six hexadecimal hash as ID for that seed.

if len(sys.argv)!=6:
	print("\tERROR: please input a .feat file, an estimated coverage, a seed pattern, a max window size,"
			" and whether to look for more patterns (default to minimal patterns,"
			" unless M as last positional argument).", file=sys.stderr)
	sys.exit(2)

feat_file = sys.argv[1]
est_cov = int(sys.argv[2])
start_pattern = sys.argv[3]
max_SVkmer_len = int(sys.argv[4])
max_search = sys.argv[5]

seed_tag = hashlib.sha1(start_pattern.encode()).hexdigest()[:6]
out_prefix = f"{os.path.splitext(feat_file)[0]}_{seed_tag}"
patterns_out = f"{out_prefix}_draft_anchors.tsv"
alignment_out = f"{out_prefix}_draft_SValign.tsv"
print(f"running for seed {start_pattern} -> {seed_tag} seed ID")


##### read in the feat file
with open(feat_file, 'r') as feat_obj:
	feat_fields = [",".join(units) for units in mod.read_in_feat_aslol(feat_obj)]
	
# count number of features for first line
num_features_present = len(feat_fields[0].split(",")[0].split("&"))


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


def make_varcall_consensus(search_pattern1, varcalls1, num_feat1):
	### align current pattern
	current_table = [["-"] * start_col + calls
		for start_col, calls in mod.align_varcall_table(search_pattern1, varcalls1, num_feat1, 0)]
	tp_current_table = list(itertools.zip_longest(*current_table, fillvalue="-"))

	### evaluate and make consensus
	## for column in original table (using transposed table)
	consensus_calls = []
	good_columns = True
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

		## make a consensus if 4 or more in column
		if col_total >= 4:
			highest_count = 0
			for call in calls:
				if calls[call] > highest_count:
					consensus_call = call
					highest_count = calls[call]
			## don't keep low consensus
			if highest_count < 0.5*col_total:
				consensus_call = "--"
			## don't keep Unmapped consensus for low cov
			if col_total <= 6 and consensus_call[-8:] == "Unmapped":
				consensus_call = "--"

		## else, empty consensus value
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
				good_columns = False
				print("\nbad column", file=sys.stderr)
				print(col, file=sys.stderr)

	### convert consensus back to string for re searching
	cons_string = ""
	for i, cons_call in enumerate(consensus_calls, start=1):
		if cons_call != "-":	
			if i % num_features_present != 0:
				cons_string += cons_call + "&"
			else:
				cons_string += cons_call + ","

	return cons_string, good_columns

def test_cov(window, feat_fields, distribution, norm_factor, num_feat):
	"""Return (pattern, new consensus string, expected cov, count) if window has expected coverage and a clean alignment, else None."""
	pre_re_string = ",".join(window)

	## using lookahead with internal capture to not consume overlapping matches
	re_string = r"(?=(" + pre_re_string + r"))"
	count = sum(len(re.findall(re_string, feat_field)) for feat_field in feat_fields)

	## if expected coverage found and alignment is good
	exp_cov = distribution[len(window)] / norm_factor
	if 0.7*exp_cov <= count <= 1.3*exp_cov:
		new_cons_string, good = make_varcall_consensus(pre_re_string, feat_fields, num_feat)
		if good:
			return pre_re_string, new_cons_string, exp_cov, count

	## Either threshold not met
	return None

##### align patterns and find new ones until can't find unique pattern or consensus is ambiguous
### align and test start pattern
start_cons_string, good_start = make_varcall_consensus(start_pattern, feat_fields, num_features_present)

#if good_start == False:
#	### write the alignment for troubleshooting
#	with open(alignment_out, 'w') as alignment_obj:
#		for start_col, calls in mod.align_varcall_table(pattern, feat_fields, num_features_present, abs_spacing):
#				alignment_obj.write("\t" * start_col + "\t".join(calls) + "\n")
#
#	### end script
#	raise Exception("\nInitial pattern contains mismatching variant calls.\n")

### setup
patterns_line = [start_pattern]
spacing_line = ["50"]

### go forward and reverse
for F in [True,False]:
	## initial setup for each direction
	contiguous = True
	current_cons_string = start_cons_string
	current_pattern = start_pattern
	print(f"search_pattern\trelative_spacing\texp_cov\tactual_cov")

	##### main search loop
	### loop until a break or not contiguous
	while contiguous == True:

		## Find pattern by re.split
		resplit_cons = re.split(r'(%s)' % current_pattern, current_cons_string[:-1]) #cons string has trailing character

		## check that pattern was found
		if len(resplit_cons) < 3:
			raise Exception("\npattern not found in consensus.\n")

		## setup for search based on direction
		## take the appropriate half of resplit string, split into list at "," and replace double dashes
		if F == True:
			half_set = [s_unit.replace("--","[^,]+") for s_unit in resplit_cons[2].split(",")]
		else:
			half_set = [s_unit.replace("--","[^,]+") for s_unit in resplit_cons[0].split(",")]

		## loop sizes first
		for size in range(1,max_SVkmer_len+1):
			## loop through number of windows available for this size (closest first if max searching)
			window_order = range(len(half_set[size:]))
			if max_search == "M":
				window_order = reversed(window_order)
			for w in window_order:

				## define window to examine based on direction
				if F == True:
					print(f"{-w-size-1}\t{-w-1}", file=sys.stderr)
					current_win = half_set[(-w-size-1):(-w-1)]
				else:
					print(f"{w}\t{w+size}", file=sys.stderr)
					current_win = half_set[w:(w+size)]

				result = test_cov(current_win, feat_fields, distribution, norm_factor, num_features_present)
				if result:
					current_pattern, current_cons_string, exp_cov, count = result
					if F:
						relative_spacing = resplit_cons[1].count(",") + len(half_set) - w - 1 - size
						patterns_line.append(current_pattern)
						spacing_line.append(str(relative_spacing))
					else:
						relative_spacing = len(half_set) - w - 1
						patterns_line.insert(0, current_pattern)
						spacing_line[0] = str(relative_spacing)
						spacing_line.insert(0, "50")

					print(f"{current_pattern}\t{relative_spacing}\t{exp_cov}\t{count}")

					break

			else: ## no break at this size
				## Continue the outer loop
				continue
			## inner loop was broken, break the outer.
			break

		else: ## no break in finding new pattern, end search in this direction
			contiguous = False #so that while loop stops



##### Output patterns file
print(spacing_line)
with open(patterns_out, 'w') as patterns_obj:
	patterns_obj.write("\t".join(patterns_line) + "\n")
	patterns_obj.write("\t".join(spacing_line) + "\n")

##### Output alignment file
with open(alignment_out, 'w') as alignment_obj:
	abs_spacing = 0
	for pattern,rel_spacing in zip(patterns_line,spacing_line):
		abs_spacing += int(rel_spacing)
		alignment_obj.write("\t"*abs_spacing*num_features_present + pattern + "\n")
		for start_col, calls in mod.align_varcall_table(pattern, feat_fields, num_features_present, abs_spacing):
				alignment_obj.write("\t" * start_col + "\t".join(calls) + "\n")