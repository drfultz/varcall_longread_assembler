#!/usr/bin/env python3
import sys
import csv
import tagged_longread_mod as mod

##### Using a regular expression string from command line, align matching reads from a .feat file at this point.
##### Outputs a SValign.tsv .

if len(sys.argv) != 3:
	print("\tERROR: please input a .feat file and a regular expression search string. "
		"Best to place the search pattern in quotes. Prints to STDOUT.", file=sys.stderr)
	sys.exit(2)
	
features_file = sys.argv[1]
search_pattern = sys.argv[2]


##### read in the feat file, then rewind to collect readnames in the same order
with open(features_file, 'r') as feat_obj:
	feat_fields = [",".join(units) for units in mod.read_in_feat_aslol(feat_obj)]
	feat_obj.seek(0)
	readnames = [row[0] for row in csv.reader(feat_obj, delimiter="\t")]

# count number of features from the first non-empty read
num_features_present = len(next(f for f in feat_fields if f).split(",")[0].split("&"))


##### Find reads matching the pattern exactly once and measure units left of it
matched_reads = []
for readname, feat_line in zip(readnames, feat_fields, strict=True):
	num_L_units = mod.count_left_units(search_pattern, feat_line)
	if num_L_units is not None:
		matched_reads.append((readname, num_L_units, feat_line.split(",")))

# Check that some match was found
if not matched_reads:
	sys.exit(f"ERROR: {search_pattern} matched no reads exactly once.")

### Print the coverage level
print(f"coverage = {len(matched_reads)}", file=sys.stderr)

### The anchor column is the farthest-left reach of the matching reads
max_L_offset = max(num_L_units for _, num_L_units, _ in matched_reads)


##### Output the feature calls horizontally
### print the search pattern above the anchor column where it matches
print("\t" * (1 + max_L_offset * num_features_present) + search_pattern)
### print the feature calls for each read, one per line, aligned to the pattern
for readname, num_L_units, units in matched_reads:
	start_col = 1 + (max_L_offset - num_L_units) * num_features_present
	print(readname + "\t" * start_col + "\t".join(units).replace("&", "\t"))
