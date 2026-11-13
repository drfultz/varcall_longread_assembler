#!/usr/bin/env python3
import sys
import tagged_longread_mod as mod

##### Using an anchors file (_anchors.tsv), align all matching reads from a .feat file at each pattern
##### along the entire contig.  Outputs a SValign.tsv that can be visualized with Excel.

if len(sys.argv)!=3:
	print("\tERROR: please input a .feat file and an anchors file. Prints to STDOUT.", file=sys.stderr)
	sys.exit(2)

feat_file = sys.argv[1]
patterns_file = sys.argv[2]

##### read in the feat file
with open(feat_file, 'r') as feat_obj:
	feat_fields = [",".join(units) for units in mod.read_in_feat_aslol(feat_obj)]
# count number of features from first line
num_features_present = len(feat_fields[0].split(",")[0].split("&"))

##### read in patterns file
with open(patterns_file, 'r') as patterns_obj:
	patterns_line = patterns_obj.readline().strip().split("\t")
	spacings_line = patterns_obj.readline().strip().split("\t")

# Output headers
print("#" + "\t#".join(patterns_line))
print("\t".join(spacings_line))

##### Place the first anchor at the farthest-left reach of its matching reads
first_L_spacings = []

### Iterate the feat field and look for the first anchor pattern's leftward spacings
for feat_line in feat_fields:
	num_L_units = mod.count_left_units(patterns_line[0], feat_line)
	if num_L_units is not None:
		first_L_spacings.append(num_L_units)

# Check that some match was found
if not first_L_spacings:
	sys.exit(f"ERROR: first anchor {patterns_line[0]} matched no reads exactly once.")

### Choose the max spacing as the starting spacing (will still add first spacing field)
abs_spacing = max(first_L_spacings)

##### Output alignment file
for pattern, rel_spacing in zip(patterns_line,spacings_line, strict=True):
	abs_spacing += int(rel_spacing)
	print("\t"*abs_spacing*num_features_present + pattern)
	for start_col, calls in mod.align_varcall_table(pattern, feat_fields, num_features_present, abs_spacing):
		print("\t" * start_col + "\t".join(calls))