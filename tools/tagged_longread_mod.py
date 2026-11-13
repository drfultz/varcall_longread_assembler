"""Shared helpers for the tools/ scripts: reading .feat files and aligning reads on a variant-call pattern."""

import csv
import re


def read_in_feat_aslol(file_obj):
	"""read in a feat file (variant calls only), parsed into a list of lists,
	removing empty units at either end and '%' marks"""
	feats = []
	for tsv_line in csv.reader(file_obj, delimiter = "\t"):

		split_units = tsv_line[1].split(",")
		trimmed_units = []
		### loop through each unit
		for u, unit in enumerate(split_units):
			## if first or last unit
			if u == 0 or u == (len(split_units)-1):
				# keep unit if any call is not Unmapped
				for subunit in unit.split("&"):
					if not subunit.endswith("Unmapped"):
						trimmed_units.append(unit.replace("%",""))
						break
			## keep all internal units
			else:
				trimmed_units.append(unit.replace("%",""))

		# put trimmed_units lists into list
		feats.append(trimmed_units)

	return feats

def count_left_units(search_pattern, line):
	"""Return the number of units left of search_pattern in line, or None unless it matches exactly once."""
	matches = list(re.finditer(search_pattern, line))
	if len(matches) != 1:
		return None
	return line.count(",", 0, matches[0].start())

def align_varcall_table(search_pattern, varcalls, num_feat, spacing):
	"""find pattern matches in variant call fields, and return a list of tuples of the spacing before the first field
	followed by flattened variant call fields."""
	aligned_fields = []
	max_L_offset = 0

	### gather any lines that match pattern and track max offset
	for line in varcalls:
		## Measure to left of pattern to position it
		num_L_units = count_left_units(search_pattern, line)

		if num_L_units is not None:
			# Keep track of the most number of L units
			max_L_offset = max(max_L_offset, num_L_units)
			# Add the num L and a flattened list of all feature calls
			aligned_fields.append((num_L_units, line.replace("&", ",").split(",")))

	### Check for sufficient L padding on non-zero spacings
	if spacing > 0 and max_L_offset > spacing:
		raise ValueError(f"Reads extend {max_L_offset} units left of {search_pattern}, more than its"
			f" spacing of {spacing}. Increase the initial buffer.")

	#prep and return
	anchor_col = spacing if spacing > 0 else max_L_offset
	expanded_fields = [((anchor_col - num_L_units) * num_feat, calls) for num_L_units, calls in aligned_fields]
	return expanded_fields