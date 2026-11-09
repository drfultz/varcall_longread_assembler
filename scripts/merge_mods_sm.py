#!/usr/bin/env python3
import csv
import sys

##### Concatenates the modkit bed lines while also adjusting the reference
##### positions to match a joined consensus.

input_beds_list = list(snakemake.input.beds)
input_cons_list = list(snakemake.input.cons)
output_bed_file = snakemake.output[0]

#input_beds = sys.argv[1]
#output_bed_file = sys.argv[2]

###Loop through the beds and output adjusted lines
with open(output_bed_file, 'w') as output_bed_obj:

	### adjust the positions of each successive bed and output
	counter = 0
	for input_bed, input_cons in zip(input_beds_list, input_cons_list, strict=True):
		with open(input_bed, 'r') as current_bed:
			current_bed_lol = list(csv.reader(current_bed, delimiter = "\t"))

			out_list = []
			for line in current_bed_lol:
				new_line = line
				new_line[0] = snakemake.wildcards.contig
				new_line[1] = str(int(line[1]) + counter)
				new_line[2] = str(int(line[2]) + counter)
				new_line[6] = str(int(line[6]) + counter)
				new_line[7] = str(int(line[7]) + counter)
				out_list.append('\t'.join(new_line) + "\n")

			output_bed_obj.writelines(out_list)

		### Increase the counter by the length of the current region
		with open(input_cons, 'r') as cons_obj:
			cons_obj.readline()
			counter += len(cons_obj.readline().strip())