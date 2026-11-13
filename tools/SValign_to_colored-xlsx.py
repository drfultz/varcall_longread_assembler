#!/usr/bin/env python3
import sys
import csv
import xlsxwriter

if len(sys.argv) != 3:
	print("\tERROR: please input a SV alignment .tsv file and a colors .tsv. Writes an .xlsx.", file=sys.stderr)
	sys.exit(2)

tsv_file = sys.argv[1]
file_name = tsv_file.rsplit(".", 1)[0]
color_file = sys.argv[2]

### read the alignment rows (ragged rows are fine)
with open(tsv_file, 'r') as tsv_obj:
	rows = [line.rstrip("\n").split("\t") for line in tsv_obj]
max_cols = max(len(values) for values in rows)

### read the colors into a dictionary, longest keys first so the most specific prefix wins
with open(color_file, 'r') as color_obj:
	color_dict = {row[0]: row[1] for row in csv.reader(color_obj, delimiter="\t") if row}
color_keys = sorted(color_dict, key=len, reverse=True)

### write the workbook directly (constant_memory streams each row to disk as it's finished)
workbook = xlsxwriter.Workbook(f'{file_name}.xlsx', {'constant_memory': True})
worksheet = workbook.add_worksheet('Sheet1')
header_format = workbook.add_format({'bold': True, 'border': 1, 'align': 'center'})

# header row of column numbers
worksheet.write_row(0, 1, range(max_cols), header_format)

format_cache = {}   # color -> format, one per color
value_formats = {}  # cell value -> format (or None), computed once per unique value

for row, values in enumerate(rows):
	worksheet.write(row+1, 0, row, header_format)  # index column
	for col, cell_value in enumerate(values):
		if not cell_value:
			continue
		if cell_value not in value_formats:
			color = next((color_dict[k] for k in color_keys if cell_value.startswith(k)), None)
			if color is not None and color not in format_cache:
				format_cache[color] = workbook.add_format({'bg_color': color})
			value_formats[cell_value] = format_cache.get(color)
		worksheet.write_string(row+1, col+1, cell_value, value_formats[cell_value])

workbook.close()