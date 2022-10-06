#!/usr/bin/env python3
import sys
import csv
import pandas
import xlsxwriter

##### converts an SV alignment tsv into a colored Excel file.
##### Colors need to be provided in a .tsv with the target string in the first field
##### and the color, in hex/HTML #RRGGBB style, in the second field of each row.

if len(sys.argv) != 3:
	print("\tERROR: please input a SV alignment .tsv file and a colors .tsv.  Writes an .xlsx.")
	sys.exit(2)
	
tsv_file = sys.argv[1]
#output name without extension
file_name = tsv_file.rsplit(".",1)[0]
color_file = sys.argv[2]

### read the tsv file into a dataframe, with empty values corrected
# Loop the data lines
with open(tsv_file, 'r') as temp_f:
    # get No of columns in each line
    col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]

# Generate column names
column_names = [i for i in range(0, max(col_count))]

# read the file into a dataframe
calls_df = pandas.read_csv(tsv_file, sep="\t", header = None, names=column_names, dtype="string")

### read the colors in
with open (color_file, 'r') as color_obj:
	color_list = [line for line in csv.reader(color_obj, delimiter = "\t")]

### Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pandas.ExcelWriter(f'{file_name}.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
calls_df.to_excel(writer, sheet_name='Sheet1')

# Get the xlsxwriter workbook and worksheet objects.
workbook  = writer.book
worksheet = writer.sheets['Sheet1']

# Get the dimensions of the dataframe.
(max_row, max_col) = calls_df.shape
#max_row = 10000
#max_col = 5000

for row in color_list:
	# Add a format.
	format1 = workbook.add_format({'bg_color': row[1]})

	# Apply a conditional format to the required cell range.
	worksheet.conditional_format(1, 1, max_row, max_col,
	                             {'type': 'text',
	                              'criteria': 'begins with',
	                              'value': row[0],
	                              'format': format1})

# Close the Pandas Excel writer and output the Excel file.
writer.close()