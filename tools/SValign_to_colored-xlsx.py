#!/usr/bin/env python3
import sys
import csv
import pandas
import xlsxwriter

if len(sys.argv) != 3:
    print("\tERROR: please input a SV alignment .tsv file and a colors .tsv. Writes an .xlsx.")
    sys.exit(2)

tsv_file = sys.argv[1]
file_name = tsv_file.rsplit(".", 1)[0]
color_file = sys.argv[2]

### read the tsv file into a dataframe, with empty values corrected
# Loop the data lines
with open(tsv_file, 'r') as temp_f:
    # get No of columns in each line
    col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]

# Generate column names
column_names = [i for i in range(0, max(col_count))]

# read the file into a dataframe
calls_df = pandas.read_csv(tsv_file, sep="\t", header = None, names=column_names, dtype="string").fillna('')

# read the colors into a dictionary
color_dict = {}
with open(color_file, 'r') as color_obj:
    color_list = [line for line in csv.reader(color_obj, delimiter="\t")]
    for row in color_list:
        color_dict[row[0]] = row[1]

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pandas.ExcelWriter(f'{file_name}.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
calls_df.to_excel(writer, sheet_name='Sheet1')

# Get the xlsxwriter workbook and worksheet objects.
workbook = writer.book
worksheet = writer.sheets['Sheet1']

# Get the dimensions of the dataframe.
(max_row, max_col) = calls_df.shape

for row in range(0, max_row):
    for col in range(0, max_col):
        cell_value = calls_df.iloc[row, col]
        color = None
        for key in color_dict:
            if cell_value.startswith(key):
                color = color_dict[key]
                break
        if color is not None:
            format1 = workbook.add_format({'bg_color': color})
            worksheet.write(row+1, col+1, cell_value, format1)
        else:
            worksheet.write(row+1, col+1, cell_value)

# Close the Pandas Excel writer and output the Excel file.
writer.close()