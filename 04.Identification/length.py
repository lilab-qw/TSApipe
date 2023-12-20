import sys
import csv

csv.field_size_limit(sys.maxsize)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)  
    sq_col_index = header.index("Sq")  

output_data = []
with open(input_file, "r") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)  
    for row in reader:
        if len(row[sq_col_index]) >= 8 and len(row[sq_col_index]) <= 11:
            output_data.append(row)

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)  
    writer.writerows(output_data)
