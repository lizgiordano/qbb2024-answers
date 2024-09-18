#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )   

for my_line in my_file:
    if "##" in my_line:
            continue
    line = line.rstrip("\.")
    line = line.rstrip("\n")
    line = line.lstrip"\."
    fields = my_line.split("\t")
    print(fields[0], fields[3], fields[4], fields[8])

my_file.close()

# for each line in my data
# separate by tabs/column (split)
# select out columns I want
# print
    
