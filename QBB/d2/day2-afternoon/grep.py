#!/usr/bin/env python3

import sys

my_file = open( sys.argv[2] )   
pattern = sys.argv[1]

for my_line in my_file:..
    my_line = my_line.rstrip("\n")
    if pattern in my_line:
        print(my_line)

my_file.close()
   
# in terminal: (qb24) cmdb@QuantBio-21 day2-afternoon % ./grep.py FIS1 ../day2-morning/gencode.v46.basic.annotation.gtf | less -S

