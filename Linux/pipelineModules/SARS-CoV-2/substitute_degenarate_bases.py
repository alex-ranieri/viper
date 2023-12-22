#!/usr/bin/env python
# coding: utf-8

# Importing libraries
from Bio import SeqIO
import re, sys

# Reading input
ivar_consensus = sys.argv[1]

# Output file name
clean_segment = sys.argv[2]

handle = open(ivar_consensus, "r")

sequences={}

for record in SeqIO.parse(handle, "fasta") :
    #replace any non-GATC characters in your sequences with N
    sequence = re.sub('[^GATC]', "N", str(record.seq).upper())
    sequences[sequence]=record.id

output_file=open(clean_segment,"w+")

for sequence in sequences:
    output_file.write(">"+sequences[sequence]+"\n"+sequence+"\n")
output_file.close()
handle.close()
