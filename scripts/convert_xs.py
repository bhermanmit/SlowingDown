#!/bin/sh env python

import sys

filename = sys.argv[1]
name = sys.argv[2]

# Read in file
with open(filename, 'r') as fh:
    lines = fh.read().splitlines()

# Organize data
ng = len(lines)
E = []
xs = []
for i in range(ng):
    sline = lines[i].split()
    E.append(sline[0])
    xs.append(sline[1])

# Create output string
filestr = "<xs_data>\n"
filestr += "  <name>{0}</name>\n".format(name)
filestr += "  <E>\n"
for i in range(ng):
    filestr += "    {0}\n".format(E[i])
filestr += "  </E>\n"
filestr += "  <xs>\n"
for i in range(ng):
    filestr += "    {0}\n".format(xs[i])
filestr += "  </xs>\n"
filestr += "</xs_data>\n"

# Write output
outfile = filename.split('/')[-1]
with open("{0}.xml".format(outfile.split('.')[0]), 'w') as fh:
    fh.write(filestr)
