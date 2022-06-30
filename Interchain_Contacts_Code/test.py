# Used for testing purposes or for short-term code (such as to print out results)

import os
import csv

cwd = os.getcwd()
f = open(cwd + "\\SKEMPI_affinities.txt", mode="a")

map = {}

with open('skempi.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            map[row[0]] = row[8]
        line_count += 1

inOrder = {}
for i in sorted(map):
    inOrder[i] = map[i]

for i in inOrder:
    f.write(i + "," + inOrder[i] + "\n")

f.close()
