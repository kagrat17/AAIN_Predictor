from Bio.PDB import *
import math
import csv
import os
import sys

cwd = os.getcwd()
f = open(cwd + "/Machine_Learning/output.txt", "r")
data = f.readlines()

o = open(cwd + "/Machine_Learning/data2.txt", 'a')

def sortAIC(line):
    return float(line.split('\t')[1])

data.sort(key = sortAIC)

for l in data:
    o.write(l)