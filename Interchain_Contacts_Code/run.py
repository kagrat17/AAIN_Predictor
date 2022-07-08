# Retreive data and run algorithms (such as score prediction) on it

import os
import csv
from contacts import *;

# calculate hydropathy scoring for the SKEMPI dataset
def calculateSKEMPI():
    with open('SKEMPI_affinities.txt') as file:
        lines = file.readlines()
        for line in lines:
            line = line.split(",")
            pdb = line[0].split("_")
            getScore(pdb[0], pdb[1], pdb[2])

# calculate hydropathy scoring for the PRODIGY dataset
def calculateProdigy():
    cwd = os.getcwd()
    for file in os.listdir(cwd + "\\PRODIGY_Dataset"):
        if file[0] != ".":
            getScore(file[0:4],"A","B")

# print input parameters and experimental affinity to data file for machine learning
def getData():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\prodigy_data.txt", 'a')
    with open(cwd + "\\PRODIGY_Dataset\\PRODIGY_dataset.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                calculate(row[0][0:4],10,"A","B",cwd + "\\Machine_Learning\\prodigy_data.txt")
                f.write(str(row[3]) + "\n")
                f.flush()
            line_count += 1
    f.close()

getData()
