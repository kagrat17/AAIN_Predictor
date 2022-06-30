# Retreive data and run algorithms (such as score prediction) on it

import os
from contactCalc import *;

def calculateSKEMPI():
    with open('SKEMPI_affinities.txt') as file:
        lines = file.readlines()
        for line in lines:
            line = line.split(",")
            pdb = line[0].split("_")
            getScore(pdb[0], pdb[1], pdb[2])

def calculateProdigy():
    cwd = os.getcwd()
    for file in os.listdir(cwd + "\\PRODIGYdataset"):
        if file[0] != ".":
            getScore(file[0:4],"A","B")

calculateSKEMPI()
