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

def calculateCoronavirus():
    cutoff = float(input("Enter the cutoff distance to calculate for: "))
    files = ["6M0J", "7EKC", "7EKF", "7EKG", "7WBP", "7WBQ"]
    chain1 = input("Receptor chain: ")
    chain2 = input("Ligand chain: ")
    for file in files:
        calculate(file, cutoff, chain1, chain2)

calculateSKEMPI()
