import os
from contactCalc import *;

def calculateScoreOne():
    cwd = os.getcwd()
    files = os.listdir(cwd + "/SKEMPI_Dataset")
    for file in files:
        if file[0] != "." and file[-1] == "b":
            print(file)
            getScore(file[0:4], 10, "A", "B")

def calculateCoronavirus():
    cutoff = float(input("Enter the cutoff distance to calculate for: "))
    files = ["6M0J", "7EKC", "7EKF", "7EKG", "7WBP", "7WBQ"]
    chain1 = input("Receptor chain: ")
    chain2 = input("Ligand chain: ")
    for file in files:
        calculate(file, cutoff, chain1, chain2)

calculateScoreOne()
