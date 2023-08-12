from Bio.PDB import *
import math
import csv
import os
import sys

def probs(pdbFile, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    struct = parser.get_structure(
        pdbFile, cwd + "/Combined_Dataset/" + pdbFile + ".pdb")
    model = struct[0]
    f = open(outputFile, 'a')

    # use first 2 chains
    residuesFirst = list(list(model)[0])
    residuesSecond = list(list(model)[1])

    all = {"GLY":0, "ALA":0, "PRO":0, "VAL":0, "ILE":0, "MET":0, "PHE":0, "LEU":0, "TRP":0, "SER":0,
           "THR":0, "CYS":0, "ASN":0, "GLN":0, "TYR":0, "HIS":0, "LYS":0, "ARG":0, "ASP":0, "GLU":0}
    
    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "TYR", "CYS"]
    polar = ["SER", "THR", "ASN", "GLN"]
    charged = ["LYS", "ARG","HIS","ASP","GLU"]

    nums = [0,0,0]

    for i in range(len(residuesFirst)):
        if residuesFirst[i].get_resname() in nonpolar:
            nums[0] += 1
        elif residuesFirst[i].get_resname() in polar:
            nums[1] += 1
        elif residuesFirst[i].get_resname() in charged:
            nums[2] += 1
    for j in range(len(residuesSecond)):
        if residuesSecond[j].get_resname() in nonpolar:
            nums[0] += 1
        elif residuesSecond[j].get_resname() in polar:
            nums[1] += 1
        elif residuesSecond[j].get_resname() in charged:
            nums[2] += 1
    
    for num in nums:
        f.write(str(num) + " ")
    f.write("\n")

'''
cwd = os.getcwd()
o = open(cwd + "/Machine_Learning/output.txt", 'a')
with open(cwd + "/Combined_Dataset/Combined141.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    o.truncate(0)
    for row in csv_reader:
        if line_count > 0:
            # getContactResidues(row[0][0:4], 4.75, cwd + "/Machine_Learning/output.txt")
            probs(row[0][0:4], cwd + "/Machine_Learning/data2.txt")
            # getHydrogenBonds(row[0][0:4], 3.5, cwd + "/Machine_Learning/allFeaturesHICombined.txt")
            # classifyHeavyByAny(row[0][0:4], 4.75, cwd + "/Machine_Learning/allFeaturesHICombined.txt")
            # o.write(row[14] + " " + row[15] + " " + row[16] + " " + row[3] + "\n") # row[14] + " " + row[15] + " " + row[16] + " " + 
            # o.write(row[4] + " " + row[5] + " " + row[6] + " ")
            # o.write("\n")
            o.flush()
            # getHI(row[0][0:4], 4.75, cwd + "/Machine_Learning/allFeaturesHICombined.txt")
            # o.write(row[3] + "\n")
            # o.flush()
        line_count += 1
'''

counts = [0,0,0,0,0,0]

cwd = os.getcwd()
o = open(cwd + "/Machine_Learning/output.txt", 'a')
with open(cwd + "/Machine_Learning/data2.txt") as data:
    lines = data.readlines()
    count = 0
    for line in lines:
        line = line.split(' ')
        apolar = int(line[0])
        polar = int(line[1])
        charged = int(line[2])
        total_res = apolar + polar + charged
        o.write(str((polar*(polar-1))/(total_res*(total_res-1))) + " ")
        o.write(str((polar*apolar*2)/(total_res*(total_res-1))) + " ")
        o.write(str((polar*charged*2)/(total_res*(total_res-1))) + " ")
        o.write(str((apolar*(apolar-1))/(total_res*(total_res-1))) + " ")
        o.write(str((apolar*charged*2)/(total_res*(total_res-1))) + " ")
        o.write(str((charged*(charged-1))/(total_res*(total_res-1))) + "\n")
