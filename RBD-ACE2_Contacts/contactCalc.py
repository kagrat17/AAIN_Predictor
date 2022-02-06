import os
from tkinter import W
from tabulate import tabulate
import numpy as np
from Bio.PDB import *

pdbFile = input("Enter the PDB file to calculate contacts for (exclude file extension): ")

cwd = os.getcwd()
f = open(cwd + "/RBD-ACE2_Contacts/Contact_Data/" + pdbFile + "_Contacts.txt", mode="w")
parser = PDBParser(PERMISSIVE=True, QUIET=True)
struct = parser.get_structure(pdbFile, cwd + "/PDB_Files/" + pdbFile + ".pdb")
model = struct.get_models()

models = list(model)
chains = list(models[0].get_chains())
residuesA = list(chains[0].get_residues())
residuesE = list(chains[1].get_residues())
# residue identifier, atom
cAlphaA = [[], []]
cAlphaE = [[], []]
for i in range(len(residuesA)):
    atoms = list(residuesA[i].get_atoms())
    for j in range(len(atoms)):
        if atoms[j].get_id() == "CA":
            cAlphaA[0].append(residuesA[i].get_resname() +
                              str(residuesA[i].get_id()[1]))
            cAlphaA[1].append(atoms[j])
            break
for i in range(len(residuesE)):
    atoms = list(residuesE[i].get_atoms())
    for j in range(len(atoms)):
        if atoms[j].get_id() == "CA":
            cAlphaE[0].append(residuesE[i].get_resname() +
                              str(residuesE[i].get_id()[1]))
            cAlphaE[1].append(atoms[j])
            break

# RBD residue, count contacts with ACE2, list of contacted residues
data = [[], [], []]

countTot = 0
for i in range(len(cAlphaE[0])):
    count = 0
    contacts = ""
    for j in range(len(cAlphaA[0])):
        if(cAlphaE[1][i] - cAlphaA[1][j]) <= 7:
            count += 1
            countTot += 1
            contacts += cAlphaA[0][j] + " " + \
                str(int((cAlphaE[1][i] - cAlphaA[1][j])*1000)/1000) + "\n"
    contacts = contacts[:len(contacts)-2]
    if count > 0:
        data[0].append(cAlphaE[0][i])
        data[1].append(count)
        data[2].append(contacts)

aceRes = dict()
for i in range(len(data[2])):
    res = data[2][i].split("\n")
    for s in res:
        if s.split(" ")[0] in aceRes:
            aceRes[s.split(" ")[0]] += 1
        else:
            aceRes[s.split(" ")[0]] = 1

dataInRows = []
for i in range(len(data[0])):
    arr = [str(data[0][i]) + "(" + str(data[1][i]) + ")", data[2][i]]
    dataInRows.append(arr)

f.write(pdbFile + ": " + struct.header["name"] + "\n\n")
f.write("Total Contacts: " + str(countTot) + "\n")
f.write("Cutoff Distance: 7 Angstroms" + "\n\n")
f.write(str(len(data[0])) + " RBD Contact Residues: ")
for i in range(len(data[0])):
    f.write(str(data[0][i]) + "(" + str(data[1][i]) + ") ")
f.write("\n")
f.write(str(len(aceRes)) + " ACE2 Contact Residues: ")
for c in aceRes:
    f.write(str(c.split(" ")[0]) + "(" + str(aceRes[c]) + ") ")
f.write("\n\n")
f.write(tabulate((dataInRows), headers=["RBD Residue", "ACE2 Residues"]))

f.flush()
f.close()