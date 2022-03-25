import os
from tabulate import tabulate
from Bio.PDB import *

pdbFile = input("Enter the PDB file to calculate contacts for (exclude file extension): ")
cutoff = float(input("Enter the cutoff distance for contacts (in Angstroms): "))

cwd = os.getcwd()
parser = PDBParser(PERMISSIVE=True, QUIET=True)
struct = parser.get_structure(pdbFile, cwd + "/PDB_Files/" + pdbFile + ".pdb")
model = struct.get_models()
f = open(cwd + "/RBD-ACE2_Contacts/" + str(cutoff) + "-Angstroms/" + pdbFile + "_Contacts_" + str(cutoff) + ".txt", mode="w")

models = list(model)
chains = list(models[0].get_chains())
residuesFirst = list(chains[0].get_residues())
residuesSecond = list(chains[1].get_residues())
# residue identifier, atom
cAlphaFirst = [[], []]
cAlphaSecond = [[], []]
for i in range(len(residuesFirst)):
    atoms = list(residuesFirst[i].get_atoms())
    for j in range(len(atoms)):
        if atoms[j].get_id() == "CA":
            cAlphaFirst[0].append(residuesFirst[i].get_resname() +
                              str(residuesFirst[i].get_id()[1]))
            cAlphaFirst[1].append(atoms[j])
            break
for i in range(len(residuesSecond)):
    atoms = list(residuesSecond[i].get_atoms())
    for j in range(len(atoms)):
        if atoms[j].get_id() == "CA":
            cAlphaSecond[0].append(residuesSecond[i].get_resname() +
                              str(residuesSecond[i].get_id()[1]))
            cAlphaSecond[1].append(atoms[j])
            break

# same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
contactTypes = [0,0,0,0,0,0,0]
# Second chain residue, count contacts with first chain, list of contacted residues
data = [[], [], []]

nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE"]
polar = ["SER", "THR", "CYS", "ASN", "GLN", "LEU", "TYR", "TRP", "HIS"]
positive = ["LYS", "ARG"]
negative = ["ASP", "GLU"]

hydroIndexes = {
    "ALA": 1.80,
    "ARG": -4.50,
    "ASN": -3.50,
    "ASP":	-3.50,
    "CYS":	2.50,
    "GLN":	-3.50,
    "GLU":	-3.50,
    "GLY":	-0.40,
    "HIS":	-3.20,
    "ILE":	4.50,
    "LEU":	3.80,
    "LYS":	-3.90,
    "MET":	1.90,
    "PHE":	2.80,
    "PRO":	1.60,
    "SER":	-0.80,
    "THR":	-0.70,
    "TRP":	-0.90,
    "TYR":	-1.30,
    "VAL":	4.20
}

countTot = 0
hiScoreMult = 0.0
for i in range(len(cAlphaSecond[0])):
    count = 0
    contacts = ""
    for j in range(len(cAlphaFirst[0])):
        if(cAlphaSecond[1][i] - cAlphaFirst[1][j]) <= cutoff:
            count += 1
            countTot += 1
            contacts += cAlphaFirst[0][j] + " " + str(int((cAlphaSecond[1][i] - cAlphaFirst[1][j])*1000)/1000) + "\n"
            # classifying contacts
            nonpolarA = cAlphaFirst[0][j][:3] in nonpolar
            nonpolarE = cAlphaSecond[0][i][:3] in nonpolar
            polarA = cAlphaFirst[0][j][:3] in polar
            polarE = cAlphaSecond[0][i][:3] in polar
            positiveA = cAlphaFirst[0][j][:3] in positive
            positiveE = cAlphaSecond[0][i][:3] in positive
            negativeA = cAlphaFirst[0][j][:3] in negative
            negativeE = cAlphaSecond[0][i][:3] in negative
            if positiveA and positiveE or negativeA and negativeE:
                contactTypes[0] += 1
            elif negativeE and positiveA or positiveA and negativeE:
                contactTypes[1] += 1
            elif (polarA or polarE) and (positiveA or positiveE or negativeA or negativeE):
                contactTypes[2] += 1
            elif (nonpolarA or nonpolarE) and (positiveA or positiveE or negativeA or negativeE):
                contactTypes[3] += 1
            elif polarA and polarE:
                contactTypes[4] += 1
            elif (polarA or polarE) and (nonpolarA or nonpolarE):
                contactTypes[5] += 1
            elif nonpolarA and nonpolarE:
                contactTypes[6] += 1
            
            # custom hydropathy scoring
            h1 = hydroIndexes[cAlphaSecond[0][i][:3]]
            h2 = hydroIndexes[cAlphaFirst[0][j][:3]]
            score = h1*h2/(cAlphaSecond[1][i] - cAlphaFirst[1][j])
            if (positiveA and positiveE or negativeA and negativeE) and score > 0:
                score *= -1
            elif (negativeE and positiveA or positiveA and negativeE) and score < 0:
                score *= -1
            hiScoreMult += score

    contacts = contacts[:len(contacts)-2]
    if count > 0:
        data[0].append(cAlphaSecond[0][i])
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
f.write("Cutoff Distance: " + str(cutoff) + " Angstroms" + "\n")
f.write("Hydropathy index score (multiplication): " + str(hiScoreMult) + "\n\n")
f.write(str(len(data[0])) + " Chain " + chains[1].id + " Contact Residues: ")
for i in range(len(data[0])):
    f.write(str(data[0][i]) + "(" + str(data[1][i]) + ") ")
f.write("\n")
f.write(str(len(aceRes)) + " Chain " + chains[0].id + " Contact Residues: ")
for c in aceRes:
    f.write(str(c.split(" ")[0]) + "(" + str(aceRes[c]) + ") ")
f.write("\n\n")
f.write(tabulate((dataInRows), headers=["Chain " + chains[1].id, "Chain " + chains[0].id]))
f.write("\n\n")

f.write("Contact Types\n")
f.write("Same charge: " + str(contactTypes[0]) + "\n")
f.write("Opposite charge: " + str(contactTypes[1]) + "\n")
f.write("Charged-polar: " + str(contactTypes[2]) + "\n")
f.write("Charged-nonpolar: " + str(contactTypes[3]) + "\n")
f.write("Polar-polar: " + str(contactTypes[4]) + "\n")
f.write("Polar-nonpolar: " + str(contactTypes[5]) + "\n")
f.write("Nonpolar-nonpolar: " + str(contactTypes[6]) + "\n")

f.flush()
f.close()