import os
from tabulate import tabulate
from Bio.PDB import *

def calculateWithInput():
    pdbFile = input("Enter the PDB file to calculate contacts for (exclude file extension): ")
    cutoff = float(input("Enter the cutoff distance for contacts (in Angstroms): "))
    chain1 = input("Receptor chain: ")
    chain2 = input("Ligand chain: ")
    calculate(pdbFile, cutoff, chain1, chain2)

def calculate(pdbFile, cutoff, chain1, chain2):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    struct = parser.get_structure(pdbFile, cwd + "/PDB_Files/" + pdbFile + ".pdb")
    model = struct.get_models()
    f = open(cwd + "/RBD-ACE2_Contacts/" + str(cutoff) + "-Angstroms-Rose/" + pdbFile + "_Contacts_" + str(cutoff) + ".txt", mode="w")
    #f = open(cwd + "/RBD-ACE2_Contacts/dataSet.txt", mode="a")

    models = list(model)
    residuesFirst = model[chain1]
    residuesSecond = model[chain2]
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

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP"]
    polar = ["SER", "THR", "CYS", "ASN", "GLN", "TYR", "HIS"]
    positive = ["LYS", "ARG"]
    negative = ["ASP", "GLU"]

    hydroIndexesKyte = {
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

    hydroIndexesEngelman = {
        "ALA": 1.60,
        "ARG": -12.30,
        "ASN": -4.80,
        "ASP":	-9.20,
        "CYS":	2.00,
        "GLN":	-4.10,
        "GLU":	-8.20,
        "GLY":	1.00,
        "HIS":	-3.00,
        "ILE":	3.10,
        "LEU":	2.80,
        "LYS":	-8.80,
        "MET":	3.40,
        "PHE":	3.70,
        "PRO":	-0.20,
        "SER":	0.60,
        "THR":  1.20,
        "TRP":	1.90,
        "TYR":	-0.70,
        "VAL":	2.60
    }

    hydroIndexesCornette = {
        "ALA": 0.20,
        "ARG": 1.40,
        "ASN": -0.50,
        "ASP":	-3.10,
        "CYS":	4.10,
        "GLN":	-2.80,
        "GLU":	-1.80,
        "GLY":	0.00,
        "HIS":	0.50,
        "ILE":	4.80,
        "LEU":	5.70,
        "LYS":	-3.10,
        "MET":	4.20,
        "PHE":	4.40,
        "PRO":	-2.20,
        "SER":	-0.50,
        "THR":  -1.90,
        "TRP":	1.00,
        "TYR":	3.20,
        "VAL":	4.70
    }

    hydroIndexesEisenberg = {
        "ALA": 0.62,
        "ARG": -2.53,
        "ASN": -0.78,
        "ASP":	-0.90,
        "CYS":	0.29,
        "GLN":	-0.85,
        "GLU":	-0.74,
        "GLY":	0.48,
        "HIS":	-0.40,
        "ILE":	1.38,
        "LEU":	1.06,
        "LYS":	-1.50,
        "MET":	0.64,
        "PHE":	1.19,
        "PRO":	0.12,
        "SER":	-0.18,
        "THR":  -0.05,
        "TRP":	0.81,
        "TYR":	0.26,
        "VAL":	1.08
    }

    hydroIndexesHoppWoods = {
        "ALA": -0.50,
        "ARG": 3.00,
        "ASN": 0.20,
        "ASP":	3.00,
        "CYS":	-1.00,
        "GLN":	0.20,
        "GLU":	3.00,
        "GLY":	0.00,
        "HIS":	-0.50,
        "ILE":	-1.80,
        "LEU":	-1.80,
        "LYS":	3.00,
        "MET":	-1.30,
        "PHE":	-2.50,
        "PRO":	0.00,
        "SER":	0.30,
        "THR":  -0.40,
        "TRP":	-3.40,
        "TYR":	-2.30,
        "VAL":	-1.50
    }

    hydroIndexesJanin = {
        "ALA": 0.30,
        "ARG": -1.30,
        "ASN": -0.50,
        "ASP":	-0.60,
        "CYS":	0.90,
        "GLN":	-0.70,
        "GLU":	-0.70,
        "GLY":	0.30,
        "HIS":	-0.10,
        "ILE":	0.70,
        "LEU":	0.50,
        "LYS":	-1.80,
        "MET":	0.40,
        "PHE":	0.50,
        "PRO":	-0.30,
        "SER":	-0.10,
        "THR":  -0.20,
        "TRP":	0.30,
        "TYR":	-0.40,
        "VAL":	0.60
    }

    hydroIndexesRose = {
        "ALA": 0.74,
        "ARG": 0.64,
        "ASN": 0.63,
        "ASP": 0.62,
        "CYS": 0.91,
        "GLN": 0.62,
        "GLU": 0.62,
        "GLY": 0.72,
        "HIS": 0.78,
        "ILE": 0.88,
        "LEU": 0.85,
        "LYS": 0.52,
        "MET": 0.85,
        "PHE": 0.88,
        "PRO": 0.64,
        "SER": 0.66,
        "THR": 0.70,
        "TRP": 0.85,
        "TYR": 0.76,
        "VAL": 0.86
    }

    countTot = 0
    hiScoreMult = 0.0
    hiScoreAdd = 0.0
    for i in range(len(cAlphaSecond[0])):
        count = 0
        contacts = ""
        for j in range(len(cAlphaFirst[0])):
            if(cAlphaSecond[1][i] - cAlphaFirst[1][j]) <= cutoff:
                count += 1
                countTot += 1
                contacts += cAlphaFirst[0][j] + " " + str(int((cAlphaSecond[1][i] - cAlphaFirst[1][j])*1000)/1000) + "\n"
                # classifying contacts
                nonpolarFirst = cAlphaFirst[0][j][:3] in nonpolar
                nonpolarSecond = cAlphaSecond[0][i][:3] in nonpolar
                polarFirst = cAlphaFirst[0][j][:3] in polar
                polarSecond = cAlphaSecond[0][i][:3] in polar
                positiveFirst = cAlphaFirst[0][j][:3] in positive
                positiveSecond = cAlphaSecond[0][i][:3] in positive
                negativeFirst = cAlphaFirst[0][j][:3] in negative
                negativeSecond = cAlphaSecond[0][i][:3] in negative
                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    contactTypes[0] += 1
                elif negativeSecond and positiveFirst or positiveFirst and negativeSecond:
                    contactTypes[1] += 1
                elif (polarFirst or polarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                    contactTypes[2] += 1
                elif (nonpolarFirst or nonpolarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                    contactTypes[3] += 1
                elif polarFirst and polarSecond:
                    contactTypes[4] += 1
                elif (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                    contactTypes[5] += 1
                elif nonpolarFirst and nonpolarSecond:
                    contactTypes[6] += 1
                
                # custom hydropathy scoring
                h1 = hydroIndexesRose[cAlphaSecond[0][i][:3]]
                h2 = hydroIndexesRose[cAlphaFirst[0][j][:3]]
                #h1 = abs(h1)
                #h2 = abs(h2)

                score = h1*h2/(cAlphaSecond[1][i] - cAlphaFirst[1][j])
                #score = abs(score)
                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    score = -1*abs(score)
                
                """
                elif nonpolarFirst and polarSecond or polarFirst and nonpolarSecond:
                    score *= -1
                elif nonpolarFirst and (positiveSecond or negativeSecond) or (positiveFirst or negativeFirst) and nonpolarSecond:
                    score *= -1
                """
                hiScoreMult += score
                
                """"
                score = (h1+h2)/(cAlphaSecond[1][i] - cAlphaFirst[1][j])
                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    score = -1*abs(score)
                hiScoreAdd += score
                """

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
    #f.write("Hydropathy index score (addition): " + str(hiScoreAdd) + "\n\n")
    f.write(str(len(data[0])) + " Chain " + chain1 + " Contact Residues: ")
    for i in range(len(data[0])):
        f.write(str(data[0][i]) + "(" + str(data[1][i]) + ") ")
    f.write("\n")
    f.write(str(len(aceRes)) + " Chain " + chain1 + " Contact Residues: ")
    for c in aceRes:
        f.write(str(c.split(" ")[0]) + "(" + str(aceRes[c]) + ") ")
    f.write("\n\n")
    f.write(tabulate((dataInRows), headers=["Chain " + chain1, "Chain " + chain2]))
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

def getScore(pdbFile, cutoff, chain1, chain2):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    struct = parser.get_structure(pdbFile, cwd + "/PRODIGYdataset/" + pdbFile + ".pdb")
    model = struct[0]
    #f = open(cwd + "/RBD-ACE2_Contacts/" + str(cutoff) + "-Angstroms-Rose/" + pdbFile + "_Contacts_" + str(cutoff) + ".txt", mode="w")
    f = open(cwd + "/RBD-ACE2_Contacts/dataSet.txt", mode="a")

    residuesFirst = list(model["A"])
    residuesSecond = list(model["B"])
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

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP"]
    polar = ["SER", "THR", "CYS", "ASN", "GLN", "TYR", "HIS"]
    positive = ["LYS", "ARG"]
    negative = ["ASP", "GLU"]

    # maxDiff = 9
    hydroIndexesKyte = {
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
    # maxDiff = 16
    hydroIndexesEngelman = {
        "ALA": 1.60,
        "ARG": -12.30,
        "ASN": -4.80,
        "ASP":	-9.20,
        "CYS":	2.00,
        "GLN":	-4.10,
        "GLU":	-8.20,
        "GLY":	1.00,
        "HIS":	-3.00,
        "ILE":	3.10,
        "LEU":	2.80,
        "LYS":	-8.80,
        "MET":	3.40,
        "PHE":	3.70,
        "PRO":	-0.20,
        "SER":	0.60,
        "THR":  1.20,
        "TRP":	1.90,
        "TYR":	-0.70,
        "VAL":	2.60
    }
    # maxDiff = 8.8
    hydroIndexesCornette = {
        "ALA": 0.20,
        "ARG": 1.40,
        "ASN": -0.50,
        "ASP":	-3.10,
        "CYS":	4.10,
        "GLN":	-2.80,
        "GLU":	-1.80,
        "GLY":	0.00,
        "HIS":	0.50,
        "ILE":	4.80,
        "LEU":	5.70,
        "LYS":	-3.10,
        "MET":	4.20,
        "PHE":	4.40,
        "PRO":	-2.20,
        "SER":	-0.50,
        "THR":  -1.90,
        "TRP":	1.00,
        "TYR":	3.20,
        "VAL":	4.70
    }
    # maxDiff = 3.91
    hydroIndexesEisenberg = {
        "ALA": 0.62,
        "ARG": -2.53,
        "ASN": -0.78,
        "ASP":	-0.90,
        "CYS":	0.29,
        "GLN":	-0.85,
        "GLU":	-0.74,
        "GLY":	0.48,
        "HIS":	-0.40,
        "ILE":	1.38,
        "LEU":	1.06,
        "LYS":	-1.50,
        "MET":	0.64,
        "PHE":	1.19,
        "PRO":	0.12,
        "SER":	-0.18,
        "THR":  -0.05,
        "TRP":	0.81,
        "TYR":	0.26,
        "VAL":	1.08
    }
    # maxDiff = 6.4
    hydroIndexesHoppWoods = {
        "ALA": -0.50,
        "ARG": 3.00,
        "ASN": 0.20,
        "ASP":	3.00,
        "CYS":	-1.00,
        "GLN":	0.20,
        "GLU":	3.00,
        "GLY":	0.00,
        "HIS":	-0.50,
        "ILE":	-1.80,
        "LEU":	-1.80,
        "LYS":	3.00,
        "MET":	-1.30,
        "PHE":	-2.50,
        "PRO":	0.00,
        "SER":	0.30,
        "THR":  -0.40,
        "TRP":	-3.40,
        "TYR":	-2.30,
        "VAL":	-1.50
    }
    # maxDiff = 2.7
    hydroIndexesJanin = {
        "ALA": 0.30,
        "ARG": -1.30,
        "ASN": -0.50,
        "ASP":	-0.60,
        "CYS":	0.90,
        "GLN":	-0.70,
        "GLU":	-0.70,
        "GLY":	0.30,
        "HIS":	-0.10,
        "ILE":	0.70,
        "LEU":	0.50,
        "LYS":	-1.80,
        "MET":	0.40,
        "PHE":	0.50,
        "PRO":	-0.30,
        "SER":	-0.10,
        "THR":  -0.20,
        "TRP":	0.30,
        "TYR":	-0.40,
        "VAL":	0.60
    }
    # maxDiff = 0.39
    hydroIndexesRose = {
        "ALA": 0.74,
        "ARG": 0.64,
        "ASN": 0.63,
        "ASP": 0.62,
        "CYS": 0.91,
        "GLN": 0.62,
        "GLU": 0.62,
        "GLY": 0.72,
        "HIS": 0.78,
        "ILE": 0.88,
        "LEU": 0.85,
        "LYS": 0.52,
        "MET": 0.85,
        "PHE": 0.88,
        "PRO": 0.64,
        "SER": 0.66,
        "THR": 0.70,
        "TRP": 0.85,
        "TYR": 0.76,
        "VAL": 0.86
    }

    countTot = 0
    hiScoreMult = 0.0
    minDist = 5.0
    maxDist = 12.0
    for i in range(len(cAlphaSecond[0])):
        count = 0
        contacts = ""
        for j in range(len(cAlphaFirst[0])):
            if(cAlphaSecond[1][i] - cAlphaFirst[1][j]) <= maxDist:
                count += 1
                countTot += 1
                contacts += cAlphaFirst[0][j] + " " + str(int((cAlphaSecond[1][i] - cAlphaFirst[1][j])*1000)/1000) + "\n"
                # classifying contacts
                nonpolarFirst = cAlphaFirst[0][j][:3] in nonpolar
                nonpolarSecond = cAlphaSecond[0][i][:3] in nonpolar
                polarFirst = cAlphaFirst[0][j][:3] in polar
                polarSecond = cAlphaSecond[0][i][:3] in polar
                positiveFirst = cAlphaFirst[0][j][:3] in positive
                positiveSecond = cAlphaSecond[0][i][:3] in positive
                negativeFirst = cAlphaFirst[0][j][:3] in negative
                negativeSecond = cAlphaSecond[0][i][:3] in negative
                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    contactTypes[0] += 1
                elif negativeSecond and positiveFirst or positiveFirst and negativeSecond:
                    contactTypes[1] += 1
                elif (polarFirst or polarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                    contactTypes[2] += 1
                elif (nonpolarFirst or nonpolarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                    contactTypes[3] += 1
                elif polarFirst and polarSecond:
                    contactTypes[4] += 1
                elif (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                    contactTypes[5] += 1
                elif nonpolarFirst and nonpolarSecond:
                    contactTypes[6] += 1
                
                # custom hydropathy scoring
                h1 = hydroIndexesRose[cAlphaSecond[0][i][:3]]
                h2 = hydroIndexesRose[cAlphaFirst[0][j][:3]]

                score = (-abs(h1-h2)/(0.195) + 1.0) * (1-min(1,max(0,(abs((cAlphaSecond[1][i]-cAlphaFirst[1][j])-minDist))/(maxDist-minDist))))

                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    score = -1*abs(score)
                elif positiveFirst and negativeSecond or negativeFirst and positiveSecond:
                    score = abs(score)

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
    
    f.write(str(hiScoreMult) + "\n")

    f.flush()
    f.close()