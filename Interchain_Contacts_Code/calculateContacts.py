from Bio.PDB import *
import math
import csv
import os
import sys

# Contacts based on CA method. Classifies contacts into 9 categories.
def calculateCA(pdbFile, cutoff, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    model = struct[0]
    f = open(outputFile, 'a')

    # use given chains
    if (specificChains):
        residuesFirst = []
        residuesSecond = []
        for i in range(len(chain1)):
            residuesFirst += list(model[chain1[i]])
        for i in range(len(chain2)):
            residuesSecond += list(model[chain2[i]])
    # use first 2 chains
    else:
        residuesFirst = list(list(model)[0])
        residuesSecond = list(list(model)[1])

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER",
           "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]

    # residue identifier, CA atom
    atomsFirst = [[], []]
    atomsSecond = [[], []]
    for i in range(len(residuesFirst)):
        atoms = list(residuesFirst[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA" and residuesFirst[i].get_resname() in all:
                atomsFirst[0].append(residuesFirst[i].get_resname() +
                                     str(residuesFirst[i].get_id()[1]))
                atomsFirst[1].append(atoms[j])
                break
    for i in range(len(residuesSecond)):
        atoms = list(residuesSecond[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA" and residuesSecond[i].get_resname() in all:
                atomsSecond[0].append(residuesSecond[i].get_resname() +
                                      str(residuesSecond[i].get_id()[1]))
                atomsSecond[1].append(atoms[j])
                break

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypes = [0, 0, 0, 0, 0, 0, 0]
    # >0, <0
    HITypes = [0,0]

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP"]
    polar = ["SER", "THR", "CYS", "ASN", "GLN", "TYR", "HIS"]
    positive = ["LYS", "ARG"]
    negative = ["ASP", "GLU"]

    # maxDiff = 9
    hydroIndexesKyte = {
        "ALA": 1.80,
        "ARG": -4.50,
        "ASN": -3.50,
        "ASP": -3.50,
        "CYS":	2.50,
        "GLN": -3.50,
        "GLU": -3.50,
        "GLY": -0.40,
        "HIS": -3.20,
        "ILE":	4.50,
        "LEU":	3.80,
        "LYS": -3.90,
        "MET":	1.90,
        "PHE":	2.80,
        "PRO":	1.60,
        "SER": -0.80,
        "THR": -0.70,
        "TRP": -0.90,
        "TYR": -1.30,
        "VAL":	4.20
    }
    
    countTot = 0
    # loop through residues in second chain
    for i in range(len(atomsSecond[0])):
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            distance = atomsSecond[1][i] - atomsFirst[1][j]
            if (distance) <= cutoff:
                countTot += 1

                HIdiff = abs(hydroIndexesKyte[atomsFirst[0][j][:3]]- hydroIndexesKyte[atomsSecond[0][i][:3]])
                HIscaledDiff = 1-(HIdiff)/(4.5)
                
                # classifying contacts
                nonpolarFirst = atomsFirst[0][j][:3] in nonpolar
                nonpolarSecond = atomsSecond[0][i][:3] in nonpolar
                polarFirst = atomsFirst[0][j][:3] in polar
                polarSecond = atomsSecond[0][i][:3] in polar
                positiveFirst = atomsFirst[0][j][:3] in positive
                positiveSecond = atomsSecond[0][i][:3] in positive
                negativeFirst = atomsFirst[0][j][:3] in negative
                negativeSecond = atomsSecond[0][i][:3] in negative

                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    contactTypes[0] += 1
                elif negativeSecond and positiveFirst or negativeFirst and positiveSecond:
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
                
                if HIscaledDiff >= 0:
                    HITypes[0] += 1
                elif HIscaledDiff < 0:
                    HITypes[1] += 1


    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    f.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " " + str(HITypes[0]) + " " + str(HITypes[1]) + "\n")
    f.close()

# Contacts based on CA method. Classifies contacts into 21 categories (3 distance ranges, 7 contact types).
def calculateCASplitDistances(pdbFile, close, mid, far, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    model = struct[0]
    f = open(outputFile, 'a')

    # use given chains
    if (specificChains):
        residuesFirst = []
        residuesSecond = []
        for i in range(len(chain1)):
            residuesFirst += list(model[chain1[i]])
        for i in range(len(chain2)):
            residuesSecond += list(model[chain2[i]])
    # use first 2 chains
    else:
        residuesFirst = list(list(model)[0])
        residuesSecond = list(list(model)[1])

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER",
           "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]

    # residue identifier, CA atom
    atomsFirst = [[], []]
    atomsSecond = [[], []]
    for i in range(len(residuesFirst)):
        atoms = list(residuesFirst[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA" and residuesFirst[i].get_resname() in all:
                atomsFirst[0].append(residuesFirst[i].get_resname() +
                                     str(residuesFirst[i].get_id()[1]))
                atomsFirst[1].append(atoms[j])
                break
    for i in range(len(residuesSecond)):
        atoms = list(residuesSecond[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA" and residuesSecond[i].get_resname() in all:
                atomsSecond[0].append(residuesSecond[i].get_resname() +
                                      str(residuesSecond[i].get_id()[1]))
                atomsSecond[1].append(atoms[j])
                break

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypesClose = [0, 0, 0, 0, 0, 0, 0]
    contactTypesMid = [0, 0, 0, 0, 0, 0, 0]
    contactTypesFar = [0, 0, 0, 0, 0, 0, 0]
    # >0, <0
    HITypes = [0,0]

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP"]
    polar = ["SER", "THR", "CYS", "ASN", "GLN", "TYR", "HIS"]
    positive = ["LYS", "ARG"]
    negative = ["ASP", "GLU"]
    
    countTot = 0
    # loop through residues in second chain
    for i in range(len(atomsSecond[0])):
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            distance = atomsSecond[1][i] - atomsFirst[1][j]
            if distance <= far:
                countTot += 1
                
                # classifying contacts
                nonpolarFirst = atomsFirst[0][j][:3] in nonpolar
                nonpolarSecond = atomsSecond[0][i][:3] in nonpolar
                polarFirst = atomsFirst[0][j][:3] in polar
                polarSecond = atomsSecond[0][i][:3] in polar
                positiveFirst = atomsFirst[0][j][:3] in positive
                positiveSecond = atomsSecond[0][i][:3] in positive
                negativeFirst = atomsFirst[0][j][:3] in negative
                negativeSecond = atomsSecond[0][i][:3] in negative

                if distance < close:
                    if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                        contactTypesClose[0] += 1
                    elif negativeSecond and positiveFirst or positiveFirst and negativeSecond:
                        contactTypesClose[1] += 1
                    elif (polarFirst or polarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                        contactTypesClose[2] += 1
                    elif (nonpolarFirst or nonpolarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                        contactTypesClose[3] += 1
                    elif polarFirst and polarSecond:
                        contactTypesClose[4] += 1
                    elif (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                        contactTypesClose[5] += 1
                    elif nonpolarFirst and nonpolarSecond:
                        contactTypesClose[6] += 1
                
                elif distance < mid:
                    if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                        contactTypesMid[0] += 1
                    elif negativeSecond and positiveFirst or positiveFirst and negativeSecond:
                        contactTypesMid[1] += 1
                    elif (polarFirst or polarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                        contactTypesMid[2] += 1
                    elif (nonpolarFirst or nonpolarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                        contactTypesMid[3] += 1
                    elif polarFirst and polarSecond:
                        contactTypesMid[4] += 1
                    elif (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                        contactTypesMid[5] += 1
                    elif nonpolarFirst and nonpolarSecond:
                        contactTypesMid[6] += 1
                
                elif distance <= far:
                    if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                        contactTypesFar[0] += 1
                    elif negativeSecond and positiveFirst or positiveFirst and negativeSecond:
                        contactTypesFar[1] += 1
                    elif (polarFirst or polarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                        contactTypesFar[2] += 1
                    elif (nonpolarFirst or nonpolarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                        contactTypesFar[3] += 1
                    elif polarFirst and polarSecond:
                        contactTypesFar[4] += 1
                    elif (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                        contactTypesFar[5] += 1
                    elif nonpolarFirst and nonpolarSecond:
                        contactTypesFar[6] += 1

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    f.write(str(contactTypesClose[0]) + " " + str(contactTypesClose[1]) + " " + str(contactTypesClose[2]) + " " + str(contactTypesClose[3]) + " " + str(contactTypesClose[4]) + " " + str(contactTypesClose[5]) + " " + str(contactTypesClose[6]) + " ")
    f.write(str(contactTypesMid[0]) + " " + str(contactTypesMid[1]) + " " + str(contactTypesMid[2]) + " " + str(contactTypesMid[3]) + " " + str(contactTypesMid[4]) + " " + str(contactTypesMid[5]) + " " + str(contactTypesMid[6]) + " ")
    f.write(str(contactTypesFar[0]) + " " + str(contactTypesFar[1]) + " " + str(contactTypesFar[2]) + " " + str(contactTypesFar[3]) + " " + str(contactTypesFar[4]) + " " + str(contactTypesFar[5]) + " " + str(contactTypesFar[6]) + "\n")
    f.close()

# Contacts based on heavy atom residue pair method. This calculates the contacts (without classifying them). This was used to generate the files in PRODIGY_contacts_by_res
def calculateHeavyByRes(pdbFile, cutoff, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    model = struct[0]
    f = open(outputFile, 'a')
    f.truncate(0)

    # use given chains
    if (specificChains):
        residuesFirst = []
        residuesSecond = []
        for i in range(len(chain1)):
            residuesFirst += list(model[chain1[i]])
        for i in range(len(chain2)):
            residuesSecond += list(model[chain2[i]])
    # use first 2 chains
    else:
        residuesFirst = list(list(model)[0])
        residuesSecond = list(list(model)[1])

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER",
           "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]

    # residue identifier, atoms in residue
    atomsFirst = [[], []]
    atomsSecond = [[], []]
    for i in range(len(residuesFirst)):
        atoms = list(residuesFirst[i].get_atoms())
        atomsFirst[0].append(residuesFirst[i].get_resname() +
                             str(residuesFirst[i].get_id()[1]))
        for j in range(len(atoms)):
            if atoms[j].get_id()[0] != "H" and residuesFirst[i].get_resname() in all:
                if len(atomsFirst[1]) <= i:
                    atomsFirst[1].append([])
                atomsFirst[1][i].append(atoms[j])
    for i in range(len(residuesSecond)):
        atoms = list(residuesSecond[i].get_atoms())
        atomsSecond[0].append(residuesSecond[i].get_resname() +
                              str(residuesSecond[i].get_id()[1]))
        for j in range(len(atoms)):
            if atoms[j].get_id()[0] != "H" and residuesSecond[i].get_resname() in all:
                if len(atomsSecond[1]) <= i:
                    atomsSecond[1].append([])
                atomsSecond[1][i].append(atoms[j])

    # loop thorough residues in second chain
    for i in range(len(atomsSecond[0])):
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            # loop through atoms in second residue
            minDist = 999999999
            for k in range(len(atomsSecond[1][i])):
                # loop through atoms in first residue
                for l in range(len(atomsFirst[1][j])):
                    distance = atomsSecond[1][i][k] - atomsFirst[1][j][l]
                    minDist = min(minDist, distance)
            if (minDist) <= cutoff:
                f.write(str(minDist) + " " + atomsFirst[1][j][l].name + " " + atomsSecond[1]
                        [i][k].name + " " + atomsFirst[0][j] + " " + atomsSecond[0][i] + " " + "\n")

# Contacts based on heavy atom any (not residue pair) method. This calculates the contacts (without classifying them). This was used to generate the files in PRODIGY_contacts_by_any
def calculateHeavyByAny(pdbFile, cutoff, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/PDBbind_PPI_used/" + pdbFile + ".ent.pdb")
    model = struct[0]
    f = open(outputFile, 'a')
    f.truncate(0)

    # use given chains
    if (specificChains):
        residuesFirst = []
        residuesSecond = []
        for i in range(len(chain1)):
            residuesFirst += list(model[chain1[i]])
        for i in range(len(chain2)):
            residuesSecond += list(model[chain2[i]])
    # use first 2 chains
    else:
        residuesFirst = list(list(model)[0])
        residuesSecond = list(list(model)[1])

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER",
           "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]

    # residue identifier, atoms in residue
    atomsFirst = [[], []]
    atomsSecond = [[], []]
    for i in range(len(residuesFirst)):
        atomsFirst = list(residuesFirst[i].get_atoms())
    for i in range(len(residuesSecond)):
        atomsSecond = list(residuesSecond[i].get_atoms())

    # loop through atoms in second residue
    for k in range(len(atomsSecond)):
        # loop through atoms in first residue
        for l in range(len(atomsFirst)):
            distance = atomsSecond[k] - atomsFirst[l]
            if distance <= cutoff:
                f.write(str(distance) + " " + atomsFirst[l].name + " " + atomsSecond[k].name)

# Contacts based on heavy atom residue pair method. Classifies contacts into 9 categories (using files in PRODIGY_contacts_by_res)
def classifyHeavyByRes(pdbFile, cutoff, outputFile):

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypes = [0, 0, 0, 0, 0, 0, 0]
    # >0, <0
    HITypes = [0,0]

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP"]
    polar = ["SER", "THR", "CYS", "ASN", "GLN", "TYR", "HIS"]
    positive = ["LYS", "ARG"]
    negative = ["ASP", "GLU"]

    # maxDiff = 9
    hydroIndexesKyte = {
        "ALA": 1.80,
        "ARG": -4.50,
        "ASN": -3.50,
        "ASP": -3.50,
        "CYS":	2.50,
        "GLN": -3.50,
        "GLU": -3.50,
        "GLY": -0.40,
        "HIS": -3.20,
        "ILE":	4.50,
        "LEU":	3.80,
        "LYS": -3.90,
        "MET":	1.90,
        "PHE":	2.80,
        "PRO":	1.60,
        "SER": -0.80,
        "THR": -0.70,
        "TRP": -0.90,
        "TYR": -1.30,
        "VAL":	4.20
    }

    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "\\Machine_Learning\\PRODIGY_contacts_by_res\\" + pdbFile + ".txt")
    lines = c.readlines()
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:

            HIdiff = abs(hydroIndexesKyte[distance[3][:3]]- hydroIndexesKyte[distance[4][:3]])
            HIscaledDiff = 1-(HIdiff)/(4.5)
                
            nonpolarFirst = distance[3][:3] in nonpolar
            nonpolarSecond = distance[4][:3] in nonpolar
            polarFirst = distance[3][:3] in polar
            polarSecond = distance[4][:3] in polar
            positiveFirst = distance[3][:3] in positive
            positiveSecond = distance[4][:3] in positive
            negativeFirst = distance[3][:3] in negative
            negativeSecond = distance[4][:3] in negative

            if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                contactTypes[0] += 1
            elif negativeSecond and positiveFirst or negativeFirst and positiveSecond:
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

            if HIscaledDiff >= 0:
                HITypes[0] += 1
            elif HIscaledDiff < 0:
                HITypes[1] += 1

    o.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " " + str(HITypes[0]) + " " + str(HITypes[1]) + "\n")
    o.close()

# Contacts based on heavy atom any (not residue pair) method. Calculates total number of contacts.
def classifyHeavyByAny(pdbFile, cutoff, outputFile):
    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "\\Machine_Learning\\PRODIGY_contacts_by_any\\" + pdbFile + ".txt")

    polar = ["O", "N", "S"]
    nonpolar = ["C"]

    contactTypes = [0,0,0]

    lines = c.readlines()
    countTot = 0
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:
            countTot += 1

            if distance[1][0] in polar and distance[2][0] in polar:
                contactTypes[0] += 1
            elif (distance[1][0] in polar and distance[2][0] in nonpolar) or (distance[1][0] in nonpolar and distance[2][0] in polar):
                contactTypes[1] += 1
            elif distance[1][0] in nonpolar and distance[2][0] in nonpolar:
                contactTypes[2] += 1

    o.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " ")
    o.close()

''' Prodigy
cwd = os.getcwd()
o = open(cwd + "\\Machine_Learning\\prodigy_data.txt", 'a')
with open(cwd + "\\PRODIGY_Dataset\\PRODIGY_dataset.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0 and row[2] == "ER":
            classifyHeavyByAny(row[0][0:4], 6, cwd + "\\Machine_Learning\\prodigy_data.txt")
            o.write(row[14] + " " + row[15] + " " + row[16] + " " + row[3] + "\n")
            o.flush()
        line_count += 1
'''

cwd = os.getcwd()
o = open(cwd + "\\Machine_Learning\\ppi_data.txt", 'a')
with open(cwd + "\\PDBbind_PPI_used\\set_4.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            calculateHeavyByAny(row[0], 6, False, "A", "B", cwd + "\\Machine_Learning\\PPI_contacts_by_any\\" + row[0] + ".txt")
        line_count += 1