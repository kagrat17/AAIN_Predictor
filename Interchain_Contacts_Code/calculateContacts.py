from Bio.PDB import *
import math
import csv
import os
import sys

sys.path.append(os.getcwd() + "/Machine_Learning")

from models import *

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

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER", "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]

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

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP","CYS","TYR"]
    polar = ["SER", "THR", "ASN", "GLN"]
    positive = ["LYS", "ARG","HIS"]
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
    f.write(str(contactTypes[0] + contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " ")
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
    
    # residue identifier, atoms in residue
    atomsFirst = [[], []]
    atomsSecond = [[], []]
    for i in range(len(residuesFirst)):
        if residuesFirst[i].get_resname() not in all or residuesFirst[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesFirst[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA":
                atomsFirst[0].append(residuesFirst[i].get_resname() +
                                    str(residuesFirst[i].get_full_id()[3][1]))
                atomsFirst[1].append(atoms[j])
                break
    for i in range(len(residuesSecond)):
        if residuesSecond[i].get_resname() not in all or residuesSecond[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesSecond[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA":
                atomsSecond[0].append(residuesSecond[i].get_resname() +
                                    str(residuesSecond[i].get_full_id()[3][1]))
                atomsSecond[1].append(atoms[j])
                break

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypesClose = [0, 0, 0, 0, 0, 0, 0]
    contactTypesMid = [0, 0, 0, 0, 0, 0, 0]
    contactTypesFar = [0, 0, 0, 0, 0, 0, 0]
    # >0, <0
    HITypes = [0,0]

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP","TYR","CYS"]
    polar = ["SER", "THR", "ASN", "GLN", ]
    positive = ["LYS", "ARG","HIS"]
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
# Combined_hres includes hydrogens in contacts. Combined_contacts_by_res does not.
def calculateHeavyByRes(pdbFile, cutoff, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/Combined_Dataset/" + pdbFile + ".pdb")
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
        if residuesFirst[i].get_resname() not in all or residuesFirst[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesFirst[i].get_atoms())
        atomsFirst[0].append(residuesFirst[i].get_resname() +
                             str(residuesFirst[i].get_full_id()[3][1]))
        atomsFirst[1].append([])
        atomsFirst[1][len(atomsFirst[1])-1].extend(atoms)
    for i in range(len(residuesSecond)):
        if residuesSecond[i].get_resname() not in all or residuesSecond[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesSecond[i].get_atoms())
        atomsSecond[0].append(residuesSecond[i].get_resname() +
                              str(residuesSecond[i].get_full_id()[3][1]))
        atomsSecond[1].append([])
        atomsSecond[1][len(atomsSecond[1])-1].extend(atoms)

    
    # loop through residues in second chain
    for i in range(len(atomsSecond[0])):
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            minDist = 999999999
            atom1 = ""
            atom2 = ""
            # loop through atoms in second residue
            for k in range(len(atomsSecond[1][i])):
                # loop through atoms in first residue
                for l in range(len(atomsFirst[1][j])):
                    distance = atomsSecond[1][i][k] - atomsFirst[1][j][l]
                    if distance < minDist and atomsFirst[1][j][l].name[0] != "H" and atomsSecond[1][i][k].name[0] != "H":
                        minDist = min(minDist, distance)
                        atom1 = atomsFirst[1][j][l].name
                        atom2 = atomsSecond[1][i][k].name
            if (minDist) <= cutoff:
                f.write(str(minDist) + " " + atom1 + " " + atom2 + " " + atomsFirst[0][j] + " " + atomsSecond[0][i] + " " + "\n")
    
    print(pdbFile + "\n")

# Contacts based on heavy atom any (not residue pair) method. This calculates the contacts (without classifying them). This was used to generate the files in PRODIGY_contacts_by_any
# Includes hydrogens
def calculateHeavyByAny(pdbFile, cutoff, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/Combined_Dataset/" + pdbFile + ".pdb")
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
        if residuesFirst[i].get_resname() not in all or residuesFirst[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesFirst[i].get_atoms())
        atomsFirst[0].append(residuesFirst[i].get_resname() +
                             str(residuesFirst[i].get_full_id()[3][1]))
        atomsFirst[1].append([])
        atomsFirst[1][len(atomsFirst[1])-1].extend(atoms)
    for i in range(len(residuesSecond)):
        if residuesSecond[i].get_resname() not in all or residuesSecond[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesSecond[i].get_atoms())
        atomsSecond[0].append(residuesSecond[i].get_resname() +
                              str(residuesSecond[i].get_full_id()[3][1]))
        atomsSecond[1].append([])
        atomsSecond[1][len(atomsSecond[1])-1].extend(atoms)

    # loop thorough residues in second chain
    for i in range(len(atomsSecond[0])):
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            # loop through atoms in second residue
            for k in range(len(atomsSecond[1][i])):
                # loop through atoms in first residue
                for l in range(len(atomsFirst[1][j])):
                    distance = atomsSecond[1][i][k] - atomsFirst[1][j][l]
                    if (distance) <= cutoff:
                        f.write(str(distance) + " " + atomsFirst[1][j][l].name + " " + atomsSecond[1][i][k].name + " " + atomsFirst[0][j] + " " + atomsSecond[0][i] + " " + "\n")
    # to track progress
    print(pdbFile + "\n")

# Contacts based on heavy atom residue pair method. Classifies contacts into 9 categories (using files in PRODIGY_contacts_by_res)
def classifyHeavyByRes(pdbFile, cutoff, outputFile):

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypes = [0, 0, 0, 0, 0, 0, 0]
    # >0, <0
    HITypes = [0,0]

    tyr = [0, 0, 0, 0, 0, 0, 0]

    nonpolar = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "TYR", "CYS"]
    polar = ["SER", "THR", "ASN", "GLN"]
    positive = ["LYS", "ARG","HIS"]
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
    c = open(cwd + "/Machine_Learning/Combined_contacts_by_res/" + pdbFile + ".txt")
    lines = c.readlines()
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:

            # o.write(distance[3] + " " + distance[4] + " " + distance[0] + "\n")

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
    # o.write(str(contactTypes[0] + contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " ")
    # o.write(str(contactTypes[0] + contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " " + str(sum(contactTypes)) + " ")
    o.write(str(sum(contactTypes)) + " ")
    # o.write(str(tyr[3]) + " " + str(tyr[5]) + " " + str(tyr[6]) + " ")
    o.close()

# Contacts based on heavy atom any (not residue pair) method. Classifies into 3 categories, no hydrogens.
def classifyHeavyByAny(pdbFile, cutoff, outputFile):
    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "/Machine_Learning/Combined_contacts_by_any/" + pdbFile + ".txt")

    polar = ["O", "N", "S"]
    nonpolar = ["C"]

    contactTypes = [0,0,0]
    count = [0,0,0]

    lines = c.readlines()
    countTot = 0
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:

            if distance[1][0] == "C":
                count[0] += 1
            elif distance[1][0] == "N":
                count[1] += 1
            elif distance[1][0] == "O":
                count[2] += 1

            if distance[2][0] == "C":
                count[0] += 1
            elif distance[2][0] == "N":
                count[1] += 1
            elif distance[2][0] == "O":
                count[2] += 1

            countTot += 1

            if distance[1][0] in polar and distance[2][0] in polar:
                contactTypes[0] += 1
            elif (distance[1][0] in polar and distance[2][0] in nonpolar) or (distance[1][0] in nonpolar and distance[2][0] in polar):
                contactTypes[1] += 1
            elif distance[1][0] in nonpolar and distance[2][0] in nonpolar:
                contactTypes[2] += 1
        else:
            break

    # o.write(str(count[0]) + " " + str(count[1]) + " " + str(count[2]) + " ")
    # o.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " ")
    # o.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(sum(contactTypes)) + " ")
    o.write(str(sum(contactTypes)) + " ")
    o.close()

# Gets total number of each type of residue in the complex
def getTotalResidues(pdbFile, outputFile):
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
    
    for i in range(len(residuesFirst)):
        if residuesFirst[i].get_resname() in all:
            all[residuesFirst[i].get_resname()] += 1
    for j in range(len(residuesSecond)):
        if residuesSecond[j].get_resname() in all:
            all[residuesSecond[j].get_resname()] += 1
    
    for res in all:
        f.write(str(all[res]) + " ")
    f.write("\n")

# Gets number of each type of residue involved in the contact
def getContactResidues(pdbFile, cutoff, outputFile):
    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "/Machine_Learning/Combined_contacts_by_res/" + pdbFile + ".txt")

    all = {"GLY":0, "ALA":0, "PRO":0, "VAL":0, "ILE":0, "MET":0, "PHE":0, "LEU":0, "TRP":0, "SER":0,
           "THR":0, "CYS":0, "ASN":0, "GLN":0, "TYR":0, "HIS":0, "LYS":0, "ARG":0, "ASP":0, "GLU":0}

    lines = c.readlines()
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:
            all[distance[3][0:3]] += 1
            all[distance[4][0:3]] += 1

    for res in all:
        o.write(str(all[res]) + " ")
    o.write("\n")
    o.close()

def getHydrogenBonds(pdbFile, cutoff, outputFile):
    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "/Machine_Learning/Combined_contacts_by_any/" + pdbFile + ".txt")

    contactTypes = [0,0,0]

    lines = c.readlines()
    countTot = 0
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:
            countTot += 1

            if distance[1][0] == "N" and distance[2][0] == "N":
                contactTypes[0] += 1
            elif (distance[1][0] == "N" and distance[2][0] == "O") or (distance[1][0] == "O" and distance[2][0] == "N"):
                contactTypes[1] += 1
            elif distance[1][0] == "O" and distance[2][0] == "O":
                contactTypes[2] += 1

    o.write(str(contactTypes[0] + contactTypes[1] + contactTypes[2]) + " ")
    o.close()

def getHI(pdbFile, cutoff, outputFile):
    # >0, <0
    HITypes = [0,0]

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

    total = 0

    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "/Machine_Learning/Combined_contacts_by_res/" + pdbFile + ".txt")
    lines = c.readlines()
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:

            # o.write(distance[3] + " " + distance[4] + " " + distance[0] + "\n")

            HIdiff = abs(hydroIndexesKyte[distance[3][:3]]- hydroIndexesKyte[distance[4][:3]])
            HIscaledDiff = 1-(HIdiff)/(4.5)
            total += HIscaledDiff

            if HIscaledDiff >= 0:
                HITypes[0] += 1
            elif HIscaledDiff < 0:
                HITypes[1] += 1

    # o.write(str(HITypes[0]) + " " + str(HITypes[1]) + " ")
    o.write(str(round(total,2)) + " ")
    o.close()

def atomByRes(pdbFile, cutoff, outputFile):
    cwd = os.getcwd()
    o = open(outputFile, 'a')
    c = open(cwd + "/Machine_Learning/Combined_contacts_by_any/" + pdbFile + ".txt")

    polar = ["O", "N", "S"]
    nonpolar = ["C"]

    nonpolarRes = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "TYR", "CYS"]
    polarRes = ["SER", "THR", "ASN", "GLN"]
    positiveRes = ["LYS", "ARG","HIS"]
    negativeRes = ["ASP", "GLU"]

    contactTypes = [0,0,0]
    count = [0,0,0]

    lines = c.readlines()
    countTot = 0
    for distance in lines:
        distance = distance.split(' ')
        if float(distance[0]) <= cutoff:
            countTot += 1

            nonpolarFirst = distance[3][:3] in nonpolarRes
            nonpolarSecond = distance[4][:3] in nonpolarRes
            polarFirst = distance[3][:3] in polarRes
            polarSecond = distance[4][:3] in polarRes
            positiveFirst = distance[3][:3] in positiveRes
            positiveSecond = distance[4][:3] in positiveRes
            negativeFirst = distance[3][:3] in negativeRes
            negativeSecond = distance[4][:3] in negativeRes

            if (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                if distance[1][0] == "C":
                    count[0] += 1
                elif distance[1][0] == "N":
                    count[1] += 1
                elif distance[1][0] == "O":
                    count[2] += 1

                if distance[2][0] == "C":
                    count[0] += 1
                elif distance[2][0] == "N":
                    count[1] += 1
                elif distance[2][0] == "O":
                    count[2] += 1
                
                if distance[1][0] in polar and distance[2][0] in polar:
                    contactTypes[0] += 1
                elif (distance[1][0] in polar and distance[2][0] in nonpolar) or (distance[1][0] in nonpolar and distance[2][0] in polar):
                    contactTypes[1] += 1
                elif distance[1][0] in nonpolar and distance[2][0] in nonpolar:
                    contactTypes[2] += 1
        else:
            break
    
    o.write(str(count[0]) + " " + str(count[1]) + " " + str(count[2]) + " ")
    # o.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " ")
    # o.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(sum(contactTypes)) + " ")
    # o.write(str(sum(contactTypes)) + " ")
    o.close()

def allAtoms(pdbFile, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(
        pdbFile, cwd + "/Combined_Dataset/" + pdbFile + ".pdb")
    model = struct[0]
    f = open(outputFile, 'a')

    residuesFirst = list(list(model)[0])
    residuesSecond = list(list(model)[1])

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER",
           "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]

    # residue identifier, atoms in residue
    atomsList = []
    for i in range(len(residuesFirst)):
        if residuesFirst[i].get_resname() not in all or residuesFirst[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesFirst[i].get_atoms())
        atomsList.extend(atoms)
    for i in range(len(residuesSecond)):
        if residuesSecond[i].get_resname() not in all or residuesSecond[i].get_full_id()[3][2] != " ":
            continue
        atoms = list(residuesSecond[i].get_atoms())
        atomsList.extend(atoms)
    
    count = [0,0,0]
    for atom in atomsList:
        if atom.name[0] == "C":
            count[0] += 1
        elif atom.name[0] == "N":
            count[1] += 1
        elif atom.name[0] == "O":
            count[2] += 1
    
    o.write(str(count[0]) + " " + str(count[1]) + " " + str(count[2]) + " ")


""" cwd = os.getcwd()
o = open(cwd + "/Machine_Learning/data.txt", 'a')
with open(cwd + "/Combined_Dataset/Combined171.csv") as csv_file:
    o.truncate(0)
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            classifyHeavyByRes(row[0][0:4], 4.75, cwd + "/Machine_Learning/data.txt")
            classifyHeavyByAny(row[0][0:4], 4.75, cwd + "/Machine_Learning/data.txt")
            getHydrogenBonds(row[0][0:4], 3.5, cwd + "/Machine_Learning/data.txt")
            o.write("\n")
            # o.write(row[4] + " " + row[5] + " " + row[6] + " " + row[2] + "\n")
            o.flush()
        line_count += 1 """

def sortFirst(line):
    return float(line.split(' ')[0])

# sort the contacts by any by distance (to make looping through them faster)
""" cwd = os.getcwd()
o = open(cwd + "/Machine_Learning/data.txt", 'a')
with open(cwd + "/Combined_Dataset/Combined171.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            # os.rename(row[0][0:4] + '.pdb', cwd + '/PDBbind_PPI_used/' + row[0][0:4] + '.pdb')
            c = open(cwd + "/Machine_Learning/Combined_contacts_by_any/" + row[0][0:4] + ".txt")
            lines = c.readlines()
            lines.sort(key = sortFirst)
            c.close()
            o = open(cwd + "/Machine_Learning/Combined_contacts_by_any/" + row[0][0:4] + ".txt", 'a')
            o.truncate(0)
            for line in lines:
                o.write(line)
            o.close()
            # classifyHeavyByRes(row[0][0:4], 5.5, cwd + "/Machine_Learning/ppi_data.txt")
            # o.write("\n")
            # o.flush()
            # o.write(row[14] + " " + row[15] + " " + row[16] + " " + row[3] + "\n")
        line_count += 1 """


# Optimize cutoff

""" cwd = os.getcwd()
o = open(cwd + "/Machine_Learning/data2.txt", 'a')
cutoff = 0
while cutoff <= 20:
    o.truncate(0)
    with open(cwd + "/Combined_Dataset/Combined141.csv") as csv_file:
        csv_reader = csv.reader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count > 0:
                classifyHeavyByRes(row[0][0:4], cutoff, cwd + "/Machine_Learning/data2.txt")
                # o.write(row[14] + " " + row[15] + " " + row[16] + " " + row[3] + "\n")
                o.write(row[3] + "\n") # row[4] + " " + row[5] + " " + row[6] + " " + 
                # o.write("\n")
                o.flush()
            line_count += 1
        train(3,81,141,0.5)
        cutoff += 0.25 """

# loop through dataset        

cwd = os.getcwd()
o = open(cwd + "/Machine_Learning/data2.txt", 'a')
with open(cwd + "/Combined_Dataset/Combined141.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    o.truncate(0)
    for row in csv_reader:
        if line_count > 0:
            classifyHeavyByAny(row[0][0:4], 10, cwd + "/Machine_Learning/data2.txt")
            o.write("\n")
            o.flush()
            """ 
            # getContactResidues(row[0][0:4], 4.75, cwd + "/Machine_Learning/output.txt")
            classifyHeavyByRes(row[0][0:4], 4.75, cwd + "/Machine_Learning/allFeaturesHIProd.txt")
            getHydrogenBonds(row[0][0:4], 3.5, cwd + "/Machine_Learning/allFeaturesHIProd.txt")
            classifyHeavyByAny(row[0][0:4], 4.75, cwd + "/Machine_Learning/allFeaturesHIProd.txt")
            # o.write(row[14] + " " + row[15] + " " + row[16] + " " + row[3] + "\n") # row[14] + " " + row[15] + " " + row[16] + " " + 
            o.write(row[14] + " " + row[15] + " " + row[16] + " ")
            getHI(row[0][0:4], 4.75, cwd + "/Machine_Learning/allFeaturesHICombined.txt")
            o.write(row[3] + "\n")
            o.flush() """
        line_count += 1