import os
from tabulate import tabulate
from Bio.PDB import *

# User-friendly way to get detailed contact info


def calculateWithInput():
    pdbFile = input(
        "Enter the PDB file to calculate contacts for (exclude file extension): ")
    cutoff = float(
        input("Enter the cutoff distance for contacts (in Angstroms): "))
    chain1 = input("Receptor chain: ")
    chain2 = input("Ligand chain: ")
    calculate(pdbFile, cutoff, chain1, chain2)

# Detailed contact data for a given pdf, cutoff, and chains


def calculateHeavy(pdbFile, hisplit, cutoff, specificChains, chain1, chain2, outputFile):
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
            if atoms[j].get_id() != "H" and residuesFirst[i].get_resname() in all:
                if len(atomsFirst[1]) <= i:
                    atomsFirst[1].append([])
                atomsFirst[1][i].append(atoms[j])
    for i in range(len(residuesSecond)):
        atoms = list(residuesSecond[i].get_atoms())
        atomsSecond[0].append(residuesSecond[i].get_resname() +
                              str(residuesSecond[i].get_id()[1]))
        for j in range(len(atoms)):
            if atoms[j].get_id() != "H" and residuesSecond[i].get_resname() in all:
                if len(atomsSecond[1]) <= i:
                    atomsSecond[1].append([])
                atomsSecond[1][i].append(atoms[j])

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypes = [0, 0, 0, 0, 0, 0, 0]
    # favorable (similar HI), nuetral, unfavorable (large HI difference or same-charges)
    contactTypesHI = [0, 0]

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
    # loop thorough residues in second chain
    for i in range(len(atomsSecond[0])):
        count = 0
        contacts = ""
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            # contact = False
            # loop through atoms in second residue
            minDist = 999999999
            for k in range(len(atomsSecond[1][i])):
                # if contact: break
                # loop through atoms in first residue
                for l in range(len(atomsFirst[1][j])):
                    distance = atomsSecond[1][i][k] - atomsFirst[1][j][l]
                    minDist = min(minDist, distance)
                    # HIdiff = abs(hydroIndexesKyte[cAlphaFirst[0][j][:3]]- hydroIndexesKyte[cAlphaSecond[0][i][:3]])
                    # HIscaledDiff = 1-(HIdiff)/(4.5)

                    # contact = True
                    # break
                    # contacts += cAlphaFirst[0][j] + " " + str(int((distance)*1000)/1000) + "\n"
            if (minDist) <= cutoff:
                f.write(str(minDist) + " " + atomsFirst[1][j][l].name + " " + atomsSecond[1]
                        [i][k].name + " " + atomsFirst[0][j] + " " + atomsSecond[0][i] + " " + "\n")
                f.flush()
                
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

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    numFavorable = contactTypes[1] + contactTypes[2] + contactTypes[4] + contactTypes[6] + contactTypes[1] + contactTypes[2] + contactTypes[4] + contactTypes[6]
    numUnfavorable = contactTypes[0] + contactTypes[3] + contactTypes[5] + contactTypes[0] + contactTypes[3] + contactTypes[5]

    '''
    for type in contactTypesHI:
        for num in type:
            f.write(str(num) + " ")
    '''

    f.write(str(numFavorable) + " " + str(numUnfavorable) + " ")
    f.close()


def calculateCA(pdbFile, hisplit, cutoff, specificChains, chain1, chain2, outputFile):
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
    # favorable (similar HI), nuetral, unfavorable (large HI difference or same-charges)
    contactTypesHI = [0, 0]

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
    # loop thorough residues in second chain
    for i in range(len(atomsSecond[0])):
        # loop through residues in first chain
        for j in range(len(atomsFirst[0])):
            distance = atomsSecond[1][i] - atomsFirst[1][j]
            if (distance) <= cutoff:
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


    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    numFavorable = contactTypes[1] + contactTypes[2] + contactTypes[4] + contactTypes[6]
    numUnfavorable = contactTypes[0] + contactTypes[3] + contactTypes[5] + contactTypes[0] + contactTypes[3] + contactTypes[5]

    '''
    for type in contactTypesHI:
        for num in type:
            f.write(str(num) + " ") 
    '''

    # f.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " ")
    f.write(str(numFavorable) + " " + str(numUnfavorable) + " ")
    f.close()
