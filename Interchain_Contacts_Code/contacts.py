import os
from tabulate import tabulate
from Bio.PDB import *

# User-friendly way to get detailed contact info
def calculateWithInput():
    pdbFile = input("Enter the PDB file to calculate contacts for (exclude file extension): ")
    cutoff = float(input("Enter the cutoff distance for contacts (in Angstroms): "))
    chain1 = input("Receptor chain: ")
    chain2 = input("Ligand chain: ")
    calculate(pdbFile, cutoff, chain1, chain2)

# Detailed contact data for a given pdf, cutoff, and chains
def calculate(pdbFile, cutoff, split, specificChains, chain1, chain2, outputFile):
    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    # struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    struct = parser.get_structure(pdbFile, cwd + "/PRODIGY_Dataset/" + pdbFile + ".pdb")
    model = struct[0]
    f = open(outputFile, 'a')

    # use given chains
    if(specificChains):
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

    # residue identifier, atom
    cAlphaFirst = [[], []]
    cAlphaSecond = [[], []]
    for i in range(len(residuesFirst)):
        atoms = list(residuesFirst[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA" and residuesFirst[i].get_resname() in all:
                cAlphaFirst[0].append(residuesFirst[i].get_resname() +
                                str(residuesFirst[i].get_id()[1]))
                cAlphaFirst[1].append(atoms[j])
                break
    for i in range(len(residuesSecond)):
        atoms = list(residuesSecond[i].get_atoms())
        for j in range(len(atoms)):
            if atoms[j].get_id() == "CA" and residuesSecond[i].get_resname() in all:
                cAlphaSecond[0].append(residuesSecond[i].get_resname() +
                                str(residuesSecond[i].get_id()[1]))
                cAlphaSecond[1].append(atoms[j])
                break

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    contactTypes = [[0,0,0,0,0,0,0], [0,0,0,0,0,0,0]]
    # favorable (similar HI), nuetral, unfavorable (large HI difference or same-charges)
    contactTypesHI = [[0,0], [0,0]]

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

    

    countTot = 0
    for i in range(len(cAlphaSecond[0])):
        count = 0
        contacts = ""
        for j in range(len(cAlphaFirst[0])):
            distance = cAlphaSecond[1][i] - cAlphaFirst[1][j]
            if(distance) <= split:
                HIdiff = abs(hydroIndexesKyte[cAlphaFirst[0][j][:3]]- hydroIndexesKyte[cAlphaSecond[0][i][:3]])
                HIscaledDiff = 1-(HIdiff)/(4.5)
                count += 1
                countTot += 1
                contacts += cAlphaFirst[0][j] + " " + str(int((distance)*1000)/1000) + "\n"
                # classifying contacts
                nonpolarFirst = cAlphaFirst[0][j][:3] in nonpolar
                nonpolarSecond = cAlphaSecond[0][i][:3] in nonpolar
                polarFirst = cAlphaFirst[0][j][:3] in polar
                polarSecond = cAlphaSecond[0][i][:3] in polar
                positiveFirst = cAlphaFirst[0][j][:3] in positive
                positiveSecond = cAlphaSecond[0][i][:3] in positive
                negativeFirst = cAlphaFirst[0][j][:3] in negative
                negativeSecond = cAlphaSecond[0][i][:3] in negative

                expression = 0 if distance <= split else 1

                if positiveFirst and positiveSecond or negativeFirst and negativeSecond:
                    contactTypes[expression][0] += 1
                elif negativeSecond and positiveFirst or positiveFirst and negativeSecond:
                    contactTypes[expression][1] += 1
                elif (polarFirst or polarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                    contactTypes[expression][2] += 1
                elif (nonpolarFirst or nonpolarSecond) and (positiveFirst or positiveSecond or negativeFirst or negativeSecond):
                    contactTypes[expression][3] += 1
                elif polarFirst and polarSecond:
                    contactTypes[expression][4] += 1
                elif (polarFirst or polarSecond) and (nonpolarFirst or nonpolarSecond):
                    contactTypes[expression][5] += 1
                elif nonpolarFirst and nonpolarSecond:
                    contactTypes[expression][6] += 1
                
                if positiveFirst and positiveSecond or negativeFirst and negativeSecond or HIscaledDiff <= 0:
                    contactTypesHI[expression][1] += 1
                else:
                    contactTypesHI[expression][0] += 1

    # same charge, opposite charge, charged-polar, charged-nonpolar, polar-polar, polar-nonpolar, nonpolar-nonpolar
    numFavorable = contactTypes[0][1] + contactTypes[0][2] + contactTypes[0][4] + contactTypes[0][6] + contactTypes[1][1] + contactTypes[1][2] + contactTypes[1][4] + contactTypes[1][6]
    numUnfavorable = contactTypes[0][0] + contactTypes[0][3] + contactTypes[0][5] + contactTypes[1][0] + contactTypes[1][3] + contactTypes[1][5] 

    f.write(str(contactTypesHI[0][0] + contactTypesHI[0][1]) + " ")

    '''
    for type in contactTypesHI:
        for num in type:
            f.write(str(num) + " ") 
    '''

    # f.write(str(contactTypes[0]) + " " + str(contactTypes[1]) + " " + str(contactTypes[2]) + " " + str(contactTypes[3]) + " " + str(contactTypes[4]) + " " + str(contactTypes[5]) + " " + str(contactTypes[6]) + " ")

    f.flush()
    f.close()