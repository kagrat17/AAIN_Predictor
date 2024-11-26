from Bio.PDB import *
import os

# downloads PDB file into pdbFiles directory to be used for prediction
def getPDB(name):
    cwd = os.getcwd()
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(name,file_format="pdb",pdir=cwd+"/pdbFiles")

# calculates the interface numbers of the six amino acids used in the model
def getAAINs(pdbFile, chains):
    # counts of tyr, gly, ser, arg, val, and ile at the interface
    counts = {"TYR":0, "GLY":0, "SER":0, "ARG":0, "VAL":0, "ILE":0}

    cwd = os.getcwd()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    struct = parser.get_structure(
        pdbFile, cwd + "/pdbFiles/pdb" + pdbFile + ".ent") # this can also be a .pdb file.
    model = struct[0]

    # use given chains
    if (chains != None):
        residuesFirst = []
        residuesSecond = []
        for i in range(len(chains[0])):
            residuesFirst += list(model[chains[0][i]])
        for i in range(len(chains[1])):
            residuesSecond += list(model[chains[1][i]])
    # use first 2 chains
    else:
        residuesFirst = list(list(model)[0])
        residuesSecond = list(list(model)[1])

    all = ["GLY", "ALA", "PRO", "VAL", "ILE", "MET", "PHE", "LEU", "TRP", "SER",
           "THR", "CYS", "ASN", "GLN", "TYR", "HIS", "LYS", "ARG", "ASP", "GLU"]
    res = ["TYR", "GLY", "SER", "ARG", "VAL", "ILE"]

    # get all the atoms in the chains
    atomsFirst = [[], []] # residue identifier, atoms in residue
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
            minDist = 99999 # arbitrary large value
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
            if (minDist) <= 4.75:
                if atomsFirst[0][j][0:3] in res:
                    counts[atomsFirst[0][j][0:3]] += 1
                if atomsSecond[0][i][0:3] in res:
                    counts[atomsSecond[0][i][0:3]] += 1
    
    return counts

# takes a dictionary of the six interface numbers and predicts Delta G
def predictFromINs(INs):
    return -0.1535*INs["TYR"] - 0.1288*INs["GLY"] - 0.0840*INs["SER"] - 0.0805*INs["ARG"] - 0.0684*INs["VAL"]  + 0.073*INs["ILE"]  - 6.46

# takes in a pdb file, downloads it, and returns a prediction
def predict(name, chains=None):
    return predictFromINs(getAAINs(name,chains))

# Example prediction
# getPDB("1AHW")
# print(predict("1AHW",("AB","C")))
