# download all the pdb files using Biopython

from Bio.PDB import *
import os
import csv

pdbl = PDBList()

cwd = os.getcwd()
with open(cwd + "\\PPI_Dataset\\SI-File-4-protein-protein-test-set-2.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
           pdbl.retrieve_pdb_file(row[0], pdir = cwd + "\\PPI_Dataset", file_format = 'pdb')
        line_count += 1