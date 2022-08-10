# Retreive data and run algorithms (such as score prediction) on it

import os
import sys
import csv
import math
from contacts import *;

sys.path.append(os.getcwd() + "\\Machine_Learning")

from models import *;

# print input parameters and experimental affinity to data file from the PRODIGY dataset for machine learning
def getProdigyData():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\prodigy_data_2.txt", 'a')
    with open(cwd + "\\PRODIGY_Dataset\\PRODIGY_dataset.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                calculate(row[0][0:4],10,5,True,"A","B",cwd + "\\Machine_Learning\\prodigy_data_2.txt")
                f.write(str(row[3]) + "\n")
                f.flush()
            line_count += 1
    f.close()

# print input parameters and experimental affinity to data file from the PPI Affinity dataset for machine learning
def getPPIData():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\ppi_data.txt", 'a')
    with open(cwd + "\\PPI_Dataset\\SI-File-4-protein-protein-test-set-2.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0 and os.path.isfile(cwd + "/PPI_Dataset/pdb" + row[0] + ".ent"):
                calculate(row[0],10,False,"A","B",cwd + "\\Machine_Learning\\ppi_data.txt")
                f.write(str(row[1]) + "\n")
                f.flush()
            line_count += 1
    f.close()

# print input parameters and experimental affinity to data file from the SKEMPI dataset for machine learning
def getSKEMPIData():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\skempi_data.txt", 'a')
    with open(cwd + "\\SKEMPI_Dataset\\skempi.csv") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        pdbs = set()
        for row in csv_reader:
            if line_count != 0 and row[0][0:4] not in pdbs:
                pdbs.add(row[0][0:4])
                calculate(row[0][0:4],10,True,row[0].split("_")[1],row[0].split("_")[2],cwd + "\\Machine_Learning\\skempi_data.txt")
                f.write(str(row[8]) + "\n")
                f.flush()
            line_count += 1
    f.close()

def loopProdigy():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\prodigy_data_2.txt", 'a')
    o = open(cwd + "\\Machine_Learning\\output.txt", 'a')
    for dist in range(1,30):
        o.write(str(dist) + "\t")
        o.flush()
        with open(cwd + "\\PRODIGY_Dataset\\PRODIGY_dataset.csv") as csv_file:
            f.truncate(0)
            csv_reader = csv.reader(csv_file, delimiter=',')
            line_count = 0
            for row in csv_reader:
                if line_count != 0:
                    calculate(row[0][0:4],0,dist,True,"A","B",cwd + "\\Machine_Learning\\prodigy_data_2.txt")
                    f.write(str(row[3]) + "\n")
                    f.flush()
                line_count += 1
            run()
    f.close()
    o.close()

def loopSKEMPI():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\skempi_data.txt", 'a')
    o = open(cwd + "\\Machine_Learning\\output.txt", 'a')
    for dist in range(1,30):
        o.write(str(dist) + "\t")
        o.flush()
        with open(cwd + "\\SKEMPI_Dataset\\skempi.csv") as csv_file:
            f.truncate(0)
            csv_reader = csv.reader(csv_file, delimiter=',')
            line_count = 0
            pdbs = set()
            for row in csv_reader:
                try:
                    if line_count != 0 and row[0][0:4] not in pdbs and float(row[8]) and row[0][0:4] != "1KBH":
                        pdbs.add(row[0][0:4])
                        calculate(row[0][0:4],0,dist,True,row[0].split("_")[1],row[0].split("_")[2],cwd + "\\Machine_Learning\\skempi_data.txt")
                        f.write(row[8])
                        # f.write(str(math.log(row[8])*) + "\n")
                        f.flush()
                except ValueError:
                    continue
                line_count += 1

        run()
    f.close()
    o.close()

def loopPPI():
    cwd = os.getcwd()
    f = open(cwd + "\\Machine_Learning\\ppi_data.txt", 'a')
    o = open(cwd + "\\Machine_Learning\\output.txt", 'a')
    for dist in range(1,30):
        o.write(str(dist) + "\t")
        o.flush()
        with open(cwd + "\\PPI_Dataset\\set_4.csv") as csv_file:
            f.truncate(0)
            csv_reader = csv.reader(csv_file, delimiter=',')
            line_count = 0
            for row in csv_reader:
                if line_count != 0 and os.path.isfile(cwd + "/PPI_Dataset/pdb" + row[0] + ".ent"):
                    calculate(row[0],0,dist,False,"A","B",cwd + "\\Machine_Learning\\ppi_data.txt")
                    f.write(str(row[1]) + "\n")
                    f.flush()
                line_count += 1
        run()
    f.close()
    o.close()








'''

cwd = os.getcwd()
count = 0
listSkempi = os.listdir(cwd + "\\SKEMPI_Dataset")
listProdigy = os.listdir(cwd + "\\PRODIGY_Dataset")

f = open(cwd + "\\Machine_Learning\\test1.txt", 'a')

with open(cwd + "\\SKEMPI_Dataset\\skempi.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    pdbs = set()
    for row in csv_reader:
        try:
            if line_count != 0 and row[0][0:4] not in pdbs and float(row[8]) and row[0][0:4] != "1KBH":
                pdbs.add(row[0][0:4])
                if str(row[0][0:4]) + ".pdb" in listProdigy:
                    calculate(row[0][0:4],0,8,True,row[0].split("_")[1],row[0].split("_")[2],cwd + "\\Machine_Learning\\test1.txt")
                    f.write(" " + str(math.log(float(row[8]))*8.314*274/4184))
                    f.write("\n")
                    f.flush()
        except ValueError:
            continue
        line_count += 1
f.close()

'''

cwd = os.getcwd()
listSkempi = os.listdir(cwd + "\\SKEMPI_Dataset")
listSkempi = [file[0:4] for file in listSkempi]
listProdigy = os.listdir(cwd + "\\PRODIGY_Dataset")
listProdigy = [file[0:4] for file in listProdigy]



f = open(cwd + "\\Machine_Learning\\test2.txt", 'a')
with open(cwd + "\\PRODIGY_Dataset\\PRODIGY_dataset.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count != 0 and row[0][0:4] not in listSkempi:
            calculate(row[0][0:4],0,8,True,"A","B",cwd + "\\Machine_Learning\\test2.txt")
            f.write(str(row[3]) + "\n")
            f.flush()
        line_count += 1
f.close()

