# Prediction of Protein-Protein Binding Affinity Using Amino Acid Counts at the Interface

This repository contains to code to make a Delta G prediction for a protein-protein complex given a PDB file. It calculates interface numbers for six amino acids and makes a prediction with a linear regression equation that was found by training models on a dataset of 141 complexes.

Paper: https://pubs.acs.org/doi/10.1021/acsomega.3c06996

## Functions

### predict(name, chains=None)

Returns a prediction given the name of a complex as a four-character string (not case-sensitive) referring to its PDB ID. The complex must already exist in the pdbFiles directory. Chains argument is optional and should be a tuple of length two specifying which chains belong to the first and second protein. Otherwise, the first two chains of the complex are used. This is essentially a wrapper for calculating the interface numbers and putting that into the regression equation.

```python
predict("1A2K")
predict("1AHW",("AB","C"))
```

### getPDB(name)

Downloads the PDB file given its four-character ID (not case-sensitive) into the pdbFiles directory. This uses the Biopython module which accesses the PDB database. Some complexes may not be found and must be added manually to the pdbFile directory. The resulting name of the file in the pdbFiles directory is of the form pdb + PDB ID + .ent.

```python
getPDB("1AHW")
```

This creates a file called pdb1ahw.ent in the pdbFiles directory.

### getAAINS(pdbFile, chains)

Given the ID of a pdbfile that exists in the pdbFiles directory, returns a dictionary of length 6 containing the interface numbers of the six amino acids used in the model. If the chains argument is None, uses the first two chains in the complex to represent the two proteins. Otherwise chains is a tuple containing the names of the chains in each protein.

```python
getAAINs("1A2K")
getAAINs("1AHW",("AB","C"))
```

Note that this takes in the same arguments as the predict(name, chains=None) function and in both cases the pdbfile must already exist in the pdbFiles directory.

### predictFromINs(INs)

Given a dictionary with the interface numbers of the six amino acids, makes a prediction.

```python
predictFromINs(getAAINs("1AHW",("AB","C")))
```
