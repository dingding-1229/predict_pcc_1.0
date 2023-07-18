## raw folder:
Contains raw chemical data of pharmaceutical patents downloaded from the CAS database, including training/valid1/valid2 sets.

## source folder:
Contains normalized data processed from raw data, including train/valid1/valid2 sets.

## similarity folder:
Contains similarity data for permuted pairs of molecules, including train/effective1/effective2 sets.

## connected folder:
Contains features and network graphs for drugs, including training/valid1/valid2 sets.

## calc_draw.py:
Python code for building networks, calculating parameters and drawing graphs.

## sim_calc.py:
Python code for computing similarity values for pairs of molecules in a patent by permutation.

## model.ipynb:
Python code for MLP & Xgboost.

## parallel.ipynb:
Jupyter notebook file responsible for parallel computation.

## preprocess.ipynb:
Jupyter notebook file responsible for preprocessing raw data.

## train.xlsx:
Excel worksheet containing the essential data (e.g., drug names, patents number, SMILES) of training set.

## valid1.xlsx:
Excel worksheet containing the essential data (e.g., drug names, patents number, SMILES) of valid1 set.

## valid2.xlsx:
excel worksheet containing the essential data (e.g., drug names, patents number, SMILES) of valid2 set.

## future.txt:
contain ideas of future works.
