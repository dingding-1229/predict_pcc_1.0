## raw folder:
contain the raw chemical data of drug patents downloaded from CAS database, including the training/valid1/valid2 sets.

## source folder:
contain the formalized data preprocessed from raw data, including the training/valid1/valid2 sets.

## similarity folder:
contain the similarity data of molecule pairs by permutation, including the training/valid1/valid2 sets.

## connected folder:
contain the features and network graphs of drugs, including the training/valid1/valid2 sets.

## calc_draw.py:
python code for building networks, calculating parameters and drawing graphs.

## sim_calc.py:
python code for calculating similarity values of molecule pairs in a patent by permutation.

## parallel.ipynb:
jupyter notebook file responsible for parallel computation.

## preprocess.ipynb:
jupyter notebook file responsible for preprocessing raw data.

## train.xlsx:
excel worksheet containing the essential data (e.g., drug names, patents number, SMILES) of training set.

## valid1.xlsx:
excel worksheet containing the essential data (e.g., drug names, patents number, SMILES) of valid1 set.

## valid2.xlsx:
excel worksheet containing the essential data (e.g., drug names, patents number, SMILES) of valid2 set.

## future.txt:
contain ideas of future works.