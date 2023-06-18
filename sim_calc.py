import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import multiprocessing
import os


sim_path='./similarity'
source_path='./source/train'

def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius = 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius = 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    return s

def combinations_generator(data):
    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            yield data[i], data[j]

def compute_sim(i):
    newname=str('sim.')+i
    if os.path.exists(sim_path+'/'+newname):
        print(newname,' already exists!')
    else:
        sim = pd.DataFrame(columns=['node1','node2','similarity'])
        drug=pd.read_csv(source_path+'/'+i)
        smiles_list = drug['Canonical SMILES'].tolist()
        cas_list = drug['CAS Registry Number'].tolist()
        smiles_combinations = combinations_generator(smiles_list)
        cas_combinations = combinations_generator(cas_list)
        sim_dict={}
        for ii, jj in zip(smiles_combinations, cas_combinations):
            sim_dict[jj] = tanimoto_calc(ii[0], ii[1])
        sim['similarity']=list(sim_dict.values())
        sim['node1'] = [x[0] for x in sim_dict.keys()]
        sim['node2'] = [x[1] for x in sim_dict.keys()]
        sim.to_csv(sim_path+'/'+newname,sep=',',index=False)
        print(newname,' calculated!')
 
        
def main():
    file_list = os.listdir(source_path)
    pool = multiprocessing.Pool()
    pool.map(compute_sim, file_list)
    pool.close()
    pool.join()
    