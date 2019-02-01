import sys

import matplotlib.pyplot as plt
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


original_spl = Chem.SDMolSupplier('../data/platinum_dataset_2017_01.sdf')
predicted_spl = Chem.SDMolSupplier(sys.argv[1])

for mol1, mol2 in zip(original_spl, predicted_spl):
    key1 = Chem.inchi.MolToInchiKey(mol1)
    key2 = ""
    if mol2 is not None:
        key2 = Chem.inchi.MolToInchiKey(mol2)
    rmsd = 0
    if key1 == key2:
        rmsd = AllChem.GetBestRMS(mol1, mol2)
    print("{},{},{}".format(key1, key2, rmsd))
