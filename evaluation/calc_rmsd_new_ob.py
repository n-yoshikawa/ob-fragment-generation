import matplotlib.pyplot as plt
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

import traceback


original_spl = Chem.SDMolSupplier('../data/platinum_dataset_2017_01.sdf')
predicted_spl = Chem.SDMolSupplier('../OB-new/platinum-ob-new.sdf')

original = []
for mol in original_spl:
    original.append(mol)

predicted = []
for mol in predicted_spl:
    predicted.append(mol)

rmsd = []
for mol1, mol2 in zip(original, predicted):
    if mol2 is None or Chem.MolToSmiles(mol1) != Chem.MolToSmiles(mol2):
        continue
    rmsd.append(AllChem.GetBestRMS(mol1, mol2))


print("No of success:", len(rmsd))
print("Mean RMSD:", np.mean(rmsd))

plt.title("Distribution of RMSD (New Open Babel)")
plt.xlabel("RMSD")
plt.ylabel("Number of molecules")
plt.hist(rmsd, bins=20)
plt.show()

