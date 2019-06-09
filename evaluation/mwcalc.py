import sys

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

refFileName = '../data/platinum_dataset_2017_01.sdf'
predFileName = sys.argv[1]

# Calculate RMSD with RDKit
refSpl = Chem.SDMolSupplier(refFileName)
predSpl = Chem.SDMolSupplier(predFileName)

mwList = []
rmsdList = []
print("MW,RMSD")
for ref, pred in zip(refSpl, predSpl):
    refEntry = ref.GetProp('_Name')
    if pred is None:  # in case of failure
        continue

    predEntry = pred.GetProp('_Name')
    assert(refEntry == predEntry)

    mw = Chem.Descriptors.MolWt(pred)
    try:
        rmsd = AllChem.GetBestRMS(ref, pred)
        print("{},{}".format(mw, rmsd))
    except:
        pass
