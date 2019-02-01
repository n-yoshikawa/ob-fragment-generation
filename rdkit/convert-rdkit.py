from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

w = Chem.SDWriter('platinum-rdkit-etkdg.sdf')

with open("../data/platinum.smi", "r") as f:
    for line in f:
        smiles = line.split()[0]
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        w.write(mol)
w.flush()
