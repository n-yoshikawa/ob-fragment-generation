from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

w = Chem.SDWriter('platinum-rdkit-etkdg.sdf')

with open("../data/platinum.smi", "r") as f:
    for line in f:
        smiles, entry = line.split()
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", entry)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        w.write(mol)
w.flush()
