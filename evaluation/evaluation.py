import sys

import numpy as np

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, TorsionFingerprints

import pybel
ob = pybel.ob
etab = ob.OBElementTable()


isDebug = False

refFileName = '../data/platinum_dataset_2017_01.sdf'
predFileName = sys.argv[1]

# Calculate RMSD with RDKit
refSpl = Chem.SDMolSupplier(refFileName)
predSpl = Chem.SDMolSupplier(predFileName)

entry2RMSD = {}
entry2TFD = {}
for ref, pred in zip(refSpl, predSpl):
    refEntry = ref.GetProp('_Name')
    if pred is None:  # in case of failure
        entry2RMSD[refEntry] = ''
        entry2TFD[refEntry] = ''
        continue

    predEntry = pred.GetProp('_Name')
    assert(refEntry == predEntry)
    try:
        rmsd = AllChem.GetBestRMS(ref, pred)
    except:
        rmsd = ''
    try:
        m = Chem.MolFromSmiles(Chem.MolToSmiles(ref))
        ref = AllChem.AssignBondOrdersFromTemplate(m, ref)
        pred = AllChem.AssignBondOrdersFromTemplate(m, pred)
        tfd = TorsionFingerprints.GetTFDBetweenMolecules(ref, pred)
    except:
        tfd = ''
    entry2RMSD[refEntry] = rmsd
    entry2TFD[refEntry] = tfd

# See https://baoilleach.blogspot.com/2010/11/automorphisms-isomorphisms-symmetry.html
print("Entry,SMILES,RMSD,Bond error,Angle error,Torsion error,TFD,Stereo correct")
for ref, pred in zip(pybel.readfile("sdf", refFileName),
                     pybel.readfile("sdf", predFileName)):
    refMol = ref.OBMol
    predMol = pred.OBMol

    refEntry = refMol.GetTitle()
    predEntry = predMol.GetTitle()
    assert refEntry == predEntry

    refSMILES = ref.write("can").split()[0]
    # Check stereochemistry
    refKey = ref.write("inchikey").rstrip()
    predKey = pred.write("inchikey").rstrip()
    if refKey != predKey:  # Wrong stereochemistry
        print("{},{},,,,,,F".format(refEntry, refSMILES))
        continue

    # Get mapping of two molecules using canonical SMILES order
    # https://github.com/openbabel/openbabel/pull/1712
    gs = ob.OBGraphSym(refMol)
    symmetry = ob.vectorUnsignedInt()
    gs.GetSymmetry(symmetry)
    canonical_labels = ob.vectorUnsignedInt()
    ob.CanonicalLabels(refMol, symmetry, canonical_labels)
    refOrder = list(canonical_labels)
    gs = ob.OBGraphSym(predMol)
    symmetry = ob.vectorUnsignedInt()
    gs.GetSymmetry(symmetry)
    canonical_labels = ob.vectorUnsignedInt()
    ob.CanonicalLabels(predMol, symmetry, canonical_labels)
    predOrder = list(canonical_labels)
    predOrderInv = {order: idx for idx, order in enumerate(predOrder)}
    idxMap = {idx+1: predOrderInv[order] + 1
              for idx, order in enumerate(refOrder)}

    for refAtom in ob.OBMolAtomIter(refMol):
        refIdx = refAtom.GetIdx()
        predIdx = idxMap[refIdx]
        assert refAtom.GetAtomicNum() == predMol.GetAtom(predIdx).GetAtomicNum()

    # Calculate bond length error
    bondLenErrors = []
    for refBond in ob.OBMolBondIter(refMol):
        refBeginIdx = refBond.GetBeginAtomIdx()
        refEndIdx = refBond.GetEndAtomIdx()
        predBeginIdx = idxMap[refBeginIdx]
        predEndIdx = idxMap[refEndIdx]
        predBond = predMol.GetBond(predBeginIdx, predEndIdx)
        bondLenErrors.append(
            abs(refBond.GetLength()-predBond.GetLength()))
        if isDebug:
            print(etab.GetSymbol(refMol.GetAtom(refBeginIdx).GetAtomicNum()),
                  etab.GetSymbol(refMol.GetAtom(refEndIdx).GetAtomicNum()),
                  refBond.GetLength())
            print(etab.GetSymbol(predMol.GetAtom(predBeginIdx).GetAtomicNum()),
                  etab.GetSymbol(predMol.GetAtom(predEndIdx).GetAtomicNum()),
                  predBond.GetLength())
            print()

    # Calculate bond angle error
    bondAngleErrors = []
    for refBondAngle in ob.OBMolAngleIter(refMol):
        predBondAngle = [idxMap[t+1] for t in refBondAngle]
        refBondAngle = [t + 1 for t in refBondAngle]
        refBondAngleDeg = refMol.GetAtom(refBondAngle[1]).GetAngle(
            refBondAngle[0], refBondAngle[2])
        predBondAngleDeg = predMol.GetAtom(predBondAngle[1]).GetAngle(
            predBondAngle[0], predBondAngle[2])
        bondAngleErrors.append(abs(refBondAngleDeg - predBondAngleDeg))
        if isDebug:
            print([etab.GetSymbol(refMol.GetAtom(i).GetAtomicNum()) for i in refBondAngle],
                  refBondAngleDeg)
            print([etab.GetSymbol(predMol.GetAtom(i).GetAtomicNum()) for i in predBondAngle],
                  predBondAngleDeg)
            print()

    # Calculate torsion angle error
    torsionErrors = []
    for refTorsion in ob.OBMolTorsionIter(refMol):
        predTorsion = [idxMap[t+1] for t in refTorsion]
        refTorsion = [t + 1 for t in refTorsion]
        refTorsionAngleDeg = refMol.GetTorsion(
            refTorsion[0], refTorsion[1], refTorsion[2], refTorsion[3])
        predTorsionAngleDeg = predMol.GetTorsion(
            predTorsion[0], predTorsion[1], predTorsion[2], predTorsion[3])
        angleError = abs(refTorsionAngleDeg - predTorsionAngleDeg)
        torsionErrors.append(min(angleError, 360.0-angleError))
        if isDebug:
            print([etab.GetSymbol(refMol.GetAtom(i).GetAtomicNum()) for i in refTorsion],
                  refTorsionAngleDeg)
            print([etab.GetSymbol(predMol.GetAtom(i).GetAtomicNum()) for i in predTorsion],
                  predTorsionAngleDeg)
            print()

    print("{},{},{},{},{},{},{},T".format(refEntry, refSMILES,
                                       entry2RMSD[refEntry],
                                       np.mean(bondLenErrors),
                                       np.mean(bondAngleErrors),
                                       np.mean(torsionErrors),
                                       entry2TFD[refEntry]))
