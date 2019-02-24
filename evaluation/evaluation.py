import sys

import numpy as np

import pybel
ob = pybel.ob
etab = ob.OBElementTable()

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

isDebug = False

refFileName = '../data/platinum_dataset_2017_01.sdf'
predFileName = sys.argv[1]

# Calculate RMSD with RDKit
refSpl = Chem.SDMolSupplier(refFileName)
predSpl = Chem.SDMolSupplier(predFileName)

entry2RMSD = {}
for ref, pred in zip(refSpl, predSpl):
    refEntry = ref.GetProp('_Name')
    predEntry = pred.GetProp('_Name')
    assert(refEntry == predEntry)
    try:
        rmsd = AllChem.GetBestRMS(ref, pred)
    except:
        rmsd = ''
    entry2RMSD[refEntry] = rmsd

# See https://baoilleach.blogspot.com/2010/11/automorphisms-isomorphisms-symmetry.html
print("Entry,SMILES,RMSD,Bond error,Angle error,Torsion error,Stereo correct")
for ref, pred in zip(pybel.readfile("sdf", refFileName),
                     pybel.readfile("sdf", predFileName)):
    refMol = ref.OBMol
    predMol = pred.OBMol

    refEntry = refMol.GetTitle()
    predEntry = predMol.GetTitle()
    assert(refEntry == predEntry)

    refSMILES = ref.write("can").split()[0]

    # Check stereochemistry
    refKey = ref.write("inchikey").rstrip()
    predKey = pred.write("inchikey").rstrip()
    if refKey != predKey:  # Wrong stereochemistry
        print("{},{},,,,,F".format(refEntry, refSMILES))
        continue

    # Get isomorphism between two molecules
    query = ob.CompileMoleculeQuery(refMol)
    mapper = ob.OBIsomorphismMapper.GetInstance(query)
    isomorphs = ob.vvpairUIntUInt()
    # Check validity (OB issue #1929)
    mapper.MapAll(pred.OBMol, isomorphs)
    for isomorph in isomorphs:
        idxMap = dict(isomorph)
        isCorrectMap = True
        for refAtom in ob.OBMolAtomIter(refMol):
            refIdx = refAtom.GetIdx()
            # obmol indices are 1-indexed while the mapper is zero indexed
            predIdx = idxMap[refIdx - 1] + 1
            if refAtom.GetAtomicNum() != predMol.GetAtom(predIdx).GetAtomicNum():
                isCorrectMap = False
                break
        if isCorrectMap is True:
            break
    # No correct isomorphism is found
    if isCorrectMap is False:
        print("{},{},{},,,,T".format(
            refEntry, refSMILES, entry2RMSD[refEntry]))
        continue

    # Calculate bond length error
    bondLenErrors = []
    for refBond in ob.OBMolBondIter(refMol):
        refBeginIdx = refBond.GetBeginAtomIdx()
        refEndIdx = refBond.GetEndAtomIdx()
        predBeginIdx = idxMap[refBeginIdx - 1] + 1
        predEndIdx = idxMap[refEndIdx - 1] + 1
        predBond = predMol.GetBond(predBeginIdx, predEndIdx)
        bondLenErrors.append(
            abs(refBond.GetLength()-predBond.GetLength())/refBond.GetLength())
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
        predBondAngle = [idxMap[t] + 1 for t in refBondAngle]
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
        predTorsion = [idxMap[t] + 1 for t in refTorsion]
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

    print("{},{},{},{},{},{},T".format(refEntry, refSMILES,
                                       entry2RMSD[refEntry],
                                       np.mean(bondLenErrors)*100.0,
                                       np.mean(bondAngleErrors),
                                       np.mean(torsionErrors)))
