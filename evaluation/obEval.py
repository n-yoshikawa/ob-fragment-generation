import sys

import matplotlib.pyplot as plt
import numpy as np

import pybel
ob = pybel.ob

org_filename = '../data/platinum_dataset_2017_01.sdf'
pred_filename = sys.argv[1]

# See https://baoilleach.blogspot.com/2010/11/automorphisms-isomorphisms-symmetry.html
for ref, pred in zip(pybel.readfile("sdf", org_filename), pybel.readfile("sdf", pred_filename)):
    refMol = ref.OBMol
    predMol = pred.OBMol
    refKey = ref.write("inchikey").rstrip()
    predKey = pred.write("inchikey").rstrip()
    if refKey != predKey:
        print("{},{},null,null,null".format(refKey, predKey))
        continue
    query = ob.CompileMoleculeQuery(refMol)
    mapper = ob.OBIsomorphismMapper.GetInstance(query)
    isomorphs = ob.vvpairUIntUInt()
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
        if isCorrectMap == True:
            break
    if isCorrectMap == False:
        print("{},{},null,null,null".format(refKey, predKey))
        continue

    # Calculate bond length error
    bondLenErrors = []
    for refBond in ob.OBMolBondIter(refMol):
        refBeginIdx = refBond.GetBeginAtomIdx()
        refEndIdx = refBond.GetEndAtomIdx()
        predBeginIdx = idxMap[refBeginIdx - 1] + 1
        predEndIdx = idxMap[refEndIdx - 1] + 1
        predBond = predMol.GetBond(predBeginIdx, predEndIdx)
        bondLenErrors.append(abs(refBond.GetLength() - predBond.GetLength()))

    # Calculate bond angle error
    bondAngleErrors = []
    for refBondAngle in ob.OBMolAngleIter(refMol):
        predBondAngle = [idxMap[t] + 1 for t in refBondAngle]
        refBondAngle = [t + 1 for t in refBondAngle]
        refBondAngle = refMol.GetAtom(refBondAngle[1]).GetAngle(refBondAngle[0], refBondAngle[2])
        predBondAngle = predMol.GetAtom(predBondAngle[1]).GetAngle(predBondAngle[0], predBondAngle[2])
        bondAngleErrors.append(abs(refBondAngle - predBondAngle))

    # Calculate torsion angle error
    torsionErrors = []
    for refTorsion in ob.OBMolTorsionIter(refMol):
        predTorsion = [idxMap[t] + 1 for t in refTorsion]
        refTorsion = [t + 1 for t in refTorsion]
        refTorsionAngle = refMol.GetTorsion(refTorsion[0], refTorsion[1], refTorsion[2], refTorsion[3])
        predTorsionAngle = predMol.GetTorsion(predTorsion[0], predTorsion[1], predTorsion[2], predTorsion[3])
        torsionErrors.append(abs(refTorsionAngle - predTorsionAngle))
    print("{},{},{},{},{}".format(refKey, predKey, 
          np.mean(bondLenErrors), np.mean(bondAngleErrors), np.mean(torsionErrors)))
