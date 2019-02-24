#include <iostream>

#include <openbabel/builder.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;
int main(int argc, char **argv) {
  OBMol mol;
  OBBuilder builder;

  ifstream ifs;
  OBConversion conv;
  conv.SetOutFormat("SDF");

  for (int i = 1; i < argc; i++) {
    OBFormat *inFormat = conv.FormatFromExt(argv[i]);
    if (inFormat == NULL || !conv.SetInFormat(inFormat)) {
      cerr << " Cannot read file format for " << argv[i] << endl;
      continue; // try next file
    }
    ifs.open(argv[i]);
    if (!ifs) {
      cerr << "Cannot read input file: " << argv[i] << endl;
      continue;
    }
    while (ifs.peek() != EOF && ifs.good()) {
      conv.Read(&mol, &ifs);
      builder.Build(mol);
      mol.AddHydrogens();

      // Cleanup by MMFF
      OBForceField *pFF = OBForceField::FindForceField("MMFF94");
      if (!pFF) continue;
      if (!pFF->Setup(mol)) {
        pFF = OBForceField::FindForceField("UFF");
        if (!pFF || !pFF->Setup(mol))
          continue; // can't use either MMFF94 or UFF
      }

      // Since we only want a rough geometry, use distance cutoffs for VDW,
      // Electrostatics
      pFF->EnableCutOff(true);
      pFF->SetVDWCutOff(10.0);
      pFF->SetElectrostaticCutOff(20.0);
      pFF->SetUpdateFrequency(10); // update non-bonded distances infrequently

      // uncomment here for evaluating fast MMFF
      /*
      // How many cleanup cycles?
      int iterations = 100;
      // Initial cleanup for every level
      pFF->ConjugateGradients(iterations, 1.0e-4);
      pFF->UpdateCoordinates(mol);
      */

      // uncomment here for evaluating medium MMFF
      /*
      // How many cleanup cycles?
      int iterations = 100;
      // Initial cleanup for every level
      pFF->ConjugateGradients(iterations, 1.0e-4);
      // permute central rotors
      pFF->FastRotorSearch(false);
      // Final cleanup and copy the new coordinates back
      pFF->ConjugateGradients(iterations, 1.0e-6);
      pFF->UpdateCoordinates(mol);
      */

      conv.Write(&mol, &cout);
    }
  }
}
