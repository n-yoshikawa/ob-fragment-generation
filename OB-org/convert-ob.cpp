#include <iostream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/forcefield.h>

using namespace std;
using namespace OpenBabel;
int main(int argc, char **argv) {
  OBMol mol;
  OBBuilder builder;

  ifstream ifs;
  OBConversion conv;
  conv.SetOutFormat("SDF");

  for (int i=1; i < argc; i++) {
    OBFormat *inFormat = conv.FormatFromExt(argv[i]);
    if(inFormat==NULL || !conv.SetInFormat(inFormat))
    {
      cerr << " Cannot read file format for " << argv[i] << endl;
      continue; // try next file
    }
    ifs.open(argv[i]);
    if (!ifs) {
      cerr << "Cannot read input file: " << argv[i] << endl;
      continue;
    }
    while(ifs.peek() != EOF && ifs.good()) {
      conv.Read(&mol, &ifs);
      builder.Build(mol);
      mol.AddHydrogens();
      conv.Write(&mol, &cout);
    }
  }
}
