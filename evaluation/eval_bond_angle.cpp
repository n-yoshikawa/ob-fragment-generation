#include <openbabel/isomorphism.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/query.h>
#include <openbabel/math/align.h>

#include <cassert>
#include <cmath>
#include <iterator>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>

using namespace std;
using namespace OpenBabel;

int main(int argc, char **argv) {
  OBConversion conv;
  //conv.SetOptions("O", conv.OUTOPTIONS);
  conv.SetInFormat("sdf");
  conv.SetOutFormat("inchikey");

  OBFormat *inFormat;
  if(argc != 3) {
    cerr << argv[0] << "<original SDF>" << "<predicted SDF>" << endl;
  }

  ifstream ifs1, ifs2;
  ifs1.open(argv[1]);
  ifs2.open(argv[2]);

  while(ifs1.peek() != EOF && ifs1.good() &&
        ifs2.peek() != EOF && ifs2.peek()) {
    // mol1: ground truth molecule
    // mol2: predicted molecule
    OBMol mol1, mol2;

    stringstream ss1, ss2;
    string key1, key2;

    conv.Read(&mol1, &ifs1);
    conv.Read(&mol2, &ifs2);

    // Make sure InChIKey of two molecules are the same
    conv.Write(&mol1, &ss1);
    conv.Write(&mol2, &ss2);

    ss1 >> key1;
    ss2 >> key2;

    if(key1 != key2) {
      cout << key1 << "," << key2 << ",0,0,0" << endl;
      continue;
    }

    // Calculate RMSD
    OBAlign aln(mol1, mol2);
    aln.Align();
    double rmsd = aln.GetRMSD();

    // Get mapping between two molecules
    OBQuery *query = CompileMoleculeQuery(&mol1);
    OBIsomorphismMapper *mapper = OBIsomorphismMapper::GetInstance(query);
    OBIsomorphismMapper::Mappings mappings;
    mapper->MapAll(&mol2, mappings);

    map<int, int> idxMapping; // mol1 index to mol2 index

    for(size_t m = 0; m < mappings.size(); ++m) {
      bool isCorrect = true;
      for(size_t i = 0; i < mappings[m].size(); ++i) {
        //obmol indices are 1-indexed while the mapper is zero indexed
        unsigned int idx1 = mappings[m][i].first + 1;
        unsigned int idx2 = mappings[m][i].second + 1;
        idxMapping[idx1] = idx2;

        if(mol1.GetAtom(idx1)->GetAtomicNum() != mol2.GetAtom(idx2)->GetAtomicNum()) {
          idxMapping.clear();
          isCorrect = false;
          break;
        }
      }
      if(isCorrect) break;
    }
    if(idxMapping.empty()) {
      cout << key1 << "," << key2 << "," << rmsd << ",null,null" << endl;
      continue;
    }


    // Display the length of each bond
    double bondError = 0;
    int bondNum = 0;
    FOR_BONDS_OF_MOL(bond1, mol1) {
      unsigned int begin1 = idxMapping[bond1->GetBeginAtomIdx()];
      unsigned int end1 = idxMapping[bond1->GetEndAtomIdx()];
      bondNum++;

      bool isFound = false;
      FOR_BONDS_OF_MOL(bond2, mol2) {
        unsigned int begin2 = bond2->GetBeginAtomIdx();
        unsigned int end2 = bond2->GetEndAtomIdx();
        if((begin1 == begin2 && end1 == end2) || (begin1 == end2 && end1 == begin2)) {
          bondError += abs(bond1->GetLength() - bond2->GetLength()) / bond1->GetLength() ;
          isFound = true;
          break;
        }
      }
      if(!isFound) {
        cout << key1 << "," << key2 << "," << rmsd << ",null,null" << endl;
        continue;
      }
    }
    bondError /= bondNum; 

    double angleError = 0;
    int angleNum = 0;
    FOR_ANGLES_OF_MOL(angle1, mol1) {
      unsigned int idxA, idxB, idxC;
      OBAtom *a1, *b1, *c1;
      double angleDeg1;
      b1 = mol1.GetAtom((*angle1)[0] + 1);
      a1 = mol1.GetAtom((*angle1)[1] + 1);
      c1 = mol1.GetAtom((*angle1)[2] + 1);
      idxA = idxMapping[a1->GetIdx()];
      idxB = idxMapping[b1->GetIdx()];
      idxC = idxMapping[c1->GetIdx()];
      angleDeg1 = a1->GetAngle(b1->GetIdx(), c1->GetIdx());
      angleNum++;
      bool isFound = false;
      FOR_ANGLES_OF_MOL(angle2, mol2) {
        OBAtom *a2, *b2, *c2;
        b2 = mol2.GetAtom((*angle2)[0] + 1);
        a2 = mol2.GetAtom((*angle2)[1] + 1);
        c2 = mol2.GetAtom((*angle2)[2] + 1);
        if(idxB == b2->GetIdx() && ((idxA == a2->GetIdx() && idxC == c2->GetIdx()) ||
                                    (idxA == c2->GetIdx() && idxC == a2->GetIdx()))) {
          double angleDeg2 = a2->GetAngle(b2->GetIdx(), c2->GetIdx());
          angleError += abs(angleDeg1 - angleDeg2) / angleDeg1;
          isFound = true;
          break;
        }
      }
      assert(isFound);
    }
    angleError /= angleNum;

    cout << key1 << "," << key2 << "," << rmsd << "," << setprecision(10) << bondError << "," << angleError << endl;
  }
}
