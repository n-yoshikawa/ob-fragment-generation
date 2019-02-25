# evaluation
This directory includes evaluation scripts and results

## Contents
- `evaluation.py`: Evaluation script. 
- `eval-*.csv`: Result of evaluation for SDFs
- `stat-result.py`: Script for getting statistics of results

You can run evaluation.py like this:
```
python evaluation.py <Prediction SDF Name>
```
RDKit and Python binding of Open Babel are required.

## Results
### This work
```
$ python stat-result.py eval-ob-new.csv 
RMSD:		 1.8063509451871946
Bond error:	 0.06509318496650597
Angle error:	 3.1895608307951284
Torsion error:	 45.93642429997671
Stereo correct:	 93.40369393139841
```

### Open Babel (v 2.4.1)
```
$ python stat-result.py eval-ob-org.csv
RMSD:		 1.7502232199721488
Bond error:	 0.05452480256123688
Angle error:	 2.400959138756347
Torsion error:	 48.23951761023113
Stereo correct:	 76.31926121372031
```

### RDKit (ETKDG)
```
$ python stat-result.py eval-rdkit-etkdg.csv
RMSD:		 1.5884792268202588
Bond error:	 0.060366936881424395
Angle error:	 2.8715841660752797
Torsion error:	 43.87665476075526
Stereo correct:	 99.45030782761654
```

### No TOC
```
$ python stat-result.py eval-ob-new-no-toc.csv
RMSD:		 1.8078555726101713
Bond error:	 0.06508412014635938
Angle error:	 3.186703011580234
Torsion error:	 45.89391180731905
Stereo correct:	 93.27176781002639
```

### No generic ring
```
$ python stat-result.py eval-ob-new-no-ring.csv
RMSD:		 2.0101836274419465
Bond error:	 0.1031472553026142
Angle error:	 4.747287183181868
Torsion error:	 51.15567624620266
Stereo correct:	 89.42392260334213
```

### Match generic rings before fragments
```
$ python stat-result.py  eval-ob-new-before-ring.csv
RMSD:		 1.7816643076264511
Bond error:	 0.06110149529280588
Angle error:	 2.7613183669483905
Torsion error:	 46.21791284536982
Stereo correct:	 92.70008795074757
```

### Fast (100 MMFF)
```
$ python stat-result.py eval-ob-new-fast.csv 
RMSD:		 1.740267099697708
Bond error:	 0.051581857402544194
Angle error:	 3.1540063845017365
Torsion error:	 44.190972013030006
Stereo correct:	 93.57959542656113
```

### Med (100 MMFF + conformer search)
```
$ python stat-result.py eval-ob-new-med.csv 
RMSD:		 1.6132163578299799
Bond error:	 0.049155235540285225
Angle error:	 2.6594449843860346
Torsion error:	 42.84335814050267
Stereo correct:	 93.97537379067722
```

### Full database (fastest)
```
$ python stat-result.py eval-ob-new-full.csv 
RMSD:		 1.8045784042271176
Bond error:	 0.0653110533921927
Angle error:	 3.2169574296681227
Torsion error:	 45.33813672863018
Stereo correct:	 93.7994722955145
```

### Full database (fast)
```
$ python stat-result.py eval-ob-new-full-fast.csv 
RMSD:		 1.7447315723284271
Bond error:	 0.051503947160110364
Angle error:	 3.1468707366577817
Torsion error:	 44.14667804445053
Stereo correct:	 93.66754617414247
```

### Full database (med)
```
$ python stat-result.py eval-ob-new-full-med.csv 
RMSD:		 1.6093877055710286
Bond error:	 0.04982789335805958
Angle error:	 2.6949024308130722
Torsion error:	 42.38347327588605
Stereo correct:	 94.26121372031663
```

### Evaluation
```
python evaluation.py ../OB-new/platinum-ob-new.sdf > eval-ob-new.csv 2>/dev/null
python evaluation.py ../OB-org/platinum-ob-org.sdf > eval-ob-org.csv 2>/dev/null
python evaluation.py ../rdkit/platinum-rdkit-etkdg.sdf > eval-rdkit-etkdg.csv 2>/dev/null
python evaluation.py ../OB-new/platinum-ob-new-no-toc.sdf > eval-ob-new-no-toc.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-no-ring.sdf > eval-ob-new-no-ring.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-before-ring.sdf > eval-ob-new-before-ring.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-fast.sdf > eval-ob-new-fast.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-med.sdf > eval-ob-new-med.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-full.sdf > eval-ob-new-full.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-full-fast.sdf > eval-ob-new-full-fast.csv 2> /dev/null
python evaluation.py ../OB-new/platinum-ob-new-full-med.sdf > eval-ob-new-full-med.csv 2> /dev/null
```
