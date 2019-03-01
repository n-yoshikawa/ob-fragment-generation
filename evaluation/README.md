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
Bond error:	 0.06508432761371834
Angle error:	 3.1888053802440015
Torsion error:	 46.58745894803941
TFD:		 0.31343946956243823
Stereo correct:	 93.40369393139841
```

### Open Babel (v 2.4.1)
```
$ python stat-result.py eval-ob-org.csv
RMSD:		 1.7502232199721488
Bond error:	 0.05452582170215406
Angle error:	 2.401368359074101
Torsion error:	 48.79599371069382
TFD:		 0.2732829132474846
Stereo correct:	 76.31926121372031
```

### RDKit (ETKDG)
```
$ python stat-result.py eval-rdkit-etkdg.csv
RMSD:		 1.5884792268202588
Bond error:	 0.060366011566018334
Angle error:	 2.8710107163776186
Torsion error:	 43.87063444926612
TFD:		 0.2145041011963098
Stereo correct:	 99.45030782761654
```

### No TOC
```
$ python stat-result.py eval-ob-new-no-toc.csv
RMSD:		 1.8078555726101713
Bond error:	 0.06507245151643559
Angle error:	 3.1866429685223694
Torsion error:	 46.530899885569305
TFD:		 0.31368021149897174
Stereo correct:	 93.27176781002639
```

### No generic ring
```
$ python stat-result.py eval-ob-new-no-ring.csv
RMSD:		 2.0101836274419465
Bond error:	 0.10313705514704626
Angle error:	 4.746831178727824
Torsion error:	 51.534919458039646
TFD:		 0.5078945447993917
Stereo correct:	 89.42392260334213
```

### Match generic rings before fragments
```
$ python stat-result.py  eval-ob-new-before-ring.csv
RMSD:		 1.7816643076264511
Bond error:	 0.0611101473094311
Angle error:	 2.761833971581092
Torsion error:	 46.85058413795302
TFD:		 0.3149685102419512
Stereo correct:	 92.70008795074757
```

### Fast (100 MMFF)
```
$ python stat-result.py eval-ob-new-fast.csv
RMSD:		 1.740267099697708
Bond error:	 0.051583898472284
Angle error:	 3.153518505423817
Torsion error:	 44.9892862520951
TFD:		 0.2844215908182494
Stereo correct:	 93.57959542656113
```

### Med (100 MMFF + conformer search)
```
$ python stat-result.py eval-ob-new-med.csv
RMSD:		 1.6132163578299799
Bond error:	 0.049156142179546614
Angle error:	 2.6598433969779447
Torsion error:	 43.790369954575596
TFD:		 0.24756291761915955
Stereo correct:	 93.97537379067722
```

### Full database (fastest)
```
$ python stat-result.py eval-ob-new-full.csv
RMSD:		 1.8045784042271176
Bond error:	 0.06531359401631358
Angle error:	 3.217401259391581
Torsion error:	 46.22521610816113
TFD:		 0.32660809258983975
Stereo correct:	 93.7994722955145
```

### Full database (fast)
```
$ python stat-result.py eval-ob-new-full-fast.csv
RMSD:		 1.7447315723284271
Bond error:	 0.05150176501449027
Angle error:	 3.146347208141264
Torsion error:	 44.93128305946052
TFD:		 0.2828497679096479
Stereo correct:	 93.66754617414247
```

### Full database (med)
```
$ python stat-result.py eval-ob-new-full-med.csv
RMSD:		 1.6093877055710286
Bond error:	 0.049827109735998784
Angle error:	 2.6938674497642143
Torsion error:	 43.282242773437375
TFD:		 0.2539131267772019
Stereo correct:	 94.26121372031663
```

### OB-new (New) (fastest)
```
RMSD:		 1.7468039868672727
Bond error:	 0.04871744570622707
Angle error:	 2.4859442266212275
Torsion error:	 44.12171262522172
TFD:		 0.265285440179089
Stereo correct:	 93.86543535620054
```
### OB-new (New) (fast)
```
$ python stat-result.py eval-ob-new-latest-fast.csv
RMSD:		 1.7215534805012254
Bond error:	 0.04920784138282149
Angle error:	 2.897018735784826
Torsion error:	 44.060873050734145
TFD:		 0.26265495609694267
Stereo correct:	 93.57959542656113
```
### OB-new (New) (Med)
```
$ python stat-result.py eval-ob-new-latest-med.csv
RMSD:		 1.5983279064336453
Bond error:	 0.04775556258397763
Angle error:	 2.5178799689573834
Torsion error:	 43.0225304074528
TFD:		 0.2313480015451921
Stereo correct:	 93.93139841688655
```

### OB-new (New)(No TOC)
```
$ python stat-result.py eval-ob-new-latest-no-toc.csv
RMSD:		 1.7884164723469813
Bond error:	 0.05617152808732662
Angle error:	 2.7718740231667924
Torsion error:	 45.527698476386284
TFD:		 0.29755567836935526
Stereo correct:	 92.83201407211962
```

### OB-new (New)(Full, fastest)
```
$ python stat-result.py eval-ob-new-7c73fa1.csv 
RMSD:		 1.7472534035787315
Bond error:	 0.04874244953799942
Angle error:	 2.487054065741934
Torsion error:	 44.14656531895596
TFD:		 0.26548377011134
Stereo correct:	 93.7994722955145
```

### OB-new (New)(Full, fast)
```
$ python stat-result.py eval-7c73fa1-fast.csv 
RMSD:		 1.7035357156797966
Bond error:	 0.04880397234219525
Angle error:	 2.7549237686503845
Torsion error:	 43.04960099653375
TFD:		 0.24407059862995575
Stereo correct:	 94.08531222515391
```

### OB-new (New)(Full, med)
```
$ python stat-result.py eval-ob-new-7c73fa1-med.csv 
RMSD:		 1.5780111350238781
Bond error:	 0.047667297385034835
Angle error:	 2.4779011947305953
Torsion error:	 42.28814134948697
TFD:		 0.2231891355802418
Stereo correct:	 94.1952506596306
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
