import sys

import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1])

stereoCorrect = df['Stereo correct'] == 'T'
print("RMSD:\t\t", df.loc[stereoCorrect, 'RMSD'].mean())
print("Bond error:\t", df.loc[stereoCorrect, 'Bond error'].mean())
print("Angle error:\t", df.loc[stereoCorrect, 'Angle error'].mean())
print("Torsion error:\t", df.loc[stereoCorrect, 'Torsion error'].mean())
print("TFD:\t\t", df.loc[stereoCorrect, 'TFD'].mean())
print("Stereo correct:\t", stereoCorrect.sum() / float(len(df)) * 100.0)

# print()
# print("good mol: {}%".format(((df['Stereo correct'] == 'T')
#                              & (df['RMSD'].astype(float) < 0.5)).sum() / float(len(df)) * 100.0))
# print("good mol (1.0): {}%".format(((df['Stereo correct'] == 'T')
#                              & (df['RMSD'].astype(float) < 1.0)).sum() / float(len(df)) * 100.0))
# print("bad mol: {}%".format(((df['Stereo correct'] == 'F')
#                             | (df['RMSD'].astype(float) > 3.5)).sum() / float(len(df)) * 100.0))
#
