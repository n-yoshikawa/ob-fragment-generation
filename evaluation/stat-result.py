import argparse
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1])

stereoCorrect = df['Stereo correct'] == 'T'
print("RMSD:\t\t", df.loc[stereoCorrect, 'RMSD'].mean())
print("Bond error:\t", df.loc[stereoCorrect, 'Bond error'].mean())
print("Angle error:\t", df.loc[stereoCorrect, 'Angle error'].mean())
print("Torsion error:\t", df.loc[stereoCorrect, 'Torsion error'].mean())
print("Stereo correct:\t", stereoCorrect.sum() / float(len(df)) * 100.0)

#print()
#print("good mol: {}%".format(((df['Stereo correct'] == 'T')
#                              & (df['RMSD'].astype(float) < 0.5)).sum() / float(len(df)) * 100.0))
#print("bad mol: {}%".format(((df['Stereo correct'] == 'F')
#                             | (df['RMSD'].astype(float) > 3.5)).sum() / float(len(df)) * 100.0))

plotRMSD = False

if plotRMSD is True:
    plt.hist(df['RMSD'].dropna(), bins=20)
    plt.xlabel('RMSD', fontsize=24)
    plt.ylabel('Number of molecules', fontsize=24)
    plt.tick_params(labelsize=18)
    plt.tight_layout()
    plt.show()
