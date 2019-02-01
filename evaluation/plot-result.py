import argparse
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np

rmsd_list_rdkit = []
rmsd_list_ob = []
bond_list = []
angle_list = []

with open(sys.argv[1], 'r') as f1:
    with open(sys.argv[2], 'r') as f2:
        reader1 = csv.reader(f1)
        reader2 = csv.reader(f2)
        for row1, row2 in zip(reader1, reader2):
            key1_rdkit, key2_rdkit, rmsd_rdkit = row1
            key1_ob, key2_ob, rmsd_ob, bond, angle = row2
            if key1_rdkit != key1_ob:
                print(key1_rdkit, key1_ob)
            if key1_rdkit == key2_rdkit:
                rmsd_list_rdkit.append(float(rmsd_rdkit))
            if key1_ob == key2_ob:
                rmsd_list_ob.append(float(rmsd_ob))
            if bond != 'null':
                bond_list.append(float(bond))
            if angle != 'null':
                angle_list.append(float(angle))

print("No of success (RDKit):", len(rmsd_list_rdkit))
print("No of success (Open Babel):", len(rmsd_list_ob))
print("Mean RMSD (RDKit):", np.mean(rmsd_list_rdkit))
print("Mean RMSD (Open Babel):", np.mean(rmsd_list_ob))
print("Mean bond error:", np.mean(bond_list))
print("Mean angle error:", np.mean(angle_list))

plt.title("Distribution of RMSD ({})".format(sys.argv[4]))
plt.xlabel("RMSD")
plt.ylabel("Number of molecules")
plt.hist(rmsd_list_rdkit, bins=20)
plt.savefig(sys.argv[3])

