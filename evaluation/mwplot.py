import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")
sns.set_context("paper", 1.5)

df = pd.read_csv(sys.argv[1])

h = sns.jointplot(df['MW'], df['RMSD'], kind="scatter", marker='.')
h.ax_joint.set_xlabel("Molecular weight (Da)")
h.ax_joint.set_ylabel("RMSD ($\AA$)")
plt.xlim([100, 650])
plt.ylim([0, 7])
plt.tight_layout()
plt.savefig("mw-rdkit.pdf")
