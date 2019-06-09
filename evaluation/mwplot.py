import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white")

df = pd.read_csv(sys.argv[1])

plt.plot(df['MW'], df['RMSD'], '.')
plt.xlabel("Molecular weight (Da)")
plt.xlim([100, 650])
plt.ylabel("RMSD ($\AA$)")
plt.ylim([0, 7])
plt.show()
