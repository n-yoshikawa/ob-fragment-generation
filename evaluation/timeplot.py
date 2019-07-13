import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

sns.set(style="white")
sns.set_context("paper", 1.8)

df = pd.read_csv(sys.argv[1])

slope, intercept, r_value, p_value, std_err = stats.linregress(df['MW'],df['time'])

g = sns.jointplot(df['MW'], df['time'], kind="reg", joint_kws={'color':'k'}, 
                  scatter=False, xlim=(100, 650), ylim=(0, 50000))
sns.kdeplot(df['MW'], df['time'], ax=g.ax_joint, cmap='Blues')
g.ax_joint.set_xlabel("Molecular weight (Da)")
g.ax_joint.set_ylabel("Execution time (s)")
#g.ax_joint.legend_.remove(
g.ax_joint.legend_.texts[0].set_text("$R^2$: {:.2f}, y={:.3f}x{:.3f}".format(r_value**2, slope, intercept))
plt.tight_layout()
plt.show()
plt.savefig("time-ob-new.pdf")
