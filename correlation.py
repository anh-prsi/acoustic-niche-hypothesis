from string import letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="white")

rc('text.latex', preamble='\usepackage{sfmath}')
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

pylab.rcParams['xtick.major.pad']='4'
pylab.rcParams['ytick.major.pad']='4'

matplotlib.rcParams['xtick.labelsize'] = 27
matplotlib.rcParams['ytick.labelsize'] = 27


corr =    [[   1,  0.86,   0.87,  0.74],
		   [0.86,   1,   0.18,  0.65],
		   [0.87,  0.18,   1,  0.95],
		   [0.74,  0.65,  0.95,   1]
		   ]

f, ax = plt.subplots(figsize=(11, 9))

cmap = sns.diverging_palette(220, 10, as_cmap=True)

sns.heatmap(corr, cmap = cmap,
            square=True, annot=True, annot_kws={"size": 22},
            ax=ax, xticklabels=(['G', 'M', 'T', 'Y']), yticklabels=(['G', 'M', 'T', 'Y']), cbar_kws={"shrink": .875})

ax.set_xlabel('Site', fontsize = 33)
ax.set_ylabel('Site', fontsize = 33)

plt.show()

#f.suptitle(r' \mbox { $ \bar S_4 $} ', fontsize = 36)

f.suptitle('[LS]', fontsize = 46)
