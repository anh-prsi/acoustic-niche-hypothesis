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


corr =    [[   1,  0.4780,   0.2240,  0.2050, 0.0660],
		   [0.4780,   1,   0.9210,  0.6880, 0.9890],
		   [0.2240,  0.9210,   1,  0.9230, 0.9110],
		   [0.2050,  0.6880,  0.9230,   1, 0.1310],
		   [0.0660,  0.9890,  0.9110, 0.1310,   1],
		   ]

f, ax = plt.subplots(figsize=(11, 9))

cmap = sns.diverging_palette(220, 10, as_cmap=True)

sns.heatmap(corr, cmap = cmap,
            square=True, annot=True, annot_kws={"size": 18},
            ax=ax, xticklabels=(['B', 'F', 'R', 'S', 'T']), yticklabels=(['B', 'F', 'R', 'S', 'T']), cbar_kws={"shrink": .875})

ax.set_xlabel('Habitat', fontsize = 33)
ax.set_ylabel('Habitat', fontsize = 33)

plt.show()

#f.suptitle(r' \mbox { $ \bar S_4 $} ', fontsize = 36)

f.suptitle('[LA]', fontsize = 46)
