import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
from matplotlib import rcParams

rc('text', usetex=True)
rc('text.latex', preamble='\usepackage{sfmath}')
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

rc('text.latex', preamble='\usepackage{sfmath}')
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

data = [[  0.577,  .505,  0.419, .499, .479,  .518, .396, 0.248, .300, .253],
        [  0.497, .335,  0.365,  .417, .464,  .468, .369, 0.320, .244, .192]
        ]

labels = [r'\mbox{ $\nu$}', r'\mbox{ $\tilde \nu$}']
data = np.transpose([[row[i] for row in data] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]])

colors = plt.cm.RdYlGn(np.linspace(0, 1, 4))
x0 = np.arange(10)/5.0

ix = [0, 1]
deltax = [.025, .1]

x_table = np.array([.0] * 10)

cell_text = []
for i,dx,l in zip(ix, deltax, labels):
    plt.bar(x0+dx, data[i], .075, color=colors[i], label = l)
    cell_text.append(['%1.2f' % (x) for x in data[i]])

ysig = [0.5096, 0.5096]
plt.plot([0, 2.05], ysig, linestyle= '--', dashes=[30, 12, 30, 12], lw = 3, color='black')
plt.fill_between([0, 2.05], [0.429, 0.429], [0.584, 0.584], facecolor='black', alpha = 0.17, interpolate=True)
plt.fill_between([0, 2.05], [0.585, 0.585], [0.8, 0.8], facecolor='lime', alpha = 0.32, interpolate=True)
#ysig = [0.5109, 0.5109]
#plt.plot([0, 1.85], ysig, linestyle= '--', dashes=[30, 12, 30, 12], lw = 3, color='black')

yticklabels = getp(gca(), 'yticklabels')
setp(yticklabels, fontsize=31, fontname = 'Arial')

plt.legend(prop={'size':26}, loc = 'upper right', borderaxespad=0.1, bbox_to_anchor = (1, 0.7))

plt.text(0.5, 0.67, 'Significantly even', fontsize = 30)

plt.text(0.8, 0.53, 'Random', fontsize = 30)

ax = pylab.axes() 
ax.xaxis.labelpad = 26

ylabel(r'Evenness Score \mbox{ $\left( E_{1/D}\right)$}', fontsize = 34)
xlabel('Dataset', fontsize = 34)

plt.text(0.103, -0.03, 'LS', ha='center', va='center', fontsize = 29)
plt.text(0.302, -0.03, 'LA', ha='center', va='center', fontsize = 29)
plt.text(0.507, -0.03, 'la', ha='center', va='center', fontsize = 29)
plt.text(0.711, -0.03, 'RC', ha='center', va='center', fontsize = 29)
plt.text(0.912, -0.03, 'rc', ha='center', va='center', fontsize = 29)
plt.text(1.113, -0.03, 'SI', ha='center', va='center', fontsize = 29)
plt.text(1.315, -0.03, 'M', ha='center', va='center', fontsize = 29)
plt.text(1.516, -0.03, 'N1', ha='center', va='center', fontsize = 29)
plt.text(1.716, -0.03, 'N2', ha='center', va='center', fontsize = 29)
plt.text(1.916, -0.03, 'N3', ha='center', va='center', fontsize = 29)

xlim([0, 1.996])

ylim([0, 0.8])

plt.xticks([])
plt.show()
