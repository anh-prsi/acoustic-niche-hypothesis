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

data = [[  0.970,  0.161  , 0.305,  .713, 0.216,  1,    0.999,   0.23],
        [  0.970,  0.459  , 0.109,  .472, 0.227, .998,  .713, 0.5]
        ]

labels = [r'\mbox{ $\bar S_3$}', r'\mbox{ $E_{\mathrm{1/D}}$}']
data = np.transpose([[row[i] for row in data] for i in [0, 1, 2, 3, 4, 5, 6, 7]])

colors = plt.cm.BuPu(np.linspace(0.2, 1, 3))
x0 = np.arange(8)/5.0

ix = [0, 1]
deltax = [.025, .1]

cell_text = []
for i,dx,l in zip(ix, deltax, labels):
    plt.bar(x0+dx, data[i], .075, color=colors[i], label = l)
    cell_text.append(['%1.2f' % (x) for x in data[i]])

ysig = [0.95, 0.95]
plt.plot([0, 1.59], ysig, linestyle= '--', dashes=[30, 12, 30, 12], lw = 3, color='black')

yticklabels = getp(gca(), 'yticklabels')
setp(yticklabels, fontsize=31, fontname = 'Arial')

plt.legend(prop={'size':22}, bbox_to_anchor = (0.414, 0.966))

ylabel('Mean percentile', fontsize = 34)

ax = pylab.axes() 
ax.xaxis.labelpad =31

xlabel('Dataset', fontsize = 34)

plt.text(0.103, -0.05, 'LS', ha='center', va='center', fontsize = 33)
plt.text(0.302, -0.05, 'LA', ha='center', va='center', fontsize = 33)
plt.text(0.507, -0.05, 'la', ha='center', va='center', fontsize = 33)
plt.text(0.711, -0.05, 'RC', ha='center', va='center', fontsize = 33)
plt.text(0.912, -0.05, 'rc', ha='center', va='center', fontsize = 33)
plt.text(1.113, -0.05, 'SI', ha='center', va='center', fontsize = 33)
plt.text(1.321, -0.05, 'M', ha='center', va='center', fontsize = 33)
plt.text(1.539, -0.05, 'N2', ha='center', va='center', fontsize = 33)

plt.xticks([])
plt.show()




null_new(f_LS, f_all, [1000, 4000], 'number_av');
null_new(f_SI, f_all, [1000, 8000], 'number_av');
null_new(f_M, f_all, [1000, 8000], 'number_av');
null_new(f_RC, f_all, [1000, 4000], 'number_av');
null_new(f_rc, f_all, [1000, 4000], 'number_av');
null_new(f_la, f_all, [1000, 8000], 'number_av');
null_new([f_LA1, f_LA2], f_all, [1000, 8000], 'number_av');
null_new(f_B, f_all, [1000, 8000], 'number_av');

    0.9500    0.9800
     1     1
    1.0000    0.8790
    0.6830    0.5260
    0.2860    0.2940
         0    0.0050
0.9960    0.9910
    0.8660    0.8580

    
null_new(f_LS, f_all, [1000, 4000], 'spatial_av');
null_new(f_SI, f_all, [1000, 8000], 'spatial_av');
null_new(f_M, f_all, [1000, 8000], 'spatial_av');
null_new(f_RC, f_all, [1000, 4000], 'spatial_av');
null_new(f_rc, f_all, [1000, 4000], 'spatial_av');
null_new(f_la, f_all, [1000, 8000], 'spatial_av');
null_new([f_LA1, f_LA2], f_all, [1000, 8000], 'spatial_av');
null_new(f_B, f_all, [1000, 8000], 'spatial_av');

    0.8350    0.9740
    1     1
    1.0000    0.8590
    0.8860    0.5080
    0.5140    0.3090
         0    0.0040
    0.9940    0.9890
    0.9750    0.8520
    
null0_new(f_LS, f_all, m_LS, m_all, [1000, 4000], 'number_av');
null0_new(f_SI, f_all, m_SI, m_all, [1000, 8000], 'number_av');
null0_new(f_M, f_all, m_M, m_all, [1000, 8000], 'number_av');
null0_new(f_RC, f_all, m_RC, m_all, [1000, 4000], 'number_av');
null0_new(f_rc, f_all, m_rc, m_all, [1000, 4000], 'number_av');
null0_new(f_la, f_all, m_la, m_all, [1000, 8000], 'number_av');

0.9470    0.9730
1.0000    0.9980
0.9980    0.7260
0.6620    0.4890
0.1540    0.1920
     0    0.0060
     
null0_new(f_LS, f_all, m_LS, m_all, [1000, 4000], 'spatial_av');
null0_new(f_SI, f_all, m_SI, m_all, [1000, 8000], 'spatial_av');
null0_new(f_M, f_all, m_M, m_all, [1000, 8000], 'spatial_av');
null0_new(f_RC, f_all, m_RC, m_all, [1000, 4000], 'spatial_av');
null0_new(f_rc, f_all, m_rc, m_all, [1000, 4000], 'spatial_av');
null0_new(f_la, f_all, m_la, m_all, [1000, 8000], 'spatial_av');

0.7940    0.9620
1.0000    0.9980
1.0000    0.7030
0.7700    0.4740
0.2800    0.2120
     0    0.0100
     
     including NT in the data pool...
0.8840    0.9800
0.9990    0.9900
0.9910    0.5630
0.9750    0.7980
0.6660    0.5200
     0    0.0190
