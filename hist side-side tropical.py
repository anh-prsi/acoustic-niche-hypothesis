import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from matplotlib import rc
from matplotlib import rcParams
rc('text', usetex=True)
rc('text.latex', preamble='\usepackage{sfmath}')
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

pylab.rcParams['xtick.major.pad']='4'
pylab.rcParams['ytick.major.pad']='4'

matplotlib.rcParams['xtick.labelsize'] = 30
matplotlib.rcParams['ytick.labelsize'] = 30
matplotlib.rcParams['axes.labelsize'] = 38

font0 = FontProperties()

number_of_bins = 30

# An example of three data sets to compare
number_of_data_points = 1000
labels = ["[LS]", "[LA]", "[la]", "[RC]", "[rc]", "[SI]", "[M]"]

data_sets = [np.array([4.26000000e+02,9.46000000e+02,1.03800000e+03,1.11600000e+03,1.22000000e+03,1.28300000e+03,1.40300000e+03,1.78200000e+03,2.02200000e+03,2.05800000e+03,2.07500000e+03,2.12000000e+03,2.21000000e+03,2.32400000e+03,2.42900000e+03,2.64200000e+03,3.09500000e+03,3.19000000e+03,3.23900000e+03,3.24200000e+03,3.32200000e+03,3.52900000e+03,3.97000000e+03,3.97900000e+03,4.14600000e+03,4.59200000e+03,6.13000000e+03,6.43000000e+03,6.46000000e+03,7.98200000e+03,8.84800000e+03,5.78500000e+02,NaN,NaN,NaN,1.25950000e+03,1.73850000e+03,1.85250000e+03,1.78800000e+03,2.46250000e+03,2.56600000e+03,2.70850000e+03,2.79600000e+03,3.67400000e+03,3.71700000e+03,3.98266667e+03,4.28500000e+02,4.74300000e+03,5.42500000e+02,5.96500000e+02,6.66000000e+02,1.97100000e+03,1.89800000e+03,5.76000000e+02,3.09000000e+02,1.60700000e+03,2.01600000e+03,3.93000000e+02,5.72000000e+02,6.46000000e+02,2.11800000e+03,3.66300000e+03,3.81500000e+03,4.40100000e+03,5.08500000e+03,1.47125000e+03,1.91450000e+03,2.06400000e+03,2.73800000e+03,3.21900000e+03,2.17100000e+03,1.04300000e+03,3.31150000e+03,2.22700000e+03,1.90000000e+03,3.74550000e+03,2.78900000e+03,1.19100000e+03,2.33800000e+03,3.70000000e+03,8.72900000e+03,NaN,NaN,3.42100000e+03,3.46600000e+03,5.53000000e+02,NaN,1.95800000e+03,2.43500000e+02,2.89300000e+03,3.98000000e+03,3.28850000e+03,1.15900000e+03,1.88100000e+03,3.04600000e+03,7.83000000e+02,4.42200000e+03,NaN,2.57300000e+03,1.51300000e+03,1.52800000e+03,2.27300000e+03,3.03000000e+03,3.08000000e+03,3.39600000e+03,3.42300000e+03,3.46000000e+03,3.66500000e+03,3.99800000e+03,4.37900000e+03,4.77500000e+03,1.71100000e+03,2.34500000e+03,3.09500000e+03,3.86200000e+03,5.60000000e+02,8.66666667e+02,9.94333333e+02,6.23000000e+02,5.30000000e+03,7.01700000e+03,2.93400000e+03,3.61250000e+03,1.30650000e+03,3.28800000e+03,2.29800000e+03,3.60000000e+03,NaN,NaN,3.09600000e+03,7.38600000e+03,NaN,1.72150000e+03,2.18150000e+03,3.51125000e+03,1.41100000e+03,1.48600000e+03,2.14500000e+03,2.26800000e+03,2.38000000e+03,2.52500000e+03,3.37000000e+03,7.57100000e+03,1.27033333e+03,2.27150000e+03]),
			 np.array([1629, 639, 3266, 4993, 3616, 4696, 6274, 2710, 4968, 2210, 1873, 4891, 2649, 5830, 1879, 3524, 2752, 2990, 7564, 2005, 1653, 2263, 6459, 2872, 2303, 2342, 4781, 3204, 2835, 7804, 2208, 2668, 3717, 6234, 5655, 6035, 2979, 2256, 5303, 1575, 3933, 2161, 2629, 1002, 886, 7166, 4524, 2669, 2376, 5583, 1564, 2554, 6288, 1977, 3740, 3021, 2750, 1716, NaN, 7003, 2148, 786, 4404, 2499, 1296, 3418, 1266, 2698, 1898, 6440, 2242]),
             np.array([7.02900000e+03,9.87700000e+03,4.15400000e+03,6.98300000e+03,3.33500000e+03,2.60300000e+03,2.37700000e+03,3.81400000e+03,2.91400000e+03,2.37300000e+03,2.29700000e+03,2.68400000e+03,2.13400000e+03,3.87000000e+03,7.93000000e+03,3.90200000e+03,1.89500000e+03,5.59400000e+03,1.72700000e+03,2.86100000e+03,4.20100000e+03,2.29200000e+03,3.67500000e+03,4.08800000e+03,4.51700000e+03,3.83800000e+03,3.40700000e+03,4.00400000e+03,2.78100000e+03,4.53200000e+03,3.45200000e+03,3.60800000e+03,3.02400000e+03,6.56800000e+03,4.66800000e+03,3.58300000e+03,3.21600000e+03,7.31700000e+03,3.08600000e+03,5.74900000e+03,3.53100000e+03,2.10500000e+03,3.68500000e+03,4.29700000e+03,2.50700000e+03,3.29200000e+03,1.64200000e+03,3.63700000e+03,1.23500000e+03,2.37100000e+03,3.27300000e+03,2.75600000e+03,8.44100000e+03,2.31100000e+03,2.40200000e+03,2.80900000e+03,2.56200000e+03,1.80900000e+03,2.56500000e+03,8.14400000e+03,3.59600000e+03,2.36900000e+03,3.46400000e+03,3.67300000e+03,2.79300000e+03,2.92400000e+03,3.19900000e+03,3.82600000e+03,2.61700000e+03,8.68700000e+03,1.32100000e+03,1.72700000e+03,3.41700000e+03,2.50700000e+03,5.46000000e+02,2.59400000e+03,2.04800000e+03,1.89000000e+03,2.37500000e+03,3.89500000e+03,7.37900000e+03,2.19200000e+03,4.01400000e+03,2.79700000e+03,3.61800000e+03,2.57400000e+03,3.88600000e+03,3.21100000e+03,1.88100000e+03,2.58400000e+03,2.84200000e+03,5.06700000e+03,1.91600000e+03,2.18200000e+03,2.34500000e+03,2.28200000e+03,1.93300000e+03,1.50700000e+03,2.75600000e+03,2.52700000e+03,2.62200000e+03,3.45500000e+03,2.30600000e+03,1.29200000e+03,2.20100000e+03,2.12500000e+03,2.51200000e+03,2.60300000e+03,2.79900000e+03,7.18000000e+02,1.55000000e+03,2.89500000e+03,2.19600000e+03,2.36900000e+03,1.95200000e+03,1.82300000e+03,1.43600000e+03,3.26300000e+03,6.75000000e+02,1.48300000e+03,1.45000000e+03,1.99100000e+03,2.79500000e+03,2.17300000e+03,2.19200000e+03,1.70100000e+03,1.24000000e+03,4.16000000e+02,2.19600000e+03,1.72300000e+03,8.76000000e+02,1.70800000e+03,2.70500000e+03,3.12500000e+03,3.69500000e+03,1.99100000e+03,1.44500000e+03,1.96900000e+03,1.55000000e+03,2.10200000e+03,2.13600000e+03,1.23000000e+03,1.85700000e+03,2.13900000e+03,2.67000000e+03,3.54000000e+02,1.24400000e+03,2.37800000e+03,1.95200000e+03,4.45000000e+02,2.12600000e+03,2.00000000e+03,2.46900000e+03,4.12000000e+02,1.72300000e+03,4.98000000e+02,7.64000000e+02,1.68000000e+03,3.83000000e+02,1.91400000e+03,1.94300000e+03,3.06700000e+03,1.96700000e+03,1.81400000e+03,1.99200000e+03,1.79000000e+03,5.17000000e+02,1.20600000e+03,1.29200000e+03,2.04800000e+03,1.27500000e+03,4.91000000e+03,1.49300000e+03,1.68400000e+03,1.38500000e+03,1.21100000e+03,1.50500000e+03,1.80900000e+03,1.31100000e+03,1.24000000e+03,2.23900000e+03,1.78000000e+03,1.23500000e+03,1.48500000e+03,1.55600000e+03,1.18200000e+03,4.02000000e+02,1.45700000e+03,9.67000000e+02,1.45500000e+03,1.09100000e+03,8.71000000e+02,1.49800000e+03,1.18700000e+03,1.72000000e+02,1.51600000e+03,3.50300000e+03]),
			 np.array([2.2465e+03, 1.3430e+03, 3.0610e+03, 2.5180e+03, 2.1390e+03, 2.8780e+03, 3.8130e+03, 2.5330e+03, 3.4890e+03, 1.4437e+03, 1.8810e+03, 2.9015e+03, 3.7370e+03, 1.7710e+03, 3.0720e+03, 1.0460e+03, 2.7800e+03, 1.9810e+03, 3.4280e+03, 3.6485e+03, 2.4360e+03, 2.2360e+03, 2.1300e+03, 3.1670e+03, 2.6840e+03, 1.5650e+03, 0.3150e+03, 2.1080e+03, 3.7700e+03, 1.1490e+03, 1.6340e+03, 2.3375e+03, 2.2660e+03, 1.6530e+03, 2.5020e+03, 2.3385e+03, 3.9920e+03, 1.9945e+03, 3.2774e+03, 4.5160e+03, 3.3460e+03, 4.0990e+03, 3.2615e+03, 0.5420e+03, 2.6860e+03, 3.5490e+03, 1.1555e+03, 2.0560e+03, 2.1517e+03, 0.4460e+03, 4.5122e+03, 4.2040e+03, 6.7355e+03, 2.6480e+03, 3.3140e+03, 2.2945e+03, 1.7150e+03, 1.7670e+03, 1.6240e+03, 4.6343e+03, 2.5270e+03, 2.2980e+03, 0.6094e+03, 0.8000e+03, 3.0830e+03, 2.8700e+03, 1.6460e+03, 2.5837e+03, 2.8805e+03, 1.3427e+03, 1.8710e+03, 1.7815e+03, 1.6435e+03, 2.9830e+03, 1.3490e+03, 1.7043e+03, 1.4100e+03, 1.4840e+03, 1.4490e+03, 3.3250e+03, 4.6860e+03, 1.2920e+03, 3.4500e+03, 1.7070e+03, 2.6710e+03, 1.1420e+03, 0.1350e+03, 2.0475e+03, 2.0780e+03, 2.0285e+03, 3.7968e+03, 1.3920e+03, 3.2137e+03, 4.0735e+03, 2.8785e+03, 5.5710e+03, 6.3325e+03, 2.3912e+03, 7.1730e+03, 3.1185e+03, 3.1910e+03, 2.6310e+03, 2.5295e+03, 3.0750e+03]),
             np.array([1750.67000000000, 1750.67000000000, 1406.67000000000, 1779.33000000000, 1148.33000000000, 803.670000000000, 631, 711, 516, 660, 1349, 1550, 889.670000000000, 1378, 1233.67000000000, 401.330000000000, 2612, 2382.33000000000, 1980.67000000000, 1205, 1406.33000000000, 1894.33000000000, 2210.33000000000, 3014, 4909, 2296.33000000000, 2468.67000000000, 2124.33000000000, 1808, 2038, 3129, 2210.33000000000, 2956.67000000000, 2526, 2066.67000000000, 2784.33000000000, 5799.33000000000, 4449.67000000000, 1492.33000000000, 1836.67000000000, 1722, 2325, 2468.33000000000, 4248.67000000000, 6028.67000000000, 3330, 4076.33000000000, 3014, 3530.67000000000, 1779.67000000000, 1176.33000000000, 2439.67000000000, 3703, 2669.67000000000, 3416, 4047.67000000000, 2956.67000000000, 2382.67000000000, 2038.33000000000, 487.330000000000, 861, 3387.67000000000, 2181.67000000000, 1779.33000000000, 2755.67000000000, 3846.67000000000, 1550, 3129, 4995.33000000000, 3875.33000000000, 2296.33000000000, 1406.67000000000, 2669.33000000000, 2928, 3990.33000000000, 2468.33000000000, 1980.67000000000, 2038.33000000000, 2956.67000000000, 2353.67000000000, 3157.67000000000, 3559.67000000000, 2152.67000000000]),
             np.array([380,440,590,1100,1120,1290,1320,1357,1400,1540,1550,1550,1600,1600,1645,1709,1791,2000,2150,2400,2430,2500,2550,2600,2700,2700,2720,2792,2900,2950,3000,3100,3130,3300,3300,3700,3741,3800,3800,3800,3963.5,4000,4060,4074,4200,4200,4300,4350,4400,4500,4960,5020,5260,5300,5500,5700,5900,6000,6020,6250,6300,6300,6300,6420,6425,6650,6700,6910,6950,6990,7080,7200,7400,7400,7440,7543,7640,8080,8350,9370,9400,10200]),
             np.array([6.29999999e+03,5.10000001e+03,2.69999999e+03,3.80000002e+03,2.40000001e+03,4.59999999e+03,5.30000000e+03,3.70000001e+03,6.29999999e+03,7.00000001e+03,7.60000003e+03,2.00000000e+03,6.40000000e+03,4.99999999e+03,3.00000001e+03,7.09999998e+03,3.10000000e+03,6.40000000e+03,5.39999999e+03,7.30000002e+03,9.40000002e+03,7.70000002e+03,3.20000000e+03,1.10000000e+03,4.40000000e+03,7.40000001e+03,4.70000001e+03,2.90000001e+03,2.50000000e+03,4.09999999e+03,5.49999999e+03,4.00000000e+03,1.70000000e+03,2.20000000e+03,1.70000000e+03,7.99999999e+03,4.09999999e+03,4.89999998e+03,7.00000001e+03,6.60000001e+03,7.99999999e+03,4.99999999e+03,2.69999999e+03,7.50000000e+03,8.69999996e+03,5.30000000e+03,3.00000001e+03,6.40000000e+03,3.70000001e+03,3.80000002e+03,1.50000000e+03,2.40000001e+03,2.69999999e+03,3.89999999e+03,9.30000001e+03,4.00000000e+03,7.99999999e+03,6.29999999e+03,6.90000000e+03,2.59999999e+03,1.10000000e+03,1.70000000e+03,6.90000000e+03,1.02000000e+04,7.50000000e+03,4.30000001e+02,1.30000000e+03,1.35000000e+03,6.20000001e+02,3.79999999e+02,1.60000000e+03,1.45000001e+03,1.55000000e+03,1.55000000e+03,1.30000000e+03])
             ]

# Computed quantities to aid plotting
hist_range = (0.1, 10)
binned_data_sets = [np.histogram(d/1000, range=hist_range, bins=number_of_bins)[0]
                    for d in data_sets]
binned_data_sets[0]=binned_data_sets[0]/float(sum(binned_data_sets[0]))
binned_data_sets[1]=binned_data_sets[1]/float(sum(binned_data_sets[1]))
binned_data_sets[2]=binned_data_sets[2]/float(sum(binned_data_sets[2]))
binned_data_sets[3]=binned_data_sets[3]/float(sum(binned_data_sets[3]))
binned_data_sets[4]=binned_data_sets[4]/float(sum(binned_data_sets[4]))
binned_data_sets[5]=binned_data_sets[5]/float(sum(binned_data_sets[5]))
binned_data_sets[6]=binned_data_sets[6]/float(sum(binned_data_sets[6]))
binned_maximums = np.max(binned_data_sets, axis=1)
x_locations = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5])

# The bin_edges are the same for all of the histograms
bin_edges = np.linspace(hist_range[0], hist_range[1], number_of_bins + 1)
centers = bin_edges[:len(bin_edges)-1] + 0.5*diff(bin_edges)[0]
heights = np.diff(bin_edges)

# Cycle through and plot each histogram
ax = plt.subplot(111)
counter = 1;
for x_loc, binned_data in zip(x_locations, binned_data_sets):
    lefts = x_loc - .5 * binned_data
    if counter <6:
		ax.barh(centers, binned_data, height=heights, left=lefts, color = 'navy')
    else:
		ax.barh(centers, binned_data, height=heights, left=lefts, color = 'limegreen')	
    counter = counter+1

plot(np.array([-0.1, 1.13]), np.array([1, 1]), color='darkred', linestyle= '--', dashes=[20, 8, 20, 8], lw = 4)
plot(np.array([-0.1, 1.13]), np.array([4, 4]), color='darkred', linestyle= '--', dashes=[20, 8, 20, 8], lw = 4)
plot(np.array([1.13, 1.6]), np.array([1, 1]), color='darkred', linestyle= '--', dashes=[20, 8, 20, 8], lw = 4)
plot(np.array([1.13, 1.6]), np.array([8, 8]), color='darkred', linestyle= '--', dashes=[20, 8, 20, 8], lw = 4)

ax.set_xticks(x_locations)
ax.set_xticklabels(labels)

ax.set_ylabel("Frequency(KHz)")
ax.set_xlabel("Data sets")

ax.set_ylim((0, 10))
ax.set_xlim((-0.1, 1.6))

plt.show()
