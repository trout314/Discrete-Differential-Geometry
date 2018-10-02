import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from random import choice
from math import sqrt

filenames = ['edge_graph_S3_16k_V.dat',
	'edge_graph_S3_16k_VG_diam6.dat',
	'edge_graph_S3_16k_VGL_L01_diam14.dat']

#filenames = ['edge_graph_S3_64k_VGL_L03_diam16.dat',
#	'edge_graph_S3_128k_VGL_L03_diam17.dat']

lists = []

for filename in filenames:
	g = nx.Graph()
	edges = nx.read_edgelist(filename)
	g.add_edges_from(edges.edges())

	distances = []

	samples = 100000
	for i in range(samples):
		node1 = choice(list(g.nodes))
		node2 = choice(list(g.nodes))
		while (node1 == node2):
			node2 = choice(list(g.nodes))
		distances.append(nx.shortest_path_length(g, node1, node2))	

	nDist = float(len(distances))
	meanDist = sum(distances) / nDist
	meanDistSq = sum([d**2 for d in distances]) / nDist
	stdDev = sqrt(meanDistSq - meanDist**2)

	print(filename, meanDist, stdDev)

	lists.append(distances)


font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

fig, axes = plt.subplots(nrows=3, ncols=1)
ax0, ax1, ax2 = axes.flatten()
plt.figure(num=1, figsize=(8, 6))

maxBins = 12
myBins = range(1, maxBins)

x = ax0.hist(lists[0], bins=myBins, align='right')
ax0.set_title("distance $d(v,w)$ for {} randomly chosen $v$, $w$".format(samples))
ax0.set_xticks(myBins)
ax0.set_ylim([0, max(x[0])*1.1])
ax0.text(0.57, 0.50, 'volume: 64,000 ' + r'($\alpha = 0.01)$', transform=ax0.transAxes, fontsize=11)

x = ax1.hist(lists[1], bins=myBins, align='right')
ax1.set_xticks(myBins)
ax1.set_ylim([0, max(x[0])*1.1])
ax1.text(0.57, 0.50, 'volume: 64,000 ' + r'($\alpha = 0.01$)' + '\nmean edge degree: 5.1 ' + r'($\beta = 0.01$)', transform=ax1.transAxes, fontsize=11)

x = ax2.hist(lists[2], bins=myBins, align='right')
ax2.set_xticks(myBins)
ax2.set_ylim([0, max(x[0])*1.1])
ax2.text(0.57, 0.50, 'volume: 64,000 ' + r'($\alpha = 0.01$)' + '\nmean edge degree: 5.1 ' + r'($\beta = 0.01$)' + '\nminimize edge deg var ' + r'($\gamma = 0.1$)', transform=ax2.transAxes, fontsize=11)

fig.text(0.02, 0.5, 'count', ha='center', va='center', rotation='vertical')


every_nth = 2
for ax in [ax0, ax1, ax2]:
	for n, label in enumerate(ax.yaxis.get_ticklabels()):
    		if n % every_nth != 0:
        		label.set_visible(False)

plt.savefig('testfig.png')
plt.show()

