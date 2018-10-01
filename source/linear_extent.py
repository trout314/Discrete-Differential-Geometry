import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from random import choice
from math import sqrt

filenames = ['edge_graph_S3_16k_VGL_L01_diam14.dat',
	'edge_graph_S3_23k_VGL_L01_diam12.dat',
	'edge_graph_S3_32k_VGL_L01_diam17.dat',
	'edge_graph_S3_64k_VGL_L01_diam16.dat',
	'edge_graph_S3_91k_VGL_L01.dat']

lists = []

for filename in filenames:
	g = nx.Graph()
	edges = nx.read_edgelist(filename)
	g.add_edges_from(edges.edges())

	distances = []

	for i in range(100000):
		node1 = choice(list(g.nodes))
		node2 = choice(list(g.nodes))
		distances.append(nx.shortest_path_length(g, node1, node2))	

	nDist = float(len(distances))
	meanDist = sum(distances) / nDist
	meanDistSq = sum([d**2 for d in distances]) / nDist
	stdDev = sqrt(meanDistSq - meanDist**2)

	print(filename, meanDist, stdDev)

	lists.append(distances)


fig, axes = plt.subplots(nrows=5, ncols=1)
ax0, ax1, ax2, ax3, ax4 = axes.flatten()

ax0.hist(lists[0], bins=range(12), align='right')
ax0.set_title("d(v, w) for random vertices v, w")

ax1.hist(lists[1], bins=range(14), align='right')
ax2.hist(lists[2], bins=range(14), align='right')
ax3.hist(lists[3], bins=range(14), align='right')
ax4.hist(lists[4], bins=range(14), align='right')

plt.savefig('testfig.png')

