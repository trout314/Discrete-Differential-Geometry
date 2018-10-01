import networkx as nx

G = nx.Graph()

edges = nx.read_edgelist('edge_graph_S3_64k_VGL_L01.dat')
G.add_edges_from(edges.edges())

print(nx.diameter(G))