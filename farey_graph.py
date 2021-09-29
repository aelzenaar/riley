import networkx as nx
import pydot
import farey
import numpy as np
import random



# Configuration options:
max_denom = 15 # Denominator to plot up to.
plot_uncoloured = False # Plot the uncoloured edges as well as the coloured ones.
only_tree = False # Colour only the edges of the Stern-Brocot tree.
file_name = 'farey_graph_coloured.png'


colour_list = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000']

G = nx.DiGraph()

for p in range(1,max_denom):
    for q in range(2,max_denom):
        if p<q and np.gcd(p,q) == 1:
            for nbr in farey.neighbours(p,q):
                G.add_edge(str(nbr),str((p,q)))

g = nx.drawing.nx_pydot.to_pydot(G)
g.set_dpi(300)

for edge in g.get_edge_list():
    edge.set('pendwidth',20)

def nodename(pair):
    return '"' + str(pair) + '"'

randcol=0
for p in range(1,2*max_denom):
    for q in range(2,2*max_denom):
        if p<q and np.gcd(p,q) == 1:
            (r1,s1),(r2,s2) = farey.neighbours(p,q)
            if only_tree:
                max_edges = 2
                edges = g.get_edge(nodename((p,q)),nodename((r1+p,s1+q)))\
                      + g.get_edge(nodename((p,q)),nodename((r2+p,s2+q)))
            else:
                max_edges = 10
                edges = g.get_edge(nodename((r1,s1)),nodename((p,q)))                     \
                      + g.get_edge(nodename((r1,s1)),nodename((r1+p,s1+q)))               \
                      + g.get_edge(nodename((p,q)),nodename((r1+p,s1+q)))                 \
                      + g.get_edge(nodename((r1+p,s1+q)), nodename((r1+2*p, s1+2*q)))     \
                      + g.get_edge(nodename((p,q)), nodename((r1+2*p, s1+2*q)))           \
                      + g.get_edge(nodename((r2,s2)),nodename((p,q)))                     \
                      + g.get_edge(nodename((r2,s2)),nodename((r2+p,s2+q)))               \
                      + g.get_edge(nodename((p,q)),nodename((r2+p,s2+q)))                 \
                      + g.get_edge(nodename((r2+p,s2+q)),nodename((r2+2*p, s2+2*q)))      \
                      + g.get_edge(nodename((p,q)),nodename((r2+2*p, s2+2*q)))

            if len(edges) < max_edges:
              continue

            for edge in edges:
                col = edge.get('color')
                if col == None:
                    col = colour_list[randcol]
                else:
                    col = col + ':' + colour_list[randcol]
                edge.set('color',col)
            randcol = (0 if randcol == len(colour_list)-1 else randcol+1)

if not plot_uncoloured:
    edges_to_delete = []
    for edge in g.get_edges():
        if edge.get('color') == None:
            edges_to_delete.append((edge.get_source(),edge.get_destination()))
    for (src,dest) in edges_to_delete:
        g.del_edge(src,dest)

g.set_suppress_disconnected(True)

g.write_png(file_name)
