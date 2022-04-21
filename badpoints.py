""" A script to draw the Farey graph and the Stern-Brocot tree.

    We colour the edges according to the `butterfly' scheme or just the edges of the Stern-Brocot tree
    according to the options set below (see [EMS21b] cited in README.md for a better explanation).

    Options:
      max_denom -- denominator to plot up to.
      plot_uncoloured -- True/False: plot the uncoloured edges as well as the coloured ones.
      only_tree -- True/False: colour only the edges of the Stern-Brocot tree.
      file_name -- filename to output image to.

"""
import farey
import numpy as np
import random
import math


# Configuration options:
max_denom = 15 # Denominator to plot up to.

nodes = {(0,1):1,(1,1):-1,(1,2):1}
print('digraph G {\nrankdir=LR;')
print('"(0, 1)" -> "(1, 1)"\n"(1, 1)" -> "(1, 2)"\n"(0, 1)" -> "(1, 2)"')

for p in range(1,max_denom):
    for q in range(3,max_denom):
        if p<q and np.gcd(p,q) == 1:
            (nbr1,nbr2) = farey.neighbours(p,q)
            print(f'"{str(nbr1)}" -> "{str((p,q))}"')
            print(f'"{str(nbr2)}" -> "{str((p,q))}"')
            nbr3 = (abs(nbr1[0]-nbr2[0]),abs(nbr1[1]-nbr2[1]))
            nodes[(p,q)] = -nodes[nbr3] - nodes[nbr2]*nodes[nbr1] - 2*(nodes[nbr2] + nodes[nbr1])
            print(f'"{str((p,q))}" [label="{str((p,q))}\n{nodes[(p,q)]}" {", color=red, bgcolor=red, style=filled, fontcolor=white" if math.sqrt(abs(nodes[(p,q)])) != math.sqrt(abs(nodes[nbr2]))*math.sqrt(abs(nodes[nbr1])) + math.sqrt(abs(nodes[nbr3])) else ""}]')

print('}')
