""" Script to plot limit sets of groups following Grandma's recipe from Chapter 8 of Indra's Pearls.
    Modify ta and tb as in Box 21 of Indra's Pearls (p.229).
"""

import mpmath as mp
mp.dps = 50
import kleinian
import matplotlib.pyplot as plt

# Filename to save the final picture to
filename = 'indra.png'

ta = 2+0.2j
tb = 2-0.2j
tab = ((ta*tb)+mp.sqrt(ta**2*tb**2-4*(ta**2+tb**2)))/2
z0 = ((tab-2)*tb)/(tb*tab-2*ta+2j*tab)

# A list of generators for the Kleinian group
generators = [mp.matrix([[ta/2,(ta*tab-2*tb+4j)/((2*tab+4)*z0)],[(ta*tab-2*tb-4j)*z0/(2*tab-4),ta/2]]),mp.matrix([[(tb-2j)/2,tb/2],[tb/2,(tb+2j)/2]])]

# Compute the limit set, with the Markov algorithm, using words
# of maximum length 8; compute 100000 words altogether.
ls = kleinian.limit_set_markov(generators,mp.matrix([0]),9,100000)

# Crop to just -2 < x < 2
ls = [t for t in ls if t[0].real > -2 and t[0].real < 2]

# Plot the limit points, with colours according to the generators
# starting each word. Each point returned by limit_set_markov is
# actually a pair (A,B) where A is the actual point in $\C$ and
# B is a small integer giving the index of the first letter of the
# word giving that point in the generator matrix (or the negative
# of the index, if the word starts with an inverse matrix),
colours = {-2: 'r', -1:'b', 1:'g', 2:'y'}
plt.scatter([t[0].real for t in ls],
            [t[0].imag for t in ls],
            c=[colours[t[1]] for t in ls],
            marker=".", s=0.1, linewidths=0)

plt.axis('equal')
plt.axis([-2,2,-3,3])
plt.tight_layout()

plt.savefig(filename,dpi=500)
plt.show()
