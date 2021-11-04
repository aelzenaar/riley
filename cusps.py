""" A script to plot the limit set of a particular cusp group.

    This is only as accurate as the riley.cusp_point() function allows (i.e. make sure you have
    sacrificed the correct number of goats before running the script...).

    The limit point seeds are the fixed points of the 0/1, 1/1, and 1/2 Farey words at the given cusp point.

    Limit set is computed all at once and ploted with pyplot so if the numbers you try to run are too big
    then you will run out of memory very quickly. Try cusps_shaded.py if this is bad.

    Options to change:
        p, q -- orders of the elliptic elements (set to mp.inf for the parabolic case)
        r, s -- slope of the desired cusp
        reps -- you should be able to run kleinian.limit_set_markov with reps set to this value.
        depth -- maximum word length to compute orbits with

    Output image filename is cusp_{r}_{s}_elliptic_{p}_{q}.png (in the current directory).
"""
from mpmath import mp
import kleinian
import matplotlib.pyplot as plt
import riley
import farey
import itertools

# Orders of elliptics
p = 3
q = 5

# Cusp slope
r = 1
s = 2

reps = 100000
depth = 20




filename = f'cusp_{r}_{s}_elliptic_{p}_{q}.png'
mu = riley.cusp_point(p,q,r,s)

alpha = mp.exp(2j*mp.pi/p)
beta = mp.exp(2j*mp.pi/q)
X = mp.matrix([[alpha,1],[0,mp.conj(alpha)]])
Y = mp.matrix([[beta,0],[mu,mp.conj(beta)]])

seeds = farey.fixed_points(0,1,mu,alpha,beta)\
      + farey.fixed_points(1,1,mu,alpha,beta)\
      + farey.fixed_points(1,2,mu,alpha,beta)

ls = kleinian.limit_set_markov([X,Y],mp.matrix(seeds),depth,reps)

colours = {-2: 'r', -1:'b', 1:'g', 2:'y'}
print((ls.rows,ls.cols))
plt.scatter([mp.re(t[0]) for t in ls],[mp.im(t[0]) for t in ls],c=[colours[t[1]] for t in ls],s=.1,alpha=.1,marker='.',linewidths=0)

plt.axis('equal')
plt.axis([-2,2,-1,2])
plt.tight_layout()
plt.show()
plt.savefig(filename,dpi=2000)
