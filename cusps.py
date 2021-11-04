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
plt.savefig('cusp12_elliptic35_large.png',dpi=2000)
