import numpy as np
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

reps = 900000
depth = 20




mu = riley.cusp_point(p,q,r,s)

alpha = np.csingle(np.exp(2j*np.pi/p))
beta = np.csingle(np.exp(2j*np.pi/q))
X = np.array([[alpha,1],[0,np.conj(alpha)]])
Y = np.array([[beta,0],[mu,np.conj(beta)]])

seeds = farey.fixed_points(0,1,mu,alpha,beta)\
      #+ farey.fixed_points(1,1,mu,alpha,beta)\
      #+ farey.fixed_points(1,2,mu,alpha,beta)

ls = kleinian.limit_set_markov([X,Y],np.array(seeds),depth,True,reps)

ls = list(itertools.chain(*[[(point, points[1]) for point in points[0] if np.abs(point)<3] for points in ls]))

colours = {-2: 'r', -1:'b', 1:'g', 2:'y'}
plt.scatter([np.real(t[0]) for t in ls],[np.imag(t[0]) for t in ls],c=[colours[t[1]] for t in ls],s=.1,alpha=.1,marker='.',linewidths=0)

plt.axis('equal')
plt.axis([-2,2,-1,2])
plt.tight_layout()
#plt.show()
plt.savefig('cusp12_elliptic33_large_fewseeds.png',dpi=2000)