""" A script to plot the cusp points of a Riley slice.

    Options to change:
      ordera, orderb -- the orders of the cone points on the surface (set both to mp.inf for the parabolic slice).
      max_denom -- we plot all the cusps with slope p/q with q <= max_denom.
"""

import mpmath as mp
import riley
import matplotlib.pyplot as plt


# OPTIONS
ordera = mp.inf
orderb = mp.inf
max_denom = 20
#


cusps = []
for q in range(1,max_denom+1):
    for p in range(1,2*q+1):
        if mp.gcd(p,q) == 1:
          try:
              cusps.append(riley.cusp_point(ordera,orderb,p,q))
          except RuntimeError as e:
              pass
              #raise RuntimeError(f'cusp computation failed at (p,q) = {(p,q)}') from e

plt.scatter([mp.re(t) for t in cusps],[mp.im(t) for t in cusps],marker=".",s=5,linewidths=0,c='k')
plt.axis('equal')
plt.axis([-4,4,-3,3])
plt.tight_layout()
plt.show()
