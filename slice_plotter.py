""" Plot a Riley slice with the given parameters

    Options:
      ordera, orderb -- cone angles (set to mp.inf for the parabolic slice)
      maxdenom -- maximum denominator of Farey polynomial to use in the approximation
"""

from mpmath import mp
import riley
import matplotlib.pyplot as plt
from datetime import datetime

ordera = 6
orderb = 8
maxdenom = 100

mp.dps = 100

print(f'Began run at {datetime.now()}')
ls = riley.riley_slice(ordera,orderb,maxdenom) # See documentation for this function in riley.py.

count = len(ls)
print(f'Finished computing points at {datetime.now()}')
print(f'Plotting {count} points.')

plt.scatter([mp.re(t) for t in ls],[mp.im(t) for t in ls],marker=".",s=1, alpha=0.1,linewidths=0,c='k')
plt.axis('equal')
plt.axis([-4,4,-3,3])
plt.tight_layout()
print(f'Finished plotting points at {datetime.now()}')
plt.savefig(f'rileyslice{maxdenom}.png',dpi=2000)
plt.show()
