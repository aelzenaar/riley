import numpy as np
import riley
import matplotlib.pyplot as plt
from datetime import datetime

maxdenom = 600

print(f'Began run at {datetime.now()}')
ls = riley.riley_slice(np.inf,np.inf,maxdenom) # See documentation for this function in riley.py.

count = len(ls)
print(f'Finished computing points at {datetime.now()}')
print(f'Plotting {count} points.')

plt.scatter([np.real(t) for t in ls],[np.imag(t) for t in ls],marker=".",s=5, alpha=0.1,linewidths=0,c='k')
plt.axis('equal')
plt.axis([-4,4,-3,3])
plt.tight_layout()
print(f'Finished plotting points at {datetime.now()}')
plt.savefig(f'rileyslice{maxdenom}.png',dpi=2000)
plt.show()
