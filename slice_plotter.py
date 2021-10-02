import numpy as np
import riley
import matplotlib.pyplot as plt

ls = riley.riley_slice(np.inf,np.inf,100) # See documentation for this function in riley.py.

plt.scatter([np.real(t) for t in ls],[np.imag(t) for t in ls],marker=".",s=5,linewidths=0,c='k')
plt.axis('equal')
plt.axis([-4,4,-3,3])
plt.tight_layout()
plt.show()
