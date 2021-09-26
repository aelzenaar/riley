import numpy as np
import riley
import matplotlib.pyplot as plt

ls = riley.riley_slice_from_fast_farey_accurate_elliptic(50,3,4)
plt.scatter([np.real(t) for t in ls],[np.imag(t) for t in ls],marker=".",s=5,linewidths=0,c='k')
plt.axis('equal')
plt.tight_layout()
plt.show()
