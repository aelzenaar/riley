import numpy as np
import riley
import matplotlib.pyplot as plt

cusps = []
max_denom = 20
for q in range(1,max_denom+1):
    for p in range(1,2*q+1):
        if np.gcd(p,q) == 1:
          try:
              cusps.append(riley.cusp_point(np.inf,np.inf,p,q))
          except RuntimeError as e:
              pass
              #raise RuntimeError(f'cusp computation failed at (p,q) = {(p,q)}') from e

plt.scatter([np.real(t) for t in cusps],[np.imag(t) for t in cusps],marker=".",s=5,linewidths=0,c='k')
plt.axis('equal')
plt.axis([-4,4,-3,3])
plt.tight_layout()
plt.show()
