import farey
import numpy as np

# We hit Python's maximum float value at 745.
max_denom = 300

with open('fareys_34.txt','w') as f:
    f.write('{\n')
    for q in range(1,max_denom+1):
        for p in range(1,q+1):
            if np.gcd(p,q) == 1:
                poly = farey.polynomial_coefficients_fast(p, q, np.exp(1j*np.pi/3), np.exp(1j*np.pi/4)) + 2
                f.write(f'{{{{{p},{q}}},')
                poly = f'{poly:ascii}'.replace('**', '^').replace(' ','').replace('*',' ')
                f.write(f'{poly}}}')
                f.write('\n' if q == max_denom and p == q-1 else ',\n')
    f.write('}\n')
