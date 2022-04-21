import farey
from math import gcd

# LaTeX preamble
print(r'''\begin{table}
  \centering
  \caption{Farey words $ \Word(p/q) $ for small $ q $. See \cref{ex:py_words}.\label{tab:words}}
  \begin{tabular}{r|l}
    $p/q$ &$\Word(p/q)$\\\hline''')

# Do the 0/1 and 1/1 cases separately.
print(f"    $0/1$ & ${''.join(farey.word(0,1))}$\\\\")
print(f"    $1/1$ & ${''.join(farey.word(1,1))}$\\\\")

# Now do the rest.
highest_q = 12
for q in range(2,highest_q+1):
  for p in range(1,q):
    if gcd(p,q) == 1:
      print(f"    ${p}/{q}$ & ${str().join(farey.word(p,q))}$\\\\")

# End the LaTeX code.
print(r'''  \end{tabular}
\end{table}''')
