""" Methods to compute globally and locally with the Riley slice.

    In this file are found methods to compute with the moduli space known as the `Riley slice'.
"""

from ast import literal_eval
import farey
from numpy.polynomial import Polynomial as P
import numpy as np
import scipy.optimize

def riley_slice(a, b, max_denom, max_iter=100, tol=1e-2):
    """ Return an accurate approximation to the Riley slice.

        More precisely, use Farey polynomials to determine an approximation and then attempt
        to further improve the results using Newton's algorithm.

        Arguments:
          a, b -- the order of the cone points represented by X and Y respectively. Use np.inf for the parabolic case (or 1, since exp(2*pi*i/1) = exp(0) = 1).
          max_denom -- the maximum denominator Farey polynomials to compute
          max_iter -- maximum number of iterations for Newton's algorithm (default 100)
          tol -- allowable error of each point (default 1e-2)
    """

    alpha = 1 if a == np.inf else np.exp(2j*np.pi/a)
    beta = 1 if b == np.inf else np.exp(2j*np.pi/b)

    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey.polynomial_coefficients_fast(p,q,alpha,beta) + 2
          roots_bad = poly.roots()
          try:
            roots = [scipy.optimize.newton(poly, root, poly.deriv(),fprime2=poly.deriv(2),maxiter=max_iter,tol=tol) for root in roots_bad]
          except RuntimeError as e:
            print(str(e))
            print(f'Newton failed to converge at {p}/{q}')
            return points
          points.extend(roots)
    return points
