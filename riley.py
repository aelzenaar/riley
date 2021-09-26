from ast import literal_eval
import farey_words
from numpy.polynomial import Polynomial as P
import numpy as np
import scipy.optimize

def riley_slice():
    with open('riley_slice.csv') as f:
      return [literal_eval(line) for line in f]

def riley_slice_from_fast_farey(max_denom):
    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey_words.farey_coefficients_fast_parabolic(p,q) + 2
          points.extend(poly.roots())
    return points

def riley_slice_from_fast_farey_accurate(max_denom):
    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey_words.farey_coefficients_fast_parabolic(p,q) + 2
          roots_bad = poly.roots()
          try:
            roots = [scipy.optimize.newton(poly, root, poly.deriv(),fprime2=poly.deriv(2),maxiter=10000,tol=1e-4) for root in roots_bad]
          except RuntimeError as e:
            print(str(e))
            print(f'Newton failed to converge at {p}/{q}')
            return points
          points.extend(roots)
    return points

def riley_slice_from_fast_farey_accurate_elliptic(max_denom,a,b):
    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey_words.farey_coefficients_fast_elliptic(p,q,a,b) + 2
          roots_bad = poly.roots()
          try:
            roots = [scipy.optimize.newton(poly, root, poly.deriv(),fprime2=poly.deriv(2),maxiter=100000,tol=1e-2) for root in roots_bad]
          except RuntimeError as e:
            print(str(e))
            print(f'Newton failed to converge at {p}/{q}')
            return points
          points.extend(roots)
    return points

def riley_slice_from_fast_farey_elliptic(a,b,max_denom):
    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey_words.farey_coefficients_fast_elliptic(p,q,a,b) + 2
          points.extend(poly.roots())
    return points
