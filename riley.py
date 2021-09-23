from ast import literal_eval
import farey_words
from numpy.polynomial import Polynomial as P
import numpy as np

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

def riley_slice_from_fast_farey_elliptic(a,b,max_denom):
    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey_words.farey_coefficients_fast_elliptic(p,q,a,b) + 2
          points.extend(poly.roots())
    return points
