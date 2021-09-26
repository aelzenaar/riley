#ClearAll["Global`*"];
#lookupTable = {{X, x}, {y, Y}};
#word[px_, qx_] := (
   #p = px/GCD[px, qx]; q = qx/GCD[px, qx];
   #length = (2 q/GCD[p, 2 q])*(3 -
       #GCD[2, 2 q/GCD[p, 2 q]]); (*be overly clever*)
   #wd = {};
   #Do[(
     #height = i*p/q;
     #height = If[Ceiling[height] == height, height + 1/2, height];
     #AppendTo[wd,
      #lookupTable[[If[Mod[i, 2] == 1, 1, 2],
       #If[Mod[Ceiling[height], 2] == 1, 1, 2]]]];
     #), {i, length}];
   #wd);

import numpy as np
from numpy.linalg import inv
from numpy.polynomial import Polynomial as P

# Return a string in X and Y representing the r/s Farey word.
def farey_word(rx,sx):
    lookup_table=[['x','X'],['Y','y']]
    r = int(rx/np.gcd(rx,sx))
    s = int(sx/np.gcd(rx,sx))
    length = 2*s
    def height(i):
        h = i*r/s
        h = h+1/2 if np.ceil(h)==h else h
        return int(np.ceil(h))
    return [ lookup_table[i%2][height(i)%2]  for i in range(1,length+1) ]

# Compute the Farey matrix for slope r/s and matrices [[alpha,1],[0,alpha^-1]], [[beta,0],[mu,beta^-1]].
def farey_matrix(r,s,mu,alpha,beta):
    word = farey_word(r,s)
    X = np.array([[alpha,1],[0,alpha**(-1)]])
    Y = np.array([[beta,0],[mu,beta**(-1)]])
    table = {'X': X,
            'Y': Y,
            'x': inv(X),
            'y': inv(Y)}
    product = np.identity(2)
    for letter in word:
        product = np.matmul(product,table[letter])

    return product

# Return the forwards Farey neighbour of p/q with maximum denominator q.
# Based on box 28 of "Indra's Pearls"
farey_next_cache = {}
def farey_next(p,q):
    if (p,q) in farey_next_cache:
      return farey_next_cache[(p,q)]
    denom = q
    p1,q1 = 0,1
    p2,q2 = 1,0
    r,s = p,q
    sign = -1
    while s != 0:
        a = np.floor(r/s)
        r_new = int(s)
        s_new = int(r - a*s)
        r = int(r_new/np.gcd(r_new,s_new))
        s = int(s_new/np.gcd(r_new,s_new))
        p2_new = int(a*p2 + p1)
        q2_new = int(a*q2 + q1)
        p1,q1 = p2,q2
        p2 = int(p2_new/np.gcd(p2_new,q2_new))
        q2 = int(q2_new/np.gcd(p2_new,q2_new))
        sign = -sign
    k = np.floor((denom - sign*q1)/denom)

    a = int(k*p + sign*p1)
    b = int(k*q + sign*q1)
    u, v = int(a/np.gcd(a,b)), int(b/np.gcd(a,b))
    farey_next_cache[(p,q)] = (u,v)
    return (u,v)

def farey_neighbours(p,q):
    r1,s1 = farey_next(p,q)
    r2 = p - r1
    s2 = q - s1
    return (r1,s1),(r2,s2)

farey_coefficients_fast_parabolic_cache = {}
def farey_coefficients_fast_parabolic(r,s):
  if (r,s) in farey_coefficients_fast_parabolic_cache:
      return farey_coefficients_fast_parabolic_cache[(r,s)]
  if r == 0 and s == 1:
      return P([2,-1])
  if r == 1 and s == 1:
      return P([2,1])
  if r == 1 and s == 2:
      return P([2,0,1])

  (p1,q1),(p2,q2) = farey_neighbours(r,s)

  p = 8 - farey_coefficients_fast_parabolic(p1,q1)*farey_coefficients_fast_parabolic(p2,q2) - farey_coefficients_fast_parabolic(np.abs(p1-p2),np.abs(q1-q2))
  farey_coefficients_fast_parabolic_cache[(r,s)] = p
  return p

# Coefficients for the r/s Farey polynomial corresponding to the (a,b)-elliptic group
def farey_coefficients_fast_elliptic(r,s,a,b):
  farey_coefficients_fast_elliptic_cache = {}
  if (r,s,a,b) in farey_coefficients_fast_elliptic_cache:
      return farey_coefficients_fast_elliptic_cache[(r,s,a,b)]
  alpha = np.exp(2j*np.pi/a)
  beta = np.exp(2j*np.pi/b)

  if r == 0 and s == 1:
      return P([alpha/beta+beta/alpha,-1])
  if r == 1 and s == 1:
      return P([alpha*beta+1/(alpha*beta),1])
  if r == 1 and s == 2:
      return P([2,1/(alpha*beta)+alpha*beta-alpha/beta-beta/alpha,1])

  (p1,q1),(p2,q2) = farey_neighbours(r,s)

  if ((q1 + q2) % 2) == 0:
    p = (4+1/alpha**2+alpha**2 + 1/beta**2 + beta**2) - farey_coefficients_fast_elliptic(p1,q1,a,b)*farey_coefficients_fast_elliptic(p2,q2,a,b) - farey_coefficients_fast_elliptic(np.abs(p1-p2),np.abs(q1-q2),a,b)
  else:
    p = 2*(alpha/beta + beta/alpha + 1/(alpha*beta) + alpha*beta) - farey_coefficients_fast_elliptic(p1,q1,a,b)*farey_coefficients_fast_elliptic(p2,q2,a,b) - farey_coefficients_fast_elliptic(np.abs(p1-p2),np.abs(q1-q2),a,b)
  farey_coefficients_fast_elliptic_cache[(r,s,a,b)] = p
  return p
