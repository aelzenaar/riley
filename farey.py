""" Methods to compute with Farey words and their traces (`Farey polynomials').

    The Farey words were introduced by Linda Keen and Caroline Series in the 1994 paper "The Riley slice of Schottky space" (Proc. LMS (3)69 p.72-90); their theory
    enables the explicit computation of a foliation of the Riley slice. In this file are found methods to compute Farey words and polynomials. The methods which deal
    with the Riley slice and its geometry are found in riley.py.
"""

from mpmath import mp
mp.dps = 1000

import math

from numpy.polynomial import Polynomial as P
import kleinian
from functools import cache


@cache
def generator(letter,alpha,beta,mu):
    """ Return the generator of the group corresponding to the given letter, with the correct values substituted in.

        For reference, the generators are X = [[alpha,1],[0,alpha^(-1)]] and Y = [[beta,0],[mu,beta^(-1)]] with respective inverses x, y.

        Arguments:
          letter - one of X, Y, x, y
          alpha - upper-left entry of X
          beta - upper-left entry of Y
          mu - lower-right entry of Y
    """
    X = mp.matrix([[alpha,1],[0,alpha**(-1)]])
    Y = mp.matrix([[beta,0],[mu,beta**(-1)]])
    table = {'X': X,
            'Y': Y,
            'x': kleinian._fast_inv(X),
            'y': kleinian._fast_inv(Y)}

    try:
        return table[letter]
    except IndexError:
        raise IndexError("Unknown generator for the group")

# Return a string in X and Y representing the r/s Farey word.
@cache
def word(r,s):
    """ Compute the Farey word of slope r/s using the cutting sequence definition.

        Arguments:
          r, s -- coprime integers such that r/s is the slope of the desired Farey word

        Returns:
          A list consisting of single-character strings representing generators of the group and their inverses (as defined in generator() above)
    """

    if math.gcd(r,s) != 1:
        raise ValueError("Arguments to farey_word should be coprime integers.")

    lookup_table=[['x','X'],['Y','y']]
    length = 2*s
    def height(i):
        h = i*r/s
        h = h+1/2 if mp.ceil(h)==h else h
        return int(mp.ceil(h))
    return [ lookup_table[i%2][height(i)%2]  for i in range(1,length+1) ]

def matrix(r,s,mu,alpha,beta):
    """ Compute the Farey word of slope r/s in matrix form.

        Arguments:
          r, s --- coprime integers representing the slope r/s
          mu, alpha, beta --- parameters of the group as defined in generator() above
    """
    w = word(r,s)
    product = mp.eye(2)
    for letter in w:
        product = product * generator(letter,alpha,beta,mu)

    return product

def fixed_points(r,s,mu,alpha,beta):
    """ Compute the fixed points of the Farey word of slope r/s.

        Arguments:
          r, s --- coprime integers representing the slope r/s
          mu, alpha, beta --- parameters of the group as defined in generator() above
    """


    m = matrix(r,s,mu,alpha,beta)
    surd = mp.sqrt((m[1,1] - m[0,0])**2 - 4*m[0,1]*m[1,0])
    trans = m[0,0] - m[1,1]

    return [(trans+surd)/(2*m[1,0]),(trans-surd)/(2*m[1,0])]

@cache
def next_neighbour(p,q):
    """ Return the larger Farey neighbour of p/q in the Farey sequence of denominator q.

        This method is based on the pseudocode found in Box 28 of the book "Indra's Pearls" by David Mumford, Caroline Series, and David Wright (Cambridge Uni. Press, 2002).

        Arguments:
          p,q -- Coprime integers representing the fraction p/q in the interval [0,1].
    """

    if math.gcd(p,q) != 1:
        raise ValueError("Arguments to farey_next should be coprime integers.")

    denom = q
    p1,q1 = 0,1
    p2,q2 = 1,0
    r,s = p,q
    sign = -1
    while s != 0:
        a = mp.floor(r/s)
        r_new = int(s)
        s_new = int(r - a*s)
        r = int(r_new/math.gcd(r_new,s_new))
        s = int(s_new/math.gcd(r_new,s_new))
        p2_new = int(a*p2 + p1)
        q2_new = int(a*q2 + q1)
        p1,q1 = p2,q2
        p2 = int(p2_new/math.gcd(p2_new,q2_new))
        q2 = int(q2_new/math.gcd(p2_new,q2_new))
        sign = -sign
    k = mp.floor((denom - sign*q1)/denom)

    a = int(k*p + sign*p1)
    b = int(k*q + sign*q1)
    u, v = int(a/math.gcd(a,b)), int(b/math.gcd(a,b))
    return (u,v)

@cache
def neighbours(p,q):
    """ Compute the two Farey neighbours of p/q.

        Arguments:
          p,q -- Coprime integers representing the fraction p/q in the interval [0,1].
    """
    r1,s1 = next_neighbour(p,q)
    r2 = p - r1
    s2 = q - s1
    if r1/s1 < r2/s2:
        return (r1,s1),(r2,s2)
    else:
        return (r2,s2),(r1,s1)

@cache
def _even_const(alpha,beta):
    return 4+2*mp.re(alpha**2)+2*mp.re(beta**2)

@cache
def _odd_const(alpha,beta):
    return 4*(mp.re(alpha/beta) + mp.re(alpha*beta))

@cache
def polynomial_coefficients_fast(r,s,alpha,beta,_=None):
    """ Return the coefficients of the Farey polynomial of slope r/s.

        The method used is the recursion algorithm.

        Arguments:
          r,s -- coprime integers representing the slope of the desired polynomial
          alpha, beta -- parameters of the group
    """

    if r == 0 and s == 1:
        return P([2*mp.re(alpha/beta),-1])
    if r == 1 and s == 1:
        return P([2*mp.re(alpha*beta),1])
    if r == 1 and s == 2:
        return P([2,-4*mp.im(alpha)*mp.im(beta),1])

    (p1,q1),(p2,q2) = neighbours(r,s)
    konstant = _even_const(alpha,beta) if ((q1 + q2) % 2) == 0 else _odd_const(alpha,beta)

    p =  konstant-(polynomial_coefficients_fast(p1,q1,alpha,beta)*polynomial_coefficients_fast(p2,q2,alpha,beta) + polynomial_coefficients_fast(mp.fabs(p1-p2),mp.fabs(q1-q2),alpha,beta))


    return p

@cache
def polynomial_coefficients_reduced(r,s):
    """ Return the coefficients of the reduced Farey polynomial (\Phi^{\infty,\infty}-2) of slope r/s.

        (This is the polynomial Q_r/s of our 2021 preprint.)

        The method used is the recursion algorithm.

        Arguments:
          r,s -- coprime integers representing the slope of the desired polynomial
    """

    if r == 0 and s == 1:
        return P([0,-1])
    if r == 1 and s == 1:
        return P([0,1])
    if r == 1 and s == 2:
        return P([0,0,1])

    (p1,q1),(p2,q2) = neighbours(r,s)
    p = - polynomial_coefficients_reduced(int(abs(p1-p2)),int(abs(q1-q2)))\
        - polynomial_coefficients_reduced(p1,q1)*polynomial_coefficients_reduced(p2,q2)\
        - 2 * (polynomial_coefficients_reduced(p1,q1) + polynomial_coefficients_reduced(p2,q2))

    return p
