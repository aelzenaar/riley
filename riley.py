""" Methods to compute globally and locally with the Riley slice.

    In this file are found methods to compute with the moduli space known as the `Riley slice'.
"""

from ast import literal_eval
import farey
from numpy.polynomial import Polynomial as P
import numpy as np
import scipy.optimize

try:
    import mpsolve
    mpsolve_avail = True
except ImportError:
    mpsolve_avail = False

def riley_slice(a, b, max_denom, use_mpsolve=mpsolve_avail, max_iter=100, tol=1e-2):
    """ Return an accurate approximation to the Riley slice.

        There are two possible solvers: the mpsolve solver (see https://numpi.dm.unipi.it/software/mpsolve), or the solver built
        in to scipy followed by an attempt to improve the results using Newton's algorithm.

        Arguments:
          a, b -- the order of the cone points represented by X and Y respectively. Use np.inf for the parabolic case (or 1, since exp(2*pi*i/1) = exp(0) = 1).
          max_denom -- the maximum denominator Farey polynomials to compute
          use_mpsolve - use the mpsolve solver if True or scipy if False (default True if mpsolve is installed, else False)

        The following arguments are only used if use_mpsolve is False:
          max_iter -- maximum number of iterations for Newton's algorithm (default 100)
          tol -- allowable error of each point (default 1e-2)
    """

    if use_mpsolve:
        if not mpsolve_avail:
            raise RuntimeError('mpsolve selected even though it seems not to be installed')
        mpsolve_ctx = mpsolve.Context()

    alpha = 1 if a == np.inf else np.exp(2j*np.pi/a)
    beta = 1 if b == np.inf else np.exp(2j*np.pi/b)

    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if np.gcd(p,q) == 1:
          poly = farey.polynomial_coefficients_fast(p,q,alpha,beta) + 2
          if use_mpsolve:
              mpsolve_poly = mpsolve.MonomialPoly(mpsolve_ctx, q)

              for d in range(0,q+1):
                  if alpha == 1 and beta == 1:
                      mpsolve_poly.set_coefficient(d, int(poly.coef[d]))
                  else:
                      mpsolve_poly.set_coefficient(d, np.real(poly.coef[d]), np.imag(poly.coef[d]))
              points.extend(mpsolve_ctx.solve(mpsolve_poly))
          else:
              roots_bad = poly.roots()
              try:
                roots = [scipy.optimize.newton(poly, root, poly.deriv(),fprime2=poly.deriv(2),maxiter=max_iter,tol=tol) for root in roots_bad]
              except RuntimeError as e:
                print(str(e))
                print(f'Newton failed to converge at {p}/{q}')
                return points
              points.extend(roots)
    return points

def riley_centre(a,b):
    """ Return an approximation to the centre of symmetry of the Riley slice.

        Arguments:
          a,b -- orders of X and Y respectively
    """
    alpha = 1 if a == np.inf else np.exp(2j*np.pi/a)
    beta = 1 if b == np.inf else np.exp(2j*np.pi/b)

    poly = farey.polynomial_coefficients_fast(1,1,alpha,beta) + 2
    roots = poly.roots()
    return (4 + roots[0])/2

def cusp_point(a,b,p,q, use_mpsolve=mpsolve_avail):
    """ Return an approximation to the p/q-cusp point on the Riley slice boundary.

        Let Phi_{p/q} be the p/q-Farey polynomial; the p/q-pleating ray is then the connected
        component of Phi_{p/q}^{-1}((-\infty,-2)) with asymptotic slope pi*p/q, and the p/q-cusp
        is the point on this branch which is the inverse image of -2.

        Arguments:
          a,b -- orders of X and Y respectively
          p,q -- integers representing the slope of the desired cusp; p and q must be coprime and nonnegative and p/q must lie in [0,2).
          use_mpsolve - use the mpsolve solver if True or scipy if False (default True if mpsolve is installed, else False)
    """

    if p/q > 1:
        return np.conj(cusp_point(a,b,2*q-p,q))
    alpha = 1 if a == np.inf else np.exp(2j*np.pi/a)
    beta = 1 if b == np.inf else np.exp(2j*np.pi/b)

    poly = farey.polynomial_coefficients_fast(p,q,alpha,beta) + 2

    if use_mpsolve:
        if not mpsolve_avail:
            raise RuntimeError('mpsolve selected even though it seems not to be installed')
        mpsolve_ctx = mpsolve.Context()
        mpsolve_poly = mpsolve.MonomialPoly(mpsolve_ctx, q)
        for d in range(0,q+1):
            if alpha == 1 and beta == 1:
                mpsolve_poly.set_coefficient(d, int(poly.coef[d]))
            else:
                mpsolve_poly.set_coefficient(d, np.real(poly.coef[d]), np.imag(poly.coef[d]))
        roots = mpsolve_ctx.solve(mpsolve_poly)
    else:
        roots = poly.roots()

    if q == 1:
        return roots[0]

    centre = riley_centre(a,b)

    if q > 1:
        (r1,s1),(r2,s2) = farey.neighbours(p,q)
        left_cusp_angle = np.angle(cusp_point(a,b,r1,s1) - centre)
        right_cusp_angle = np.angle(cusp_point(a,b,r2,s2) - centre)

        desired_cusp = None
        for root in roots:
            argument = np.angle(root - centre)
            if left_cusp_angle < argument and argument < right_cusp_angle:
                desired_cusp = root
                break
        if desired_cusp == None:
            raise RuntimeError('failed to find cusp')
        return desired_cusp
    else:
        raise ValueError('q < 1?')
