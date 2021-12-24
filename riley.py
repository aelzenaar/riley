""" Methods to compute globally and locally with the Riley slice.

    In this file are found methods to compute with the moduli space known as the `Riley slice'.
"""

import farey
from mpmath import mp
mp.dps = 100
import math
import scipy.optimize
from numpy.polynomial import Polynomial as P

try:
    import mpsolve
    mpsolve_avail = True
except ImportError:
    mpsolve_avail = False

try:
    import sympy
    sympy_avail = True
except ImportError:
    sympy_avail = False

try:
    import matlab.engine
    matlab_avail = True
    matlab_eng = None
except ImportError:
    matlab_avail = False

def poly_solve(poly, solver='mpsolve' if mpsolve_avail else 'scipy', max_iter=100, tol=1e-2, try_int=False):
    """ Solve a polynomial numerically.

        There are three possible solvers: the mpsolve solver (see https://numpi.dm.unipi.it/software/mpsolve), the solver built
        in to scipy followed by an attempt to improve the results using Newton's algorithm, or sympy.

        Arguments:
          poly --- a scipy polynomial
          solver -- one of 'mpsolve', 'scipy', 'sympy'

        The following arguments are only used by the scipy solver:
          max_iter -- maximum number of iterations for Newton's algorithm (default 100)
          tol -- allowable error of each point (default 1e-2)

        The following arguments are only used by the mpsolve solver:
          try_int -- assume that the coefficients of poly are integral (default False)
    """

    q = poly.degree()

    if solver == 'mpsolve':
        if not mpsolve_avail:
            raise RuntimeError('mpsolve selected even though it seems not to be installed')
        mpsolve_ctx = mpsolve.Context()
        mpsolve_poly = mpsolve.MonomialPoly(mpsolve_ctx, q)

        for d in range(0,q+1):
            if try_int:
                mpsolve_poly.set_coefficient(d, int(poly.coef[d]))
            else:
                mpsolve_poly.set_coefficient(d, float(mp.re(poly.coef[d])))
        return mpsolve_ctx.solve(mpsolve_poly)

    elif solver == 'sympy':
        if not sympy_avail:
            raise RuntimeError('sympy selected even though it seems not to be installed')
        ZZ = sympy.Symbol('Z')
        sympy_poly = 0
        for d in range(0,q+1):
            sympy_poly = sympy_poly + poly.coef[d] * ZZ**d
        return [sympy.N(root) for root in sympy.solve(sympy_poly,ZZ)]

    elif solver == 'scipy':
        poly = P(list(map(float,poly.coef))) # Need to cast to float in case we have mpmath but not mpsolve
        roots_bad = poly.roots()
        return [scipy.optimize.newton(poly, root, poly.deriv(),fprime2=poly.deriv(2),maxiter=max_iter,tol=tol) for root in roots_bad]

    elif solver == 'matlab':
        if not matlab_avail:
            raise RuntimeError('matlab selected even though it seems not to be installed')
        global matlab_eng
        if matlab_eng == None:
            matlab_eng = matlab.engine.start_matlab()

        roots = matlab_eng.roots(matlab.double([complex(x) for x in poly.coef],is_complex=True))
        try:
            return list(roots)
        except TypeError:
            return [roots]

    raise RuntimeError(f'unknown solver {solver}')

def riley_slice(a, b, max_denom, solver='mpsolve' if mpsolve_avail else 'scipy', **kwargs):
    """ Return an accurate approximation to the Riley slice.

        There are three possible solvers: the mpsolve solver (see https://numpi.dm.unipi.it/software/mpsolve), the solver built
        in to scipy followed by an attempt to improve the results using Newton's algorithm, or sympy.

        Arguments:
          a, b -- the order of the cone points represented by X and Y respectively. Use mp.inf for the parabolic case (or 1, since exp(2*pi*i/1) = exp(0) = 1).
          max_denom -- the maximum denominator Farey polynomials to compute
          solver -- one of 'mpsolve', 'scipy', 'sympy'

        Further keyword arguments are passed directly to poly_solve(), i.e. tol and max_iter for scipy.
    """

    alpha = 1 if a == mp.inf else mp.exp(2j*mp.pi/a)
    beta = 1 if b == mp.inf else mp.exp(2j*mp.pi/b)

    points = []
    for q in range(1,max_denom+1):
      for p in range(1,q+1):
        if math.gcd(p,q) == 1:
          poly = farey.polynomial_coefficients_fast(p, q, alpha, beta) + 2
          try_int = True if alpha == 1 and beta == 1 else False
          try:
            points.extend(poly_solve(poly, solver, try_int=try_int, **kwargs))
          except Exception as e:
            raise RuntimeError(f'p = {p} q={q}') from e
    return points

def riley_centre(a,b):
    """ Return an approximation to the centre of symmetry of the Riley slice.

        Arguments:
          a,b -- orders of X and Y respectively
    """
    alpha = 1 if a == mp.inf else mp.exp(2j*mp.pi/a)
    beta = 1 if b == mp.inf else mp.exp(2j*mp.pi/b)

    poly = farey.polynomial_coefficients_fast(1, 1, alpha, beta, int if (alpha == 1 and beta == 1) else mp.mpf) + 2
    roots = poly.roots()
    return (4 + roots[0])/2

def cusp_point(a, b, p, q, solver='mpsolve' if mpsolve_avail else 'scipy', **kwargs):
    """ Return an approximation to the p/q-cusp point on the Riley slice boundary.

        Let Phi_{p/q} be the p/q-Farey polynomial; the p/q-pleating ray is then the connected
        component of Phi_{p/q}^{-1}((-\\infty,-2)) with asymptotic slope pi*p/q, and the p/q-cusp
        is the point on this branch which is the inverse image of -2.

        Arguments:
          a,b -- orders of X and Y respectively
          p,q -- integers representing the slope of the desired cusp; p and q must be coprime and nonnegative and p/q must lie in [0,2).
          solver -- one of 'mpsolve', 'scipy', 'sympy'

        Further keyword arguments are passed directly to poly_solve(), i.e. tol and max_iter for scipy.
    """

    if p/q > 1:
        return mp.conj(cusp_point(a,b,2*q-p,q))
    alpha = 1 if a == mp.inf else mp.exp(2j*mp.pi/a)
    beta = 1 if b == mp.inf else mp.exp(2j*mp.pi/b)

    poly = farey.polynomial_coefficients_fast(p, q, alpha, beta, int if (alpha == 1 and beta == 1) else mp.mpf) + 2
    roots = poly_solve(poly, solver, **kwargs)

    if q == 1:
        return roots[0]

    centre = riley_centre(a,b)

    if q > 1:
        (r1,s1),(r2,s2) = farey.neighbours(p,q)
        left_cusp_angle = mp.arg(cusp_point(a,b,r1,s1) - centre)
        right_cusp_angle = mp.arg(cusp_point(a,b,r2,s2) - centre)

        right_argument_roots = []
        for root in roots:
            argument = mp.arg(root - centre)
            #print(f'{left_cusp_angle} < {argument} < {right_cusp_angle}?')
            if left_cusp_angle < argument and argument < right_cusp_angle:
                right_argument_roots.append(root)

        if right_argument_roots == []:
            raise RuntimeError('failed to find cusp: couldn\'t bound argument')

        return max(right_argument_roots, key=mp.fabs)
    else:
        raise ValueError('q < 1?')
