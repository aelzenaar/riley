# Riley slice computational package

This is a computational package for dealing with the Riley slice of Schottky space [[KS94](#KS94), [KS98](#KS98)]  and its elliptic
generalisations, along with the associated combinatorial group theory and geometry. In order to generate Farey polynomials, we use the
results obtained in our paper [[EMS21b](#EMS21b)]. We also have some further recent results on the Riley slice which are of interest
from a computational point of view [[EMS21a](#EMS21a)]; in a later version of this software we will incorporate some of the insights
from this paper. The program [graphical_limits.py](graphical_limits.py) was inspired by the [schottky](https://github.com/dannycalegari/schottky)
software written by Danny Calegari and Alden Walker.

![The parabolic Riley slice](riley_slice.png?raw=true "The parabolic Riley slice")

* [Authorship](#authorship)
* [Background material and related work](#background-material-and-related-work)
* [The software included here](#the-software-included-here)
  + [Fun programs to run](#fun-programs-to-run)
  + [Python library](#python-library)
* [Dependencies](#dependencies)
* [References](#references)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Authorship
This software was written by [Alex Elzenaar](https://aelzenaar.github.io) under the supervision of Gaven Martin (NZ Institute of Advanced Study, Massey University) and Jeroen Schillewaert (The University of Auckland).


## Background material and related work

For background material in Kleinian groups which underpins the mathematics we study using this software, see [[B83](#B83),[M87](#M87)]. For a more
practical introduction to some of the computational geometry and some very nice pictures, see [[MSW02](#MSW02)] and [its associated website](http://klein.math.okstate.edu/IndrasPearls/).

## The software included here

### Fun programs to run

 * [graphical_limits.py](graphical_limits.py) -- this is a piece of GUI software which displays the various Riley slices and allows
   the user to click around to view the limit set of a particular group.
 * [limit_plotter.py](limit_plotter.py) -- a script to plot limit sets of general finitely-generated Kleinian groups (in colour!).
    Example output: [the parabolic 1/2-cusp limit set](cusp12_3.png) (max word length 15, 50000 points)
 * [slice_plotter.py](slice_plotter.py) -- a script to plot Riley slices without all the extra machinery in graphical_limits.py.
 * [farey_graph.py](farey_graph.py) -- a script to plot the Farey addition graph with the butterfly colouring. Example
    output: [the addition graph](farey_graph_coloured.png), [the Stern-Brocot tree](farey_graph_tree.png)
 * [cusps.py](cusps.py) -- a script to plot the limit set at a given cusp point. If you have a powerful computer and want much better pictures, you probably
    want to use [cusps_shaded.py](cusps_shaded.py) instead which attempts to use [datashader](https://datashader.org/).
 * [generate_polynomials.py](generate_polynomials.py) -- generate a file of Mathematica code containing lots of Farey polynomials.
 * [limit_set_with_circles.py](limit_set_with_circles.py) -- plot limit sets together with the isometric circles.

## Future features? Some easy, some (very) hard
 * Draw the associated surface for a point (somehow) together with the corresponding foliation
 * Plot pleating ray of given slope
 * Check computationally if a point is in or out
 * Zooming
 * Trace curve & give animation
 * Compute Teichmuller distance, draw Teichmuller geodesics. Perhaps this is best done by studying the associated foliations.
 * Plot 2-bridge knot groups and other features of interest in the exterior. See for instance [[ASWY07](#ASWY07)].

### Python library
There are three files containing general Python code which can be called in the Python interpreter or used in Python scripts.

 * [kleinian.py](kleinian.py) -- methods for general Kleinian groups (e.g. limit set calculations)
 * [farey.py](farey.py) -- methods for working with Farey words and polynomials
 * [riley.py](riley.py) -- methods for working with the Riley slices


## Dependencies
Disclaimer: this is [academic software](https://academia.stackexchange.com/questions/37370/should-i-share-my-horrible-software) so may require some fiddling to get it to work on your machine.

 * Python 3 (tested on 3.9.7)
 * scipy (all)
 * [mpmath](https://mpmath.org/) (all)
 * [mpsolve](https://numpi.dm.unipi.it/software/mpsolve) (optional, falls back to numpy if not installed - but produces much better results)
 * matplotlib (limit_plotter.py, slice_plotter.py)
 * tkinter (graphical_limits.py)
 * pydot, networkx (farey_graph.py)
 * [datashader](https://datashader.org/), [pandas](https://pandas.pydata.org/), and [dask](https://dask.org/) (cusps_shaded.py)

## References
<a id="ASWY07">[ASWY07]</a>
Hirotaka Akiyoshi, Makoto Sakuma, Masaaki Wada, and Yasushi Yamashita. _Punctured torus groups and 2-bridge knot groups I_. Lecture Notes in Mathematics 1909. Springer, 2007

<a id="B88">[B83]</a>
Alan F. Beardon. *The geometry of discrete groups*. Graduate Texts in Mathematics 91. Springer-Verlag, 1983.

<a id="EMS21a">[EMS21a]</a>
Alex Elzenaar, Gaven Martin, and Jeroen Schillewaert. “Approximations of the Riley slice”. November 2021. [arXiv:2111.03230](https://arxiv.org/abs/2111.03230) [math.GT].

<a id="EMS21b">[EMS21b]</a>
Alex Elzenaar, Gaven Martin, and Jeroen Schillewaert. “The combinatorics of Farey words and their traces”. In preparation.

<a id="KS94">[KS94]</a>
Linda Keen and Caroline Series. “The Riley slice of Schottky space”. In: *Proceedings of the London Mathematics Society* 3.1 (69 1994), pp. 72–90.

<a id="KS98">[KS98]</a>
Yohei Komori and Caroline Series. “The Riley slice revisited”. In: *The Epstein birthday schrift*. Vol. 1. Geometry and Topology Monographs. 1998, pp. 303–316.

<a id="M87">[M87]</a>
Bernard Maskit. *Kleinian groups*. Grundlehren der mathematischen Wissenshaften 287. Springer-Verlag, 1987.

<a id="MSW02">[MSW02]</a>
David Mumford, Caroline Series, and David Wright. *Indra’s pearls: The vision of Felix Klein*. Cambridge University Press, 2002.
