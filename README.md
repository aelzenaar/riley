# Riley slice computational package

This is a computational package for dealing with the Riley slice of Schottky space [[KS94](#KS94), [KS98](#KS98)]  and its elliptic
generalisations, along with the associated combinatorial group theory and geometry.

![The parabolic Riley slice](riley_slice.png?raw=true "The parabolic Riley slice")


## Background material

## The software included here

### Fun programs to run

 * [graphical_limits.py](graphical_limits.py) -- this is a piece of GUI software which displays the various Riley slices and allows
   the user to click around to view the limit set of a particular group.
 * [limit_plotter.py](limit_plotter.py) -- a script to plot limit sets of general finitely-generated Kleinian groups (in colour!).
    Example output: [the parabolic 1/2-cusp limit set](cusp12_3.png) (max word length 15, 50000 points)
 * [slice_plotter.py](slice_plotter.py) -- a script to plot Riley slices without all the extra machinery in graphical_limits.py.
 * [farey_graph.py](farey_graph.py) -- a script to plot the Farey addition graph with the butterfly colouring. Example
    output: [the addition graph](farey_graph_coloured.png), [the Stern-Brocot tree](farey_graph_tree.png)
 * [cusps.py](cusps.py) -- a script to plot the limit set at a given cusp point. See also [cusps_shaded.py](cusps_shaded.py) which attempts to use [datashader](https://datashader.org/).
 * [generate_polynomials.py](generate_polynomials.py) -- generate a file of Mathematica code containing lots of Farey polynomials.

### Python library
There are three files containing general Python code which can be called in the Python interpreter or used in Python scripts.

 * [kleinian.py](kleinian.py) -- methods for general Kleinian groups (e.g. limit set calculations)
 * [farey.py](farey.py) -- methods for working with Farey words and polynomials
 * [riley.py](riley.py) -- methods for working with the Riley slices


## Dependencies

 * Python 3 (tested on 3.9.7)
 * scipy (all)
 * [mpmath](https://mpmath.org/) (all)
 * matplotlib (limit_plotter.py, slice_plotter.py)
 * tkinter (graphical_limits.py)
 * pydot, networkx (farey_graph.py)
 * [mpsolve](https://numpi.dm.unipi.it/software/mpsolve) (optional, falls back to numpy if not installed - but produces much better results)

## References
<a id="KS94">[KS94]</a>
Linda Keen and Caroline Series. “The Riley slice of Schottky space”. In: *Proceedings of the London Mathematics Society* 3.1 (69 1994), pp. 72–90.

<a id="KS98">[KS98]</a>
Yohei Komori and Caroline Series. “The Riley slice revisited”. In: *The Epstein birthday schrift*. Vol. 1. Geometry and Topology Monographs. 1998, pp. 303–316.
