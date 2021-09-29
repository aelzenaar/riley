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

### Python library
There are three files containing general Python code which can be called in the Python interpreter or used in Python scripts.

 * [kleinian.py](kleinian.py) -- methods for general Kleinian groups (e.g. limit set calculations)
 * [farey.py](farey.py) -- methods for working with Farey words and polynomials
 * [riley.py](riley.py) -- methods for working with the Riley slices


## References
<a id="KS94">[KS94]</a>
Linda Keen and Caroline Series. “The Riley slice of Schottky space”. In: *Proceedings of the London Mathematics Society* 3.1 (69 1994), pp. 72–90.

<a id="KS98">[KS98]</a>
Yohei Komori and Caroline Series. “The Riley slice revisited”. In: *The Epstein birthday schrift*. Vol. 1. Geometry and Topology Monographs. 1998, pp. 303–316.
