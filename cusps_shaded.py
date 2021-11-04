""" A script to plot the limit set of a particular cusp group, with maximal power and efficiency.

    This is only as accurate as the riley.cusp_point() function allows (i.e. make sure you have
    sacrificed the correct number of goats before running the script...).

    The limit point seeds are the fixed points of the 0/1, 1/1, and 1/2 Farey words at the given cusp point.

    Limit set is plotted using datashader and dask --- thus this script will not run into any of the memory
    issues that cusps.py has.

    Options to change:
        p, q -- orders of the elliptic elements (set to mp.inf for the parabolic case)
        r, s -- slope of the desired cusp
        per_batch, batches -- you should be able to run kleinian.limit_set_markov with reps set to `per_batch'; we do this `batches' times (so we don't run out of RAM).
        depth -- maximum word length to compute orbits with

    Output image filename is cusp_{r}_{s}_elliptic_{p}_{q}_shaded.png (in the current directory).
"""


import mpmath as mp
import kleinian
import riley
import farey
import datashader as ds
import datashader.transfer_functions as tf
import pandas
import dask
import dask.dataframe as dd
from dask.delayed import delayed
from datashader.utils import export_image

# Orders of elliptics
p = 3
q = 5

# Cusp slope
r = 1
s = 2

filename = f'cusp_{r}_{s}_elliptic_{p}_{q}_shaded'

per_batch = 100000
batches = 20
depth = 15

mu = riley.cusp_point(p,q,r,s)

alpha = mp.exp(2j*mp.pi/p)
beta = mp.exp(2j*mp.pi/q)
X = mp.matrix([[alpha,1],[0,mp.conj(alpha)]])
Y = mp.matrix([[beta,0],[mu,mp.conj(beta)]])

seeds = farey.fixed_points(0,1,mu,alpha,beta)\
      + farey.fixed_points(1,1,mu,alpha,beta)\
      + farey.fixed_points(1,2,mu,alpha,beta)

print("Found fixed points.",flush=True)

dask.config.set(scheduler='single-threaded') # We are already running in parallel in each batch; if we don't do this then Dask launches
                                             # many copies of one_batch() and we get killed by the OOM killer.
def one_batch(batch):
  print(f"Running batch {batch+1}/{batches}")
  ls = kleinian.limit_set_markov([X,Y],seeds,depth,per_batch)
  df = pandas.DataFrame(data=[(float(mp.re(point[0])), float(mp.im(point[0])), point[1]) for point in ls], columns=['x','y','colour'], copy=False)
  df['colour']=df['colour'].astype("category")
  return df

dfs = [delayed(one_batch)(batch) for batch in range(batches)]
df = dd.from_delayed(dfs)
cvs = ds.Canvas(plot_width=4000,plot_height=4000,x_range=(-4,4), y_range=(-4,4), x_axis_type='linear', y_axis_type='linear')

agg = cvs.points(df,'x','y')
img = tf.shade(agg, cmap="black", min_alpha=0)

#aggc = cvs.points(df,'x','y',ds.by('colour', ds.count()))
#colours = {-2: 'red', -1:'blue', 1:'green', 2:'purple'}
#img =  tf.shade(aggc)

export_image(img, filename, background="white", export_path=".")
