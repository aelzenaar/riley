import numpy as np
import kleinian
import matplotlib.pyplot as plt
import riley
import farey
import itertools
import datashader as ds, pandas as pd, colorcet as cc


# Orders of elliptics
p = 3
q = 5

# Cusp slope
r = 1
s = 2

reps = 10000
depth = 10




mu = riley.cusp_point(p,q,r,s)

alpha = np.csingle(np.exp(2j*np.pi/p))
beta = np.csingle(np.exp(2j*np.pi/q))
X = np.array([[alpha,1],[0,np.conj(alpha)]])
Y = np.array([[beta,0],[mu,np.conj(beta)]])

seeds = farey.fixed_points(0,1,mu,alpha,beta)\
      #+ farey.fixed_points(1,1,mu,alpha,beta)\
      #+ farey.fixed_points(1,2,mu,alpha,beta)

ls = kleinian.limit_set_markov([X,Y],np.array(seeds),depth,True,reps)

df = pd.DataFrame(data=list(itertools.chain(*[[(np.real(point), np.imag(point), points[1]) for point in points[0]] for points in ls])), columns=['x','y','colour'])
df['colour']=df['colour'].astype("category")
cvs = ds.Canvas(x_range=(-8,8), y_range=(-8,8), x_axis_type='linear', y_axis_type='linear')
aggc = cvs.points(df,'x','y',ds.by('colour', ds.count()))

colours = {-2: 'red', -1:'blue', 1:'green', 2:'purple'}
img = ds.tf.Image(ds.tf.shade(aggc, color_key=colours))

plt.imshow(img)
plt.savefig('cusp12_elliptic33_shader.png',dpi=2000)
