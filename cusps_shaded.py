import numpy as np
import kleinian
import matplotlib.pyplot as plt
import riley
import farey
import datashader as ds
import datashader.transfer_functions as tf
import pandas

# Orders of elliptics
p = 3
q = 5

# Cusp slope
r = 1
s = 2

reps = 1000000
depth = 15



mu = np.cdouble(riley.cusp_point(p,q,r,s))

alpha = np.cdouble(np.exp(2j*np.pi/p))
beta = np.cdouble(np.exp(2j*np.pi/q))
X = np.array([[alpha,1],[0,np.conj(alpha)]])
Y = np.array([[beta,0],[mu,np.conj(beta)]])

seeds = farey.fixed_points(0,1,mu,alpha,beta)\
      + farey.fixed_points(1,1,mu,alpha,beta)\
      + farey.fixed_points(1,2,mu,alpha,beta)

print("Found fixed points.",flush=True)

df = pandas.DataFrame(data=[(np.real(point[0]), np.imag(point[0]), point[1]) for point in kleinian.limit_set_markov([X,Y],np.cdouble(seeds),depth,reps)], columns=['x','y','colour'], copy=False)
df['colour']=df['colour'].astype("category")
cvs = ds.Canvas(plot_width=2000,plot_height=2000,x_range=(-4,4), y_range=(-4,4), x_axis_type='linear', y_axis_type='linear')
aggc = cvs.points(df,'x','y',ds.by('colour', ds.count()))

#colours = {-2: 'red', -1:'blue', 1:'green', 2:'purple'}
img =  tf.shade(aggc)

plt.imshow(img)
plt.savefig('cusp12_elliptic35_shader2.png',dpi=2000)
