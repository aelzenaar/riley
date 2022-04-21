import mpmath as mp
import farey
import riley
import kleinian
import datashader as ds
import datashader.transfer_functions as tf
import pandas
import dask
import dask.dataframe as dd
from dask.delayed import delayed
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from datashader.mpl_ext import dsshow, alpha_colormap

from PyQt5.QtWidgets import QLabel
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class LimitSetView(QLabel):
    selected_changed = pyqtSignal(complex, name='selectedPointChanged')
    def __init__(self, parent=None):
        super(LimitSetView, self).__init__(parent)
        self.resize(640,480)
        self.setMinimumSize(640, 480);

    def redrawLimitSet(self,pOrder,qOrder,mu):
        per_batch = 100
        batches = 10
        depth = 3
        alpha = mp.exp(2j*mp.pi/pOrder)
        beta = mp.exp(2j*mp.pi/qOrder)
        X = farey.generator('X',alpha,beta,mu)
        Y = farey.generator('Y',alpha,beta,mu)
        seeds = farey.fixed_points(0,1,mu,alpha,beta)\
              + farey.fixed_points(1,1,mu,alpha,beta)\
              + farey.fixed_points(1,2,mu,alpha,beta)

        seeds = [seeds[0]]

        print("Found fixed points.",flush=True)

        dask.config.set(scheduler='single-threaded') # We are already running in parallel in each batch; if we don't do this then Dask launches
                                                    # many copies of one_batch() and we get killed by the OOM killer.
        def one_batch(batch):
            #print(f"Running batch {batch+1}/{batches}")
            ls = kleinian.limit_set_markov([X,Y],seeds,depth,per_batch)
            df = pandas.DataFrame(data=[(float(mp.re(point[0])), float(mp.im(point[0])), point[1]) for point in ls], columns=['x','y','colour'], copy=False)
            df['colour']=df['colour'].astype("category")
            return df

        dfs = [delayed(one_batch)(batch) for batch in range(batches)]
        df = dd.from_delayed(dfs)

        fig = plt.figure(figsize=(8, 6), dpi=160)
        canvas = FigureCanvas(fig)
        ax = plt.subplot(111)
        artist = dsshow(df.compute(), ds.Point('x', 'y'), ds.count_cat('colour'), vmax=1000, aspect='equal',ax=ax)
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])

        faces = [plt.Circle((mp.re(1/mu), mp.im(1/mu)), 1/abs(mu), color='r',fill=False),
                plt.Circle((-mp.re(1/mu), -mp.im(1/mu)), 1/abs(mu), color='r',fill=False),
                plt.Rectangle((-1/2,-2),0.01,4),
                plt.Rectangle((1/2,-2),0.01,4)]
        for f in faces:
            ax.add_patch(f)

        canvas.draw()
        width, height = canvas.get_width_height()[0], canvas.get_width_height()[1]
        self.resize(width,height)
        self.image = QImage(canvas.buffer_rgba(), width, height, QImage.Format_ARGB32)
        print((width,height))


    def resizeEvent(self,e):
        super().resizeEvent(e)
        self.setPixmap(QPixmap(self.image.scaled(self.size())))

