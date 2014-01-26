# some general tools for plotting of data 

import numpy as np
import matplotlib.pyplot as plt
from bandtools import *
from mytools import *
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap


def mycolplot(filename,numc):
  """make a color plot of a data file in this format:
  x1 y1 z1
  x2 y2 z2
  ...
  numc is the number of columns of the (x,y) grid."""
  dat = readmat(filename)
  x = np.zeros((numc,numc))
  y = np.zeros((numc,numc))
  p = np.zeros((numc,numc))
  for i in range(numc):
    for j in range(numc):
      x[i,j] = dat[i*numc+j,0]
      y[i,j] = dat[i*numc+j,1]
      p[i,j] = dat[i*numc+j,2]
  plt.figure()
  plt.axis('equal')
  plt.xlabel("x ", size='x-large')
  plt.ylabel("y ", size='x-large')
  plt.xticks(size='x-large')
  plt.yticks(size='x-large')
  #plt.xlim(0,210)
  #plt.ylim(0,210)
  plt.pcolor(x,y,p)
  plt.colorbar()
  plt.show()
  plt.savefig("color_plot.png")

def mycolplot1(dat,numr,numc):
  """make a color plot of a data in this format:
  x1 y1 z1
  x2 y2 z2
  ...
  numc is the number of columns of the (x,y) grid.
  numr is the number of rows of the (x,y) grid."""
  x = np.zeros((numr,numc))
  y = np.zeros((numr,numc))
  p = np.zeros((numr,numc))
  for i in range(numr):
    for j in range(numc):
      x[i,j] = dat[i*numc+j,0]
      y[i,j] = dat[i*numc+j,1]
      p[i,j] = dat[i*numc+j,2]
  plt.figure()
  #plt.axis('equal')
  plt.xlabel("x ", size='x-large')
  plt.ylabel("y ", size='x-large')
  plt.xticks(size='x-large')
  plt.yticks(size='x-large')
  #plt.xlim(0,210)
  #plt.ylim(0,210)
  plt.pcolor(x,y,p)
  plt.colorbar()
  plt.show()
  plt.savefig("color_plot.png")


