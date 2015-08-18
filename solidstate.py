import numpy as np
from numpy import linalg as LA
import geometry as geo

class Atom:
  def __init__(self,symbol,pos):
    self.name = symbol
    self.pos = pos
  def x(self):
    return self.pos[0]
  def y(self):
    return self.pos[1]
  def z(self):
    return self.pos[2]

def dist(a1,a2):
  return LA.norm(a1.pos - a2.pos)

#translate a list of Atoms
# t = numpy array
def trans(t,aa):
  ans = []
  for a in aa:
    p1 = t+a.pos
    ans.append(Atom(a.name,p1))
  return ans

#rotate a list of Atoms
# udir = numpy array
def rot(udir,theta,aa):
  ans = []
  for a in aa:
    p1 = np.dot(rotmat(udir,theta) * a.pos)
    ans.append(Atom(a.name,p1))
  return ans

#print a list of Atoms
def printatoms(filename,aa):
  fout = open(filename,"w")
  for a in aa:
    fout.write("%s %15.10f %15.10f %15.10f \n" % (a.name,a.pos[0],a.pos[1],a.pos[2]))
  fout.close()
 

