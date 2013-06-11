import numpy as np
from numpy import linalg as LA

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
