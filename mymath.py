import numpy as np

#tensor product of matrices
def tensprd(m1,m2):
  p1,q1 = m1.shape
  p2,q2 = m2.shape
  ans = np.zeros((p1*p2,q1*q2),dtype=complex)
  for i1 in range(p1):
    for j1 in range(q1):
      for i2 in range(p2):
        for j2 in range(q2):
          ans[i1*p2+i2,j1*q2+j2] = m1[i1,j1]*m2[i2,j2]
  return ans

#pauli matrices
pauli0 = np.array([[1,0],[0,1]])
paulix = np.array([[0,1],[1,0]])
pauliy = np.array([[0,-1.0j],[1.0j,0]])
pauliz = np.array([[1,0],[0,-1]])

#uniform random complex number with norm < 1
def unirand_cplx():
  while(True):
    x = np.random.random()*2 -1.0
    y = np.random.random()*2 -1.0
    z = x+1.0j*y
    if(np.abs(z)<1):
      return z
    
#random complex vector with norm=1
def randvec_cplx(n):
  while(True):
    xs = np.array([np.random.random()*2-1.0 for i in range(n)])
    ys = np.array([np.random.random()*2-1.0 for i in range(n)])
    norm2 = np.sum(xs*xs) + np.sum(ys*ys)
    if norm2<1:
      return (xs + 1.0j*ys)/np.sqrt(norm2)

     
