import numpy as np

# tools for analyzing bandstructure files

#returns 3d np array with indices in the order:
#kpoint idx, band idx, idx of (kx,ky,e) 
# The canonical file structure is:
# kx1 ky1 e <--band1
# kx1 ky1 e <--band2
# ...
# kx2 ky1 e <--band1
# kx2 ky2 e <--band2
def readbands(filename):
  #first count the number of bands
  fin = open(filename,"r")
  nbands = 0
  while True:
    l = fin.readline()
    xs = l.split()
    if(nbands==0 or (prev[0]==xs[0] and prev[1]==xs[1])):
      nbands = nbands + 1
      prev = xs[0:2]
    else:
      break
  fin.close()
  # store bands in a datastructure
  fin = open(filename,"r")
  bandi = 0
  ans = []
  for l in fin:
    xs = [float(x) for x in l.split()]
    if bandi == 0:
      temp = []
    temp.append(xs)
    if bandi == nbands - 1:
      ans.append(temp)
    bandi = (bandi + 1) % nbands
  fin.close()
  return np.array(ans)

#Lorentzian function
def lor(x,c,b):
  return (1.0/np.pi)*b/((x-c)*(x-c)+b*b)

#imag to real Kramers kronig transform
# wps should be shifted relative to ws
# in order to avoid divergence
def kkim2re(ima,wps,ws):
  ima1 = np.array(ima)
  wps1 = np.array(wps)
  ans = []
  for w in ws:
    f = ima1/(w-wps1)
    ans.append(np.trapz(f,x=wps))
  return np.array(ans)/np.pi

#calculate density of states from bandstructure
# the input is the output of readbands
def make_dos(bands,emin,emax,ne,eta):
  fout = open("dos","w")
  bs = bands[:,:,2].ravel()
  es = np.linspace(emin,emax,num=ne)
  ans = []
  for e in es:
    dose = sum(lor(bs,e,eta))
    fout.write("%s %s\n" % (e,dose))
  fout.close()

#calculate density of states from bandstructure
# the input is the output of readbands
# returns dos as numpy array in the format [[e,dos(e)]]
def calc_dos(bands,emin,emax,ne,eta):
  bs = bands[:,:,2].ravel()
  es = np.linspace(emin,emax,num=ne)
  ans = []
  for e in es:
    dose = sum(lor(bs,e,eta))
    ans.append([e,dose])
  return np.array(ans)

def broaden(xs0,ys0,xmin,xmax,nx,eta):
  """lorenztian broaden an x vs y, returns list of points as numpy array"""
  xs = np.linspace(xmin,xmax,num=nx)
  ans = []
  for x in xs:
    lorfac = lor(xs0,x,eta)
    ans.append([x,np.dot(lorfac,ys0)])
  return np.array(ans)


# read a general matrix from file
def readmat(filename):
  fin = open(filename,"r")
  ans = []
  for l in fin:
    xs = [float(x) for x in l.split()]
    ans.append(xs)
  fin.close()
  return np.array(ans)

# read a general vector from file
def readvec(filename):
  fin = open(filename,"r")
  ans = []
  for l in fin:
    ans.append(float(l.strip()))
  fin.close()
  return np.array(ans)

#print a 2d array to file
def printmat(filename,mat):
  fout = open(filename,"w")
  for r in mat:
    for x in r:
      fout.write("%s " % x)
    fout.write("\n")
  fout.close()

# print a vector to file
def printvec(filename,v):
  fout = open(filename,"w")
  for x in v:
    fout.write("%s\n" % x)
  fout.close()
