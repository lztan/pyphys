import numpy as np

# tools for analyzing bandstructure files

#returns 3d np array with indices in the order:
#kpoint idx, band idx, idx of (kx,ky,e) 
#Works for arbitrary dimenion of k
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
    ndim = len(xs)-1
    if(nbands==0 or prev==xs[0:ndim]):
      nbands = nbands + 1
      prev = xs[0:ndim]
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

#theta function
def thetafn(x):
  if x>0:
    return 1.0
  elif x<0:
    return 0.0
  else:
    return 0.5

#fermi dirac function
def fermidirac(x):
  return 1.0/(1.0+np.exp(x))

def boseeinstein(x):
  return 1.0/(np.exp(x)-1)

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

#print a tensor of arbitary rank to file
def printtens_h(fileh,tens):
  rank = len(tens.shape)
  if rank==2:
    for r in tens:
      for x in r:
        fileh.write("%s " % x)
      fileh.write("\n")
  else:
    for r in tens:
      printtens_h(fileh,r)
      fileh.write("\n")

def printtens(filename,tens):
  fout = open(filename,"w")
  printtens_h(fout,tens)
  fout.close()

#read a tensor of arbitrary rank to file
def readtens(filename):
  rank = 1 #current rank of ans
  fin = open(filename,"r")
  ans = []
  ans2 = [] #this is a matrix
  c = 0 # counts number of consec blank lines
  for l in fin:
    if len(l.strip())>0:
      xs = [float(x) for x in l.split()]
      ans2.append(xs)
      if c>1:
        empt = []
        for i in range(c-2):
          empt = [empt]
        currans = ans
        for i in range(rank-c):
          currans = currans[-1]
        currans.append(empt)
        #print "checkpt 1"
        #print ans
      c=0
    else:
      c+=1
      if c==1:
        currans = ans
        for i in range(rank-1):
          currans = currans[-1]
        currans.append(ans2)
        ans2 = []
        #print "checkpt 2"
        #print ans
      else:
        if c>rank:
          rank += 1
          empt = []
          for i in range(rank-2):
            empt = [empt]
          ans = [ans,empt]
          c=0
          #print "checkpt 3"
          #print ans

  fin.close()
  return np.array(ans)
  #return ans
    

    
