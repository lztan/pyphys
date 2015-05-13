#some tools for helping with quantum espresso

import numpy as np
import matplotlib.pyplot as plt
from pylab import cm
import re
import mytools
import bandtools as bt
import solidstate

class scfindata:
  """information for scf calclations. 
  from scf input file"""
  def __init__(self,filename):
    """init from scf input file"""
    f = open(filename,"r")
    fs = f.read()
    f.close()
    #crystal cell axes, units of bohr
    a = []
    ibrav_s = re.search(r"ibrav\s*=\s*(\d)",fs)
    ibrav = int(ibrav_s.group(1))
    celldm1_s = re.search(r"celldm\(\s*1\s*\)\s*=\s*([\d.]+)",fs)
    celldm2_s = re.search(r"celldm\(\s*2\s*\)\s*=\s*([\d.]+)",fs)
    celldm3_s = re.search(r"celldm\(\s*3\s*\)\s*=\s*([\d.]+)",fs)
    if ibrav==4:
      celldm1 = float(celldm1_s.group(1))
      alat = celldm1
      celldm3 = float(celldm3_s.group(1))
      a.append([celldm1,0.0,0.0])
      a.append([-1.0*celldm1*0.5,np.sqrt(3.0)*0.5*celldm1,0.0])
      a.append([0.0,0.0,alat*celldm3])
    if ibrav==2:
      celldm1 = float(celldm1_s.group(1))
      alat = celldm1
      a.append([-1.0*alat*0.5,0.0,alat*0.5])
      a.append([0.0,alat*0.5,alat*0.5])
      a.append([-1.0*alat*0.5,alat*0.5,0.0])
    if ibrav == 0:
      celldm1 = float(celldm1_s.group(1))
      alat = celldm1
      #cellcard = re.search(r"CELL_PARAMETERS.*$(\s*[\d.\-]+\s+){9}",fs,re.MULTILINE)
      cellcard = re.search(r"CELL_PARAMETERS.*$\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+",fs,re.MULTILINE)
      cellparams = [float(cellcard.group(i).strip()) for i in range(1,10)]
      cellparams1 = alat*np.array(cellparams)
      a.append(cellparams1[0:3])
      a.append(cellparams1[3:6])
      a.append(cellparams1[6:9])

    self.a = np.array(a)
    self.alat= alat
    #atoms positions, switch to angstroms
    #if specified in crystal coords, switch to angstrom
    atomcard = re.search(r"ATOMIC_POSITIONS.*K_POINTS",fs,re.DOTALL)
    atomunits_s = re.search(r"ATOMIC_POSITIONS\s+(\w+)",atomcard.group(0))
    atomunits = atomunits_s.group(1)
    atoms = re.findall(r"^\s*(\w+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)",atomcard.group(0),re.MULTILINE)
    atpos = []
    for m in atoms:
      p = [float(x) for x in m[1:]]
      if atomunits=="crystal":
        atpos.append([m[0],np.dot(np.array(p),self.a)*0.52917721092])
      if atomunits=="alat":
        atpos.append([m[0],np.array(p)*alat*0.52917721092])
    self.atpos = atpos
  

  def printwin(self,kpts):
    """print wannier90 input file.
    input: kpts = (n1,n2,n3) monkhorst-pack grid.
    """
    f = open("x.win","w")
    f.write("begin unit_cell_cart\n")
    f.write("Bohr\n")
    aa = self.a 
    for i in range(3):
      f.write("%19.27g %19.27g %19.27g\n" % (aa[i,0],aa[i,1],aa[i,2]))
    f.write("end unit_cell_cart\n\n")
    f.write("begin atoms_cart\n")
    f.write("Ang\n")
    for ap in self.atpos:
      p = ap[1] 
      f.write("%s %19.27g %19.27g %19.27g\n" % (ap[0],p[0],p[1],p[2]))
    f.write("end atoms_cart\n\n")
    f.write("mp_grid : %s %s %s\n\n" % kpts)
    f.write("begin kpoints\n")
    k1s = np.linspace(0.0,1.0,kpts[0],endpoint=False)
    k2s = np.linspace(0.0,1.0,kpts[1],endpoint=False)
    k3s = np.linspace(0.0,1.0,kpts[2],endpoint=False)
    for k1 in k1s:
      for k2 in k2s:
        for k3 in k3s:
          f.write("%s %s %s\n" % (k1,k2,k3))
    f.write("end kpoints\n\n")





class scfdata:
  """information for scf calclations.
  WARNING: There may be round-off errors because this info is 
  read from scf output file."""
  def __init__(self,filename):
    """init from scf output file"""
    f = open(filename,"r")
    fs = f.read()
    f.close()
    #crystal cell axes, units of alat
    a = []
    m = re.search(r"a\(1\)\s*=\s*\(([+\s\d.\-]+)\)",fs)
    a1s = m.group(1).split()
    a1 = [float(x) for x in a1s]
    a.append(a1)
    m = re.search(r"a\(2\)\s*=\s*\(([+\s\d.\-]+)\)",fs)
    a2s = m.group(1).split()
    a2 = [float(x) for x in a2s]
    a.append(a2)
    m = re.search(r"a\(3\)\s*=\s*\(([+\s\d.\-]+)\)",fs)
    a3s = m.group(1).split()
    a3 = [float(x) for x in a3s]
    a.append(a3)
    self.a = np.transpose(np.array(a))
    #reciprocal axes, units of 2pi/alat
    b = []
    m = re.search(r"b\(1\)\s*=\s*\(([+\s\d.\-]+)\)",fs)
    b1s = m.group(1).split()
    b1 = [float(x) for x in b1s]
    b.append(b1)
    m = re.search(r"b\(2\)\s*=\s*\(([+\s\d.\-]+)\)",fs)
    b2s = m.group(1).split()
    b2 = [float(x) for x in b2s]
    b.append(b2)
    m = re.search(r"b\(3\)\s*=\s*\(([+\s\d.\-]+)\)",fs)
    b3s = m.group(1).split()
    b3 = [float(x) for x in b3s]
    b.append(b3)
    self.b = np.transpose(np.array(b))
    #kpts, units of 2pi/alat
    matches = re.findall(r"k\([\s\d]+\)\s*=\s*\(([+\s\d.\-]+)\),\s*wk\s*=",fs)
    kcart = []
    for m in matches:
      ks = m.split()
      k = [float(x) for x in ks]
      kcart.append(k)
    self.kcart = np.array(kcart)
    #kpts, crystal coords
    kcrys = []
    ginv = np.linalg.inv(self.b)
    for k in self.kcart:
      kcrys.append(np.dot(ginv,k))
    self.kcrys = np.array(kcrys)
    #alat (bohr)
    m = re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)\s*a\.u\.",fs)
    self.alat = float(m.group(1))
    #atoms positions, alat
    matches = re.findall(r"\s*\d+\s*(\w+)\s+tau\(\s*\d+\)\s*=\s*\(([\s\d.-]+)\)",fs)
    atpos = []
    for m in matches:
      pstr = m[1].split()
      p = [float(x) for x in pstr]
      atpos.append([m[0],np.array(p)])
    self.atpos = atpos
  

  def printwin(self,kpts):
    """print wannier90 input file.
    input: kpts = (n1,n2,n3) monkhorst-pack grid.
    WARNING: There may be round-off errors because this info is 
    read from scf output file."""
    f = open("x.win","w")
    f.write("begin unit_cell_cart\n")
    f.write("Bohr\n")
    aa = self.a * self.alat
    for i in range(3):
      f.write("%s %s %s\n" % (aa[0,i],aa[1,i],aa[2,i]))
    f.write("end unit_cell_cart\n\n")
    f.write("begin atoms_cart\n")
    f.write("Bohr\n")
    for ap in self.atpos:
      p = ap[1] * self.alat
      f.write("%s %s %s %s\n" % (ap[0],p[0],p[1],p[2]))
    f.write("end atoms_cart\n\n")
    f.write("mp_grid : %s %s %s\n\n" % kpts)
    f.write("begin kpoints\n")
    k1s = np.linspace(0.0,1.0,kpts[0],endpoint=False)
    k2s = np.linspace(0.0,1.0,kpts[1],endpoint=False)
    k3s = np.linspace(0.0,1.0,kpts[2],endpoint=False)
    for k1 in k1s:
      for k2 in k2s:
        for k3 in k3s:
          f.write("%s %s %s\n" % (k1,k2,k3))
    f.write("end kpoints\n\n")



def molspectra(filename,emin,emax):
  """takes an scf output file and creates a plot of
  the energy eigenstates. Gamma pt only. """
 
  # get data
  f = open(filename,"r")
  m = re.search(r"k = .+bands \(ev\):([\s\d.\-]+)the Fermi energy is([\s\d.\-]+)ev",f.read())
  f.close()
  bands_str = m.group(1).split() 
  ef_str = m.group(2).split() 
  bands = [float(x) for x in bands_str]
  ef = float(ef_str[0])

  # make figure
  plt.figure()
  plt.subplot(1,1,1,aspect=2.0)
  x = np.linspace(emin,emax,2,endpoint=True)
  for b in bands:
    y = np.linspace(b,b,2,endpoint=True)
    plt.plot(x,y,color="blue",linestyle="-",linewidth=2.0)
  y = np.linspace(ef,ef,2,endpoint=True)
  plt.plot(x,y,color="black",linestyle="--")
  plt.xlim(emin,emax)
  plt.ylim(emin,emax)
  plt.xticks([])
  plt.ylabel("E (eV)",size="x-large")
  ax = plt.gca()
  for lab in ax.get_yticklabels():
    lab.set_fontsize(18)
  plt.savefig("molspectra.png")
  plt.show()

#def bandplot(filename,emin,emax):
#  """takes scf file and makes a bandstructure plot"""

class spaghetti:
  """band structures"""
  def __init__(self,filename):
    """reads nscf file"""
    f = open(filename,"r")
    fstr= f.read()
    matches = re.findall(r"k =([\s\d.\-]+)\([\s\w]+\)\s+bands \(ev\):([\s\d.\-]+)",fstr)
    #matches = re.findall(r"k =([\s\d.\-]+)band energies \(ev\):([\s\d.\-]+)",fstr)
    efmatch = re.search(r"the Fermi energy is([\s\d.\-]+)ev",fstr)
    ef1match = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+([\d.\-]+)\s+([\d.\-]+)",fstr)
    kpts = []
    bands = []
    for m in matches:
      ks0 = m[0].strip()
      ks1 = ks0.replace("-"," -")
      ks2 = ks1.split()
      kpts.append(map(float,ks2))
      b0 = m[1].strip()
      b1 = b0.replace("-"," -")
      b2 = b1.split()
      bands.append(list(map(float,b2)))
    self.kpts = np.array(kpts)
    self.bands = np.transpose(np.array(bands))
    if efmatch:
      self.ef = float(efmatch.group(1).strip())
    elif ef1match:
      self.ef = (float(ef1match.group(1)) + float(ef1match.group(2)) )/2.0

  def plot(self,emin,emax,sympts,symlabels):
    """creates a plot of the bandstructure. 
    sympts are the indices [n] of high-symmetry points,
    symlabels are their ["labels"] ."""
    plt.figure()
    plt.subplot(1,1,1)
    x = range(1,len(self.kpts)+1)
    for b in self.bands:
      plt.plot(x,b,color='blue',linestyle='-',linewidth=2.0)
    try:
      y = np.linspace(self.ef,self.ef,len(self.kpts),endpoint=True)
      plt.plot(x,y,color='black',linestyle='--')
    except AttributeError:
      z=1+1
    plt.ylim(emin,emax)
    plt.xticks([])
    plt.ylabel("E (eV)",size="x-large")
    for s in sympts:
      plt.vlines(s,emin,emax,color='black')
    ax = plt.gca()
    ax.set_xticks(sympts)
    ax.set_xticklabels(symlabels)
    for lab in ax.get_yticklabels():
      lab.set_fontsize(18)
    for lab in ax.get_xticklabels():
      lab.set_fontsize(18)
    plt.gcf().subplots_adjust(left=0.2)
    plt.savefig("bands.png")
    plt.show()

#bands[iband,ik]
class spaghettiv:
  """band structures from verbose nscf file"""
  def __init__(self,filename):
    """reads verbose nscf file"""
    f = open(filename,"r")
    fstr= f.read()
    #matches = re.findall(r"k =([\s\d.\-]+)\([\s\w]+\)\s+bands \(ev\):([\s\d.\-]+)",fstr)
    matches = re.findall(r"k =([\s\d.\-]+)band energies \(ev\):([\s\d.\-]+)",fstr)
    efmatch = re.search(r"the Fermi energy is([\s\d.\-]+)ev",fstr)
    ef1match = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+([\d.\-]+)\s+([\d.\-]+)",fstr)
    kpts = []
    bands = []
    for m in matches:
      ks0 = m[0].strip()
      ks1 = ks0.replace("-"," -")
      ks2 = ks1.split()
      kpts.append(list(map(float,ks2)))
      b0 = m[1].strip()
      b1 = b0.replace("-"," -")
      b2 = b1.split()
      bands.append(list(map(float,b2)))
    self.kpts = np.array(kpts)
    self.bands = np.transpose(np.array(bands))
    if efmatch:
      self.ef = float(efmatch.group(1).strip())
    elif ef1match:
      self.ef = (float(ef1match.group(1)) + float(ef1match.group(2)) )/2.0

  def plot(self,emin,emax,sympts,symlabels):
    """creates a plot of the bandstructure. 
    sympts are the indices [n] of high-symmetry points,
    symlabels are their ["labels"] ."""
    plt.figure()
    plt.subplot(1,1,1)
    x = range(1,len(self.kpts)+1)
    for b in self.bands:
      #plt.plot(x,b,color='blue',linestyle='-',linewidth=2.0)
      plt.plot(x,b,'.',color='blue',linewidth=2.0)
    y = np.linspace(self.ef,self.ef,len(self.kpts),endpoint=True)
    plt.plot(x,y,color='black',linestyle='--')
    plt.ylim(emin,emax)
    plt.xticks([])
    plt.ylabel("E (eV)",size="x-large")
    for s in sympts:
      plt.vlines(s,emin,emax,color='black')
    ax = plt.gca()
    ax.set_xticks(sympts)
    ax.set_xticklabels(symlabels)
    for lab in ax.get_yticklabels():
      lab.set_fontsize(18)
    for lab in ax.get_xticklabels():
      lab.set_fontsize(18)
    plt.gcf().subplots_adjust(left=0.2)
    plt.savefig("bands.png")
    plt.show()



class spaghetti_scf:
  """band structures"""
  def __init__(self,filename):
    """reads scf file"""
    f = open(filename,"r")
    fstr= f.read()
    matches = re.findall(r"k =([\s\d.\-]+).*bands \(ev\):([\s\d.\-]+)",fstr)
    efmatch = re.search(r"the Fermi energy is([\s\d.\-]+)ev",fstr)
    kpts = []
    bands = []
    for m in matches:
      #ks0 = m[0].strip()
      #ks1 = ks0.split()
      #kpts.append(map(float,ks1))
      b0 = m[1].strip()
      b1 = b0.split()
      bands.append(map(float,b1))
    #self.kpts = np.array(kpts)
    self.bands = np.transpose(np.array(bands))
    self.ef = float(efmatch.group(1).strip())

  def plot(self,emin,emax,sympts,symlabels):
    """creates a plot of the bandstructure. 
    sympts are the indices [n] of high-symmetry points,
    symlabels are their ["labels"] ."""
    plt.figure()
    plt.subplot(1,1,1)
    x = range(1,len(self.kpts)+1)
    for b in self.bands:
      plt.plot(x,b,color='blue',linestyle='-',linewidth=2.0)
    y = np.linspace(self.ef,self.ef,len(self.kpts),endpoint=True)
    plt.plot(x,y,color='black',linestyle='--')
    plt.ylim(emin,emax)
    plt.xticks([])
    plt.ylabel("E (eV)",size="x-large")
    for s in sympts:
      plt.vlines(s,emin,emax,color='black')
    ax = plt.gca()
    ax.set_xticks(sympts)
    ax.set_xticklabels(symlabels)
    for lab in ax.get_yticklabels():
      lab.set_fontsize(18)
    for lab in ax.get_xticklabels():
      lab.set_fontsize(18)
    plt.gcf().subplots_adjust(left=0.2)
    plt.savefig("bands.png")
    plt.show()




class XSF:
  """xcrysden data files"""
  atoms = []

  def __init__(self,fname):
    """reads xcrysden xsf file for a 2d datagrid
    dimensions stored in n1,n2; 
    extent stored in vectors a1,a2, data stored in data.
    format is data[y_index,x_index]"""
    f = open(fname,"r")
    m = re.search(r"DATAGRID_2D_UNKNOWN\n\s+(\d+)\s+(\d+)\n(.+\n)(.+\n)(.+\n)([E\s\d.\-\+]+)END_DATAGRID_2D",f.read())
    f.close()
    self.n1 = int(m.group(1))
    self.n2 = int(m.group(2))
    self.a1 = np.array([float(d) for d in m.group(4).split()])
    self.a2 = np.array([float(d) for d in m.group(5).split()])
    dataflat = [float(d) for d in m.group(6).split()]
    self.data = np.array(mytools.partn(dataflat,self.n1))

  def makeimage(self):
    plt.figure()
    plt.subplot(1,1,1)
    plt.imshow(self.data,cmap=cm.gray)
    plt.xticks([])
    plt.yticks([])
    plt.savefig("xsfdata.png")
    plt.show()

class XSF_3D:
  """xcrysden data files"""
  atoms = []

  def __init__(self,fname):
    """reads xcrysden xsf file for a 3d datagrid
    dimensions stored in n1,n2,n3; 
    extent stored in vectors a1,a2,a3, data stored in data.
    format is data[z_index,y_index,x_index]"""
    f = open(fname,"r")
    m = re.search(r"DATAGRID_3D_UNKNOWN\n\s+(\d+)\s+(\d+)\s+(\d+)\n(.+\n)(.+\n)(.+\n)(.+\n)([E\s\d.\-\+]+)END_DATAGRID_3D",f.read())
    f.close()
    self.n1 = int(m.group(1))
    self.n2 = int(m.group(2))
    self.n3 = int(m.group(3))
    self.a1 = np.array([float(d) for d in m.group(5).split()])
    self.a2 = np.array([float(d) for d in m.group(6).split()])
    self.a3 = np.array([float(d) for d in m.group(7).split()])
    dataflat = [float(d) for d in m.group(8).split()]
    temp = mytools.partn(dataflat,self.n1)
    self.data = np.array(mytools.partn(temp,self.n2))




class matdyn:
  """dynamical matrix files produced by ph.x"""
  def __init__(self,filename):
    f = open(filename,"r")
    fstr = f.read()
    #this gets the eigen vectors
    # evecs[nu,i,j] = jth component of displacement of ith atom of mode nu
    matches = re.findall(r"omega.*\s+((?:\([\s\.\d\-]+\)\s+)+)",fstr)
    #print len(matches)
    evecs = []
    for m0 in matches:
      m0s = m0.strip().split('\n')
      m0sl = [x.strip('() \n') for x in m0s]
      temp = []
      for m in m0sl:
        m1 = [float(x) for x in m.split()]
        temp.append([m1[0],m1[2],m1[4]])
      evecs.append(temp)
    self.evecs = np.array(evecs)
      
    #this gets the phonon frequencies, convert to eV
    matches = re.findall(r"omega.*THz\D*([\d.]+)\s+\[cm-1\]",fstr)
    #matches = re.findall(r"omega.*THz\D*([\d.]+)",fstr)
    self.freqs = np.array([float(m)*0.1241/1000 for m in matches])
    
    f.close()


class relaxed:
  """takes an espresso relax outputfile and finds the final structure
  as a list of Atom objects"""
  def __init__(self,filename):
    """init from scf output file"""
    f = open(filename,"r")
    fs = f.read()
    f.close()

    #final structure
    m = re.search(r"Begin final coordinates\s+ATOMIC_POSITIONS[^\n]*\n(.*)End final coordinates",fs,re.DOTALL)
    ss = m.group(1).strip().split('\n')
    sss = [s.split() for s in ss]
    coords = [solidstate.Atom(x[0],np.array([float(x[1]),float(x[2]),float(x[3])])) for x in sss]
    self.coords = coords
    self.coords_str = m.group(1).strip()

    #final energy
    m = re.search(r"Final energy\s+=\s+([-\d.]+)\s*Ry",fs,re.DOTALL)
    en = float(m.group(1))
    #print(en)
    self.en = en


class mmn:
  """reads .mmn file which is output from wannier90.
  These are the overlap matrices.
  overlaps: matrices, dict. indexed by (i,j);
  where i, j are the indices of the kpoints.
  i,j indexed from 1.
  celldiffs: cell differences between i, j
  nkpts: # of kpoints
  nbnds: # of bands
  nnei : # of neighbors each kpoint has
  """

  def __init__(self,filename):

    f = open(filename,"r")
    f.readline()
    info = f.readline()
    infos = [int(x) for x in info.split()]
    self.nbnds = infos[0]
    self.nkpts = infos[1]
    self.nnei = infos[2]
    self.celldiffs = dict([])
    self.overlaps = dict([])

    for kidx in range(self.nkpts):
      for nidx in range(self.nnei):
        head = f.readline()
        heads = [int(x) for x in head.split()]
        cd = np.array(heads[2:])
        accu = []
        for uidx in range(self.nbnds * self.nbnds):
          us = [float(x) for x in f.readline().split()]
          u = complex(*us)
          accu.append(u)
        
        if(len(heads)>0):
          self.overlaps[(heads[0],heads[1])] = np.reshape(np.array(accu), (self.nbnds,self.nbnds))
          self.celldiffs[(heads[0],heads[1])] = cd
        else:
          break

    f.close()


  
class amn:
  """reads .amn file which is output from wannier90.
  These are transformation matrices between bloch states
  and starting guess in wannier90 input file.
  Stored in amat. """

  def __init__(self,filename):

    f = open(filename,"r")
    f.readline()
    info = f.readline()
    infos = [int(x) for x in info.split()]
    self.nbnds = infos[0]
    self.nkpts = infos[1]
    self.nwann = infos[2]
    self.amat = np.zeros((self.nkpts,self.nbnds,self.nwann), dtype=complex)

    for kidx in range(self.nkpts):
      for nwidx in range(self.nwann):
        for nbidx in range(self.nbnds):
          head = f.readline().split()
          headsi = [int(head[i]) for i in range(3)]
          val = float(head[3])+ 1.0j*float(head[4])
          if( headsi[0]-1 != nbidx or headsi[1]-1 != nwidx or headsi[2]-1 != kidx):
            print("reading .amn error: ", headsi, nbidx, nwidx, kidx)
          self.amat[kidx,nbidx,nwidx] = val

    f.close()

class wanneig:
  """reads .eig file which is output from wannier90.
  These are the eigenvalues for different bands and kpoints.
  Stored in eigs. """

  def __init__(self,filename):

    mat = bt.readmat(filename)
    nb = int(mat[-1,0])
    nk = int(mat[-1,1])
    eigs = mat[:,2]
    self.eigs = np.reshape(eigs,(nb,nk), order='F')


